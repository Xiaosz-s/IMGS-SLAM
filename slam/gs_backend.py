from parameters import module_settings

settings = module_settings("gs_backend")
from pathlib import Path
from parameters import get_config
import random
import numpy as np
import torch
import torch.multiprocessing as mp
from munch import munchify
from lietorch import SE3, SO3
from gaussian.utils.loss_utils import l1_loss, ssim
from gaussian.scene.gaussian_model import GaussianModel
from gaussian.utils.graphics_utils import getProjectionMatrix2
from gaussian.utils.slam_utils import (
    update_pose,
    to_se3_vec,
    get_loss_normal,
    get_loss_mapping_rgbd,
    get_loss_mapping_rgb,
)
from gaussian.utils.camera_utils import Camera
from gaussian.renderer import render_fastgs
from gaussian.utils.fastgs_utils import sampling_cameras, compute_gaussian_score_fastgs
from quadtree import generate_quadtree_depth_mask


class GSBackEnd(mp.Process):

    def __init__(self, config, save_dir):
        super().__init__()
        self.config = config
        self.iteration_count = 0
        self.viewpoints = {}
        self.current_window = []
        self.initialized = False
        self.save_dir = save_dir
        self.non_kf_cache = []
        self.mapping_candidate_history = []
        self.mapping_candidate_tstamps = set()
        self.mapping_frame_tstamps = set()
        self.last_mapping_tstamp = None
        self.opt_params = munchify(config["opt_params"])
        self.gaussians = GaussianModel(
            sh_degree=settings["sh_degree"], config=self.config
        )
        self.gaussians.init_lr(settings["spatial_lr_scale"])
        self.gaussians.training_setup(self.opt_params)
        self.background = torch.tensor(
            settings["background"], dtype=torch.float32, device="cuda"
        )
        self.cameras_extent = settings["GSBackEnd.__init__.cameras_extent"]
        self.set_hyperparams()

    def set_hyperparams(self):
        self.init_itr_num = self.config["Training"]["init_itr_num"]
        self.init_gaussian_update = self.config["Training"]["init_gaussian_update"]
        self.init_gaussian_reset = self.config["Training"]["init_gaussian_reset"]
        self.init_gaussian_th = self.config["Training"]["init_gaussian_th"]
        self.init_gaussian_extent = (
            self.cameras_extent * self.config["Training"]["init_gaussian_extent"]
        )
        self.gaussian_update_every = self.config["Training"]["gaussian_update_every"]
        self.gaussian_update_offset = self.config["Training"]["gaussian_update_offset"]
        self.gaussian_th = self.config["Training"]["gaussian_th"]
        self.gaussian_extent = (
            self.cameras_extent * self.config["Training"]["gaussian_extent"]
        )
        self.gaussian_reset = self.config["Training"]["gaussian_reset"]
        self.size_threshold = self.config["Training"]["size_threshold"]
        self.window_size = self.config["Training"]["window_size"]
        self.lambda_dnormal = self.config["Training"]["lambda_dnormal"]
        self.track_keyframes_add_gaussians = self.config["Training"][
            "tracking_keyframes_add_gaussians"
        ]
        self.mapping_min_interval = self.config["Training"]["mapping_min_interval"]
        self.mapping_diff_lambda = self.config["Training"]["mapping_diff_lambda"]
        self.mapping_score_epsilon = self.config["Training"]["mapping_score_epsilon"]
        self.covis_temporal_window = self.config["Training"]["covis_temporal_window"]
        self.covis_topk = self.config["Training"]["covis_topk"]

    def _ensure_projection(self, packet):
        if hasattr(self, "projection_matrix"):
            return
        (H, W) = packet["images"].shape[-2:]
        self.K = K = list(packet["intrinsics"][0]) + [W, H]
        self.projection_matrix = (
            getProjectionMatrix2(
                znear=settings[
                    "GSBackEnd._ensure_projection.getProjectionMatrix2.znear"
                ],
                zfar=settings["GSBackEnd._ensure_projection.getProjectionMatrix2.zfar"],
                fx=K[0],
                fy=K[1],
                cx=K[2],
                cy=K[3],
                W=W,
                H=H,
            )
            .transpose(0, 1)
            .cuda()
        )

    def _build_viewpoint_from_packet(
        self,
        packet,
        item_idx=settings["GSBackEnd._build_viewpoint_from_packet.item_idx"],
    ):
        self._ensure_projection(packet)
        tstamp = int(packet["tstamp"][item_idx].item())
        pose = SE3(packet["poses"][item_idx : item_idx + 1]).matrix().cuda()[0]
        viewpoint = Camera.init_from_tracking(
            packet["images"][item_idx] / 255.0,
            packet["depths"][item_idx],
            packet["normals"][item_idx],
            pose,
            tstamp,
            self.projection_matrix,
            self.K,
            tstamp,
        )
        return (tstamp, viewpoint)

    def _apply_pose_updates_to_tstamps(self, tstamps, pose_updates, scale_updates=None):
        if tstamps.numel() == 0 or self.gaussians.unique_kfIDs.numel() == 0:
            return
        with torch.no_grad():
            match_matrix = tstamps.unsqueeze(
                1
            ) == self.gaussians.unique_kfIDs.unsqueeze(0)
            match_indices = match_matrix.nonzero()
            if match_indices.numel() == 0:
                return
            tstamp_idx = match_indices[:, 0]
            gaussian_idx = match_indices[:, 1]
            updates = pose_updates.cuda()[tstamp_idx]
            if scale_updates is None:
                updates_scale = torch.ones((updates.shape[0], 1), device="cuda")
            else:
                updates_scale = scale_updates.cuda()[tstamp_idx]
            xyz = self.gaussians.get_xyz[gaussian_idx]
            xyz = updates * xyz / updates_scale
            self.gaussians._xyz[gaussian_idx] = xyz
            scale = self.gaussians.get_scaling[gaussian_idx]
            scale = scale / updates_scale
            self.gaussians._scaling[gaussian_idx] = (
                self.gaussians.scaling_inverse_activation(scale)
            )
            rot = SO3(self.gaussians.get_rotation[gaussian_idx])
            rot = SO3(updates.data[:, 3:]) * rot
            self.gaussians._rotation[gaussian_idx] = rot.data

    def _register_mapping_frame(self, tstamp, viewpoint, depth_map):
        if tstamp in self.mapping_frame_tstamps:
            return False
        if tstamp not in self.viewpoints and all(
            (view.tstamp != tstamp for view in self.non_kf_cache)
        ):
            self.non_kf_cache.append(viewpoint)
        self.add_next_kf(tstamp, viewpoint, depth_map=depth_map)
        self.mapping_frame_tstamps.add(tstamp)
        self.last_mapping_tstamp = tstamp
        if tstamp in self.viewpoints:
            self.viewpoints[tstamp].detail_mask = viewpoint.detail_mask
        return True

    def _compute_mapping_score(self, viewpoint):
        render_pkg = render_fastgs(
            viewpoint, self.gaussians, self.background, get_flag=None, metric_map=None
        )
        render_depth = render_pkg["depth"]
        if render_depth.dim() == 3:
            render_depth = render_depth.squeeze(0)
        prior_depth = viewpoint.depth
        if not torch.is_tensor(prior_depth):
            prior_depth = torch.from_numpy(prior_depth)
        prior_depth = prior_depth.to(device=render_depth.device, dtype=torch.float32)
        if prior_depth.dim() == 3:
            prior_depth = prior_depth.squeeze(0)
        valid_render = render_depth > 0
        valid_prior = prior_depth > 0
        valid = valid_render & valid_prior
        e_cov = 1.0 - valid_render.float().mean()
        if valid.any():
            e_cons = (
                torch.abs(render_depth[valid] - prior_depth[valid])
                / (prior_depth[valid] + settings["score_epsilon"])
            ).mean()
        else:
            e_cons = torch.tensor(0.0, device=render_depth.device)
        score = e_cons + self.mapping_diff_lambda * e_cov
        return (float(score.item()), float(e_cons.item()), float(e_cov.item()))

    def _promote_mapping_candidate(self, candidate):
        tstamp = candidate["tstamp"]
        if tstamp in self.mapping_frame_tstamps:
            return False
        if (
            self.last_mapping_tstamp is not None
            and tstamp - self.last_mapping_tstamp < self.mapping_min_interval
        ):
            return False
        depth_map = candidate["packet"]["depths"][0]
        if torch.is_tensor(depth_map):
            depth_map = depth_map.cpu().numpy()
        promoted = self._register_mapping_frame(
            tstamp, candidate["viewpoint"], depth_map
        )
        if promoted:
            pass
        return promoted

    def flush_mapping_candidates(self):
        if not self.mapping_candidate_history:
            return
        remaining = sorted(
            self.mapping_candidate_history,
            key=lambda item: (-item["score"], item["tstamp"]),
        )
        for candidate in remaining:
            self._promote_mapping_candidate(candidate)
        self.mapping_candidate_history = []

    def process_mapping_candidate(self, packet):
        if not self.initialized:
            return None
        (tstamp, viewpoint) = self._build_viewpoint_from_packet(packet)
        if (
            tstamp in self.mapping_candidate_tstamps
            or tstamp in self.mapping_frame_tstamps
        ):
            return None
        (score, e_cons, e_cov) = self._compute_mapping_score(viewpoint)
        candidate = {
            "tstamp": tstamp,
            "packet": packet,
            "viewpoint": viewpoint,
            "score": score,
            "e_cons": e_cons,
            "e_cov": e_cov,
        }
        self.mapping_candidate_tstamps.add(tstamp)
        self.mapping_candidate_history.append(candidate)
        if len(self.mapping_candidate_history) >= 3:
            (prev_prev, prev, curr) = self.mapping_candidate_history[-3:]
            if (
                prev["score"] >= prev_prev["score"] - self.mapping_score_epsilon
                and prev["score"] >= curr["score"] - self.mapping_score_epsilon
            ):
                self._promote_mapping_candidate(prev)
            self.mapping_candidate_history = self.mapping_candidate_history[-2:]
        return score

    def process_track_data(
        self, packet, is_kf=settings["GSBackEnd.process_track_data.is_kf"]
    ):
        self._ensure_projection(packet)
        w2c = SE3(packet["poses"]).matrix().cuda()
        for i, idx in enumerate(packet["viz_idx"]):
            idx = idx.item()
            idx = packet["tstamp"][i].item()
            tstamp = packet["tstamp"][i].item()
            viewpoint = Camera.init_from_tracking(
                packet["images"][i] / 255.0,
                packet["depths"][i],
                packet["normals"][i],
                w2c[i],
                idx,
                self.projection_matrix,
                self.K,
                tstamp,
            )
            if idx not in self.current_window:
                if is_kf:
                    self.current_window = (
                        [idx] + self.current_window[:-1]
                        if len(self.current_window) > self.window_size
                        else [idx] + self.current_window
                    )
                if not self.initialized:
                    self.reset()
                    if is_kf:
                        self.viewpoints[idx] = viewpoint
                    else:
                        self.non_kf_cache.append(viewpoint)
                    self.add_next_kf(
                        0,
                        viewpoint,
                        depth_map=packet["depths"][0].cpu().numpy(),
                        init=True,
                    )
                    self.mapping_frame_tstamps.add(idx)
                    self.last_mapping_tstamp = idx
                    self.initialize_map(0, viewpoint)
                    self.initialized = True
                elif is_kf:
                    self.viewpoints[idx] = viewpoint
                    if (
                        self.track_keyframes_add_gaussians
                        and idx not in self.mapping_frame_tstamps
                    ):
                        self.add_next_kf(
                            idx, viewpoint, depth_map=packet["depths"][i].cpu().numpy()
                        )
                        self.mapping_frame_tstamps.add(idx)
                elif idx not in self.mapping_frame_tstamps and all(
                    (view.tstamp != tstamp for view in self.non_kf_cache)
                ):
                    self.non_kf_cache.append(viewpoint)
            elif is_kf:
                self.viewpoints[idx] = viewpoint
        if packet["pose_updates"] is not None:
            with torch.no_grad():
                self._apply_pose_updates_to_tstamps(
                    packet["tstamp"], packet["pose_updates"], packet["scale_updates"]
                )
        if is_kf:
            self.map(
                self.current_window,
                iters=settings["GSBackEnd.process_track_data.map.iters"],
            )

    def finalize(self, slam_tstamps=None):
        self.flush_mapping_candidates()
        self.color_refinement(iteration_total=self.gaussians.max_steps)
        self.gaussians.save_ply(
            str(Path(self.save_dir) / get_config()["runtime"]["map_name"])
        )
        slam_tstamps_set = (
            set(slam_tstamps.tolist()) if slam_tstamps is not None else None
        )
        poses_cw = []
        for view in self.viewpoints.values():
            if slam_tstamps_set is not None and view.tstamp not in slam_tstamps_set:
                continue
            T_w2c = np.eye(4)
            T_w2c[0:3, 0:3] = view.R.cpu().numpy()
            T_w2c[0:3, 3] = view.T.cpu().numpy()
            poses_cw.append(np.hstack(([view.tstamp], to_se3_vec(T_w2c))))
        poses_cw.sort(key=lambda x: x[0])
        return np.stack(poses_cw)

    def correct_non_kf_poses(self, traj_full):
        if len(self.non_kf_cache) == 0:
            return
        for viewpoint in self.non_kf_cache:
            tstamp = int(viewpoint.tstamp)
            old_w2c = torch.eye(4, device="cuda")
            old_w2c[:3, :3] = viewpoint.R
            old_w2c[:3, 3] = viewpoint.T
            new_w2c = traj_full[tstamp].matrix().cuda()
            old_vec = (
                torch.from_numpy(to_se3_vec(old_w2c.detach().cpu().numpy()))
                .float()
                .unsqueeze(0)
            )
            new_vec = (
                torch.from_numpy(to_se3_vec(new_w2c.detach().cpu().numpy()))
                .float()
                .unsqueeze(0)
            )
            pose_update = (SE3(new_vec.cuda()).inv() * SE3(old_vec.cuda())).to(
                device="cpu"
            )
            self._apply_pose_updates_to_tstamps(
                torch.tensor([tstamp], dtype=torch.float32),
                pose_update,
                torch.ones(1, 1),
            )
            viewpoint.update_RT(new_w2c[:3, :3], new_w2c[:3, 3])
            self.viewpoints[tstamp] = viewpoint
        self.non_kf_cache = []

    def _compute_overlap_hits(
        self,
        anchor_view,
        cand_view,
        step=settings["GSBackEnd._compute_overlap_hits.step"],
    ):
        depth = anchor_view.depth
        depth_np = depth.cpu().numpy() if torch.is_tensor(depth) else depth
        (H, W) = depth_np.shape
        (ys, xs) = np.mgrid[0:H:step, 0:W:step]
        d = depth_np[ys, xs]
        mask = d > 0
        if mask.sum() == 0:
            return (0.0, None)
        xs = xs[mask].astype(np.float32)
        ys = ys[mask].astype(np.float32)
        d = d[mask].astype(np.float32)
        (fx, fy, cx, cy) = (
            float(anchor_view.fx),
            float(anchor_view.fy),
            float(anchor_view.cx),
            float(anchor_view.cy),
        )
        Xc = (xs - cx) / fx * d
        Yc = (ys - cy) / fy * d
        Zc = d
        pts = np.stack([Xc, Yc, Zc], axis=0)
        Ra = anchor_view.R.cpu().numpy()
        Ta = anchor_view.T.cpu().numpy()
        Rc = cand_view.R.cpu().numpy()
        Tc = cand_view.T.cpu().numpy()
        Xw = Ra.T @ (pts - Ta[:, None])
        pc2 = Rc @ Xw + Tc[:, None]
        z2 = pc2[2]
        eps = settings["GSBackEnd._compute_overlap_hits.eps"]
        valid_z = z2 > eps
        if valid_z.sum() == 0:
            return (0.0, None)
        pc2 = pc2[:, valid_z]
        z2 = pc2[2]
        u = float(cand_view.fx) * pc2[0] / z2 + float(cand_view.cx)
        v = float(cand_view.fy) * pc2[1] / z2 + float(cand_view.cy)
        valid = (
            (u >= 0)
            & (u < cand_view.image_width)
            & (v >= 0)
            & (v < cand_view.image_height)
        )
        if valid.sum() == 0:
            return (0.0, None)
        u = u[valid]
        v = v[valid]
        overlap = valid.sum() / d.shape[0]
        return (float(overlap), np.stack([u, v], axis=1))

    def _sample_covisible_cameras(
        self,
        num_cams,
        neighbor_L=None,
        topk=None,
        step=settings["GSBackEnd._sample_covisible_cameras.step"],
    ):
        viewpoint_keys = sorted(self.viewpoints.keys())
        if len(viewpoint_keys) <= num_cams:
            return [self.viewpoints[idx] for idx in viewpoint_keys]
        neighbor_L = self.covis_temporal_window if neighbor_L is None else neighbor_L
        topk = self.covis_topk if topk is None else topk
        selected = []
        selected_ids = set()
        attempts = 0
        max_attempts = max(
            settings["GSBackEnd._sample_covisible_cameras.attempts_multiplier"]
            * len(viewpoint_keys),
            settings["GSBackEnd._sample_covisible_cameras.minimum_attempts"],
        )
        while len(selected) < num_cams and attempts < max_attempts:
            attempts += 1
            anchor_idx = random.choice(viewpoint_keys)
            if anchor_idx in selected_ids and len(selected) > 0:
                continue
            anchor_pos = viewpoint_keys.index(anchor_idx)
            cand_ids = viewpoint_keys[
                max(0, anchor_pos - neighbor_L) : min(
                    len(viewpoint_keys), anchor_pos + neighbor_L + 1
                )
            ]
            cand_ids = [idx for idx in cand_ids if idx != anchor_idx]
            overlaps = []
            anchor_view = self.viewpoints[anchor_idx]
            for cand_idx in cand_ids:
                (overlap, _) = self._compute_overlap_hits(
                    anchor_view, self.viewpoints[cand_idx], step=step
                )
                overlaps.append((cand_idx, overlap))
            overlaps.sort(key=lambda item: item[1], reverse=True)
            group = [anchor_idx] + [
                cand_idx for (cand_idx, overlap) in overlaps[:topk] if overlap > 0
            ]
            for cand_idx in group:
                if cand_idx in selected_ids:
                    continue
                selected_ids.add(cand_idx)
                selected.append(self.viewpoints[cand_idx])
                if len(selected) >= num_cams:
                    break
        if len(selected) < num_cams:
            for viewpoint in sampling_cameras(self.viewpoints, num_cams=num_cams):
                if viewpoint.uid in selected_ids:
                    continue
                selected_ids.add(viewpoint.uid)
                selected.append(viewpoint)
                if len(selected) >= num_cams:
                    break
        return selected

    def get_pose_finalize(self, slam_tstamps=None):
        slam_tstamps_set = (
            set(slam_tstamps.tolist()) if slam_tstamps is not None else None
        )
        poses_cw = []
        for view in self.viewpoints.values():
            if slam_tstamps_set is not None and view.tstamp not in slam_tstamps_set:
                continue
            T_w2c = np.eye(4)
            T_w2c[0:3, 0:3] = view.R.cpu().numpy()
            T_w2c[0:3, 3] = view.T.cpu().numpy()
            poses_cw.append(np.hstack(([view.tstamp], to_se3_vec(T_w2c))))
        poses_cw.sort(key=lambda x: x[0])
        return np.stack(poses_cw)

    def add_next_kf(
        self,
        frame_idx,
        viewpoint,
        init=settings["GSBackEnd.add_next_kf.init"],
        scale=settings["GSBackEnd.add_next_kf.scale"],
        depth_map=None,
    ):
        rgb_image = viewpoint.original_image
        depth_tensor = torch.from_numpy(depth_map.astype(np.float32))
        quadtree_max_area = self.config["Training"]["quadtree_max_area"]
        detail_area_thresh = self.config["Training"]["quadtree_detail_area_threshold"]
        (quadtree_mask, qtree, _) = generate_quadtree_depth_mask(
            image_tensor=rgb_image,
            depth_tensor=depth_tensor,
            threshold=settings[
                "GSBackEnd.add_next_kf.generate_quadtree_depth_mask.threshold"
            ],
            min_pixel_size=settings[
                "GSBackEnd.add_next_kf.generate_quadtree_depth_mask.min_pixel_size"
            ],
            max_area=quadtree_max_area,
            pixel_dev_threshold=None,
            use_lab=False,
            device=rgb_image.device,
            return_qtree=True,
        )
        mask_np = quadtree_mask.cpu().numpy()
        detail_mask_np = np.zeros_like(mask_np, dtype=bool)
        for child in qtree.all_children:
            if child.width * child.height <= detail_area_thresh:
                detail_mask_np[
                    child.y0 : child.y0 + child.height,
                    child.x0 : child.x0 + child.width,
                ] = True
        if not detail_mask_np.any():
            areas = np.array(
                [child.width * child.height for child in qtree.all_children],
                dtype=np.int32,
            )
            min_area = areas.min() if areas.size > 0 else 1
            for child in qtree.all_children:
                if child.width * child.height <= min_area:
                    detail_mask_np[
                        child.y0 : child.y0 + child.height,
                        child.x0 : child.x0 + child.width,
                    ] = True
        viewpoint.detail_mask = torch.from_numpy(detail_mask_np).to(
            device=rgb_image.device
        )
        depth_map_filtered = depth_map.copy()
        depth_map_filtered[~mask_np] = np.nan
        self.gaussians.extend_from_pcd_seq(
            viewpoint,
            kf_id=frame_idx,
            init=init,
            scale=scale,
            depthmap=depth_map_filtered,
        )

    def reset(self):
        self.iteration_count = 0
        self.current_window = []
        self.initialized = False
        self.mapping_candidate_history = []
        self.mapping_candidate_tstamps = set()
        self.mapping_frame_tstamps = set()
        self.last_mapping_tstamp = None
        self.gaussians.prune_points(self.gaussians.unique_kfIDs >= 0)

    def initialize_map(self, cur_frame_idx, viewpoint):
        for mapping_iteration in range(self.init_itr_num):
            self.iteration_count += 1
            render_pkg = render_fastgs(
                viewpoint,
                self.gaussians,
                self.background,
                get_flag=None,
                metric_map=None,
            )
            (
                image,
                viewspace_point_tensor,
                visibility_filter,
                radii,
                depth,
                n_touched,
            ) = (
                render_pkg["render"],
                render_pkg["viewspace_points"],
                render_pkg["visibility_filter"],
                render_pkg["radii"],
                render_pkg["depth"],
                render_pkg["n_touched"],
            )
            loss_init = get_loss_mapping_rgb(self.config, image, depth, viewpoint)
            loss_init.backward()
            with torch.no_grad():
                self.gaussians.max_radii2D[visibility_filter] = torch.max(
                    self.gaussians.max_radii2D[visibility_filter],
                    radii[visibility_filter],
                )
                self.gaussians.add_densification_stats(
                    viewspace_point_tensor, visibility_filter
                )
                if mapping_iteration % self.init_gaussian_update == 0:
                    self.gaussians.densify_and_prune(
                        self.opt_params.densify_grad_threshold,
                        self.init_gaussian_th,
                        self.init_gaussian_extent,
                        None,
                    )
                if self.iteration_count == self.init_gaussian_reset:
                    self.gaussians.reset_opacity()
                self.gaussians.optimizer.step()
                self.gaussians.optimizer.zero_grad(set_to_none=True)
        return render_pkg

    def map(self, current_window, iters, prune=settings["GSBackEnd.map.prune"]):
        if len(current_window) == 0:
            return
        viewpoint_stack = [self.viewpoints[kf_idx] for kf_idx in current_window]
        random_viewpoint_stack = []
        current_window_set = set(current_window)
        for cam_idx, viewpoint in self.viewpoints.items():
            if cam_idx not in current_window_set:
                random_viewpoint_stack.append(viewpoint)
        for _ in range(iters):
            self.iteration_count += 1
            loss_mapping = 0
            viewspace_point_tensor_acm = []
            visibility_filter_acm = []
            radii_acm = []
            viewpoints = viewpoint_stack + [
                random_viewpoint_stack[idx]
                for idx in torch.randperm(len(random_viewpoint_stack))[
                    : settings["historical_views"]
                ]
            ]
            for viewpoint in viewpoints:
                render_pkg = render_fastgs(
                    viewpoint,
                    self.gaussians,
                    self.background,
                    get_flag=None,
                    metric_map=None,
                )
                (
                    image,
                    viewspace_point_tensor,
                    visibility_filter,
                    radii,
                    depth,
                    n_touched,
                ) = (
                    render_pkg["render"],
                    render_pkg["viewspace_points"],
                    render_pkg["visibility_filter"],
                    render_pkg["radii"],
                    render_pkg["depth"],
                    render_pkg["n_touched"],
                )
                loss_mapping += (
                    self.lambda_dnormal
                    * get_loss_normal(depth, viewpoint)
                    / settings["online_normal_divisor"]
                )
                loss_mapping += get_loss_mapping_rgb(
                    self.config, image, depth, viewpoint
                )
                viewspace_point_tensor_acm.append(viewspace_point_tensor)
                visibility_filter_acm.append(visibility_filter)
                radii_acm.append(radii)
            scaling = self.gaussians.get_scaling
            isotropic_loss = torch.abs(scaling - scaling.mean(dim=1).view(-1, 1))
            loss_mapping += (
                settings["GSBackEnd.map.scale_regularization_weight"]
                * isotropic_loss.mean()
            )
            if torch.isnan(loss_mapping) or torch.isinf(loss_mapping):
                self.gaussians.optimizer.zero_grad(set_to_none=True)
                continue
            loss_mapping.backward()
            with torch.no_grad():
                for idx in range(len(viewspace_point_tensor_acm)):
                    self.gaussians.max_radii2D[visibility_filter_acm[idx]] = torch.max(
                        self.gaussians.max_radii2D[visibility_filter_acm[idx]],
                        radii_acm[idx][visibility_filter_acm[idx]],
                    )
                    self.gaussians.add_densification_stats(
                        viewspace_point_tensor_acm[idx], visibility_filter_acm[idx]
                    )
                update_gaussian = (
                    self.iteration_count % self.gaussian_update_every
                    == self.gaussian_update_offset
                )
                if update_gaussian:
                    self.gaussians.densify_and_prune(
                        self.opt_params.densify_grad_threshold,
                        self.gaussian_th,
                        self.gaussian_extent,
                        self.size_threshold,
                    )
                self.gaussian_reset = settings["GSBackEnd.map.gaussian_reset"]
                if self.iteration_count % self.gaussian_reset == 0 and (
                    not update_gaussian
                ):
                    self.gaussians.reset_opacity_nonvisible(visibility_filter_acm)
                self.gaussians.optimizer.step()
                self.gaussians.optimizer.zero_grad(set_to_none=True)

    def color_refinement(self, iteration_total):
        opt_params = []
        for view in self.viewpoints.values():
            opt_params.append(
                {
                    "params": [view.cam_rot_delta],
                    "lr": self.config["opt_params"]["pose_lr"],
                    "name": "rot_{}".format(view.uid),
                }
            )
            opt_params.append(
                {
                    "params": [view.cam_trans_delta],
                    "lr": self.config["opt_params"]["pose_lr"],
                    "name": "trans_{}".format(view.uid),
                }
            )
            if self.config["Training"]["compensate_exposure"]:
                opt_params.append(
                    {
                        "params": [view.exposure_a],
                        "lr": self.config["opt_params"]["exposure_lr"],
                        "name": "exposure_a_{}".format(view.uid),
                    }
                )
                opt_params.append(
                    {
                        "params": [view.exposure_b],
                        "lr": self.config["opt_params"]["exposure_lr"],
                        "name": "exposure_b_{}".format(view.uid),
                    }
                )
        self.keyframe_optimizers = torch.optim.Adam(opt_params)
        use_fastgs = self.config["Training"]["use_fastgs_refinement"]
        fastgs_densify_interval = self.config["Training"]["fastgs_densify_interval"]
        fastgs_densify_until = self.config["Training"]["fastgs_densify_until"]
        fastgs_prune_interval = self.config["Training"]["fastgs_prune_interval"]
        fastgs_prune_start = self.config["Training"]["fastgs_prune_start"]
        fastgs_prune_end = self.config["Training"]["fastgs_prune_end"]
        fastgs_loss_thresh = self.config["Training"]["fastgs_loss_thresh"]
        fastgs_importance_thresh = self.config["Training"]["fastgs_importance_thresh"]
        fastgs_prune_thresh = self.config["Training"]["fastgs_prune_thresh"]
        fastgs_grad_thresh = self.config["Training"]["fastgs_grad_thresh"]
        fastgs_grad_abs_thresh = self.config["Training"]["fastgs_grad_abs_thresh"]
        fastgs_dense_percent = self.config["Training"]["fastgs_dense_percent"]
        fastgs_num_sample_cams = self.config["Training"]["fastgs_num_sample_cams"]
        if use_fastgs:
            pass
        for iteration in range(1, iteration_total + 1):
            viewpoint_idx_stack = list(self.viewpoints.keys())
            viewpoint_cam_idx = viewpoint_idx_stack.pop(
                random.randint(0, len(viewpoint_idx_stack) - 1)
            )
            viewpoint_cam = self.viewpoints[viewpoint_cam_idx]
            render_pkg = render_fastgs(
                viewpoint_cam,
                self.gaussians,
                self.background,
                get_flag=None,
                metric_map=None,
            )
            (image, depth) = (render_pkg["render"], render_pkg["depth"])
            viewspace_point_tensor = render_pkg["viewspace_points"]
            visibility_filter = render_pkg["visibility_filter"]
            radii = render_pkg["radii"]
            image = (
                torch.exp(viewpoint_cam.exposure_a) * image + viewpoint_cam.exposure_b
            )
            gt_image = viewpoint_cam.original_image.cuda()
            loss = (1.0 - self.opt_params.lambda_dssim) * l1_loss(
                image, gt_image
            ) + self.opt_params.lambda_dssim * (1.0 - ssim(image, gt_image))
            loss += get_loss_mapping_rgbd(self.config, image, depth, viewpoint_cam)
            if iteration < settings["normal_decay_start"]:
                loss += self.lambda_dnormal * get_loss_normal(depth, viewpoint_cam)
            else:
                loss += (
                    self.lambda_dnormal
                    * get_loss_normal(depth, viewpoint_cam)
                    / settings["normal_decay_divisor"]
                )
            loss.backward()
            with torch.no_grad():
                if use_fastgs and iteration < fastgs_densify_until:
                    self.gaussians.max_radii2D[visibility_filter] = torch.max(
                        self.gaussians.max_radii2D[visibility_filter],
                        radii[visibility_filter],
                    )
                    self.gaussians.add_densification_stats_fastgs(
                        viewspace_point_tensor, visibility_filter
                    )
                    if iteration % fastgs_densify_interval == 0:
                        camlist = self._sample_covisible_cameras(
                            num_cams=min(fastgs_num_sample_cams, len(self.viewpoints))
                        )
                        (importance_score, pruning_score) = (
                            compute_gaussian_score_fastgs(
                                camlist,
                                self.gaussians,
                                self.background,
                                loss_thresh=fastgs_loss_thresh,
                                lambda_dssim=self.opt_params.lambda_dssim,
                                DENSIFY=True,
                            )
                        )
                        if iteration % self.gaussian_reset == 0:
                            self.gaussians.reset_opacity()
                        size_threshold = (
                            settings["screen_prune_radius"]
                            if iteration > settings["screen_prune_start"]
                            else None
                        )
                        self.gaussians.densify_and_prune_fastgs(
                            max_screen_size=size_threshold,
                            min_opacity=settings[
                                "GSBackEnd.color_refinement.gaussians.densify_and_prune_fastgs.min_opacity"
                            ],
                            extent=self.cameras_extent,
                            importance_score=importance_score,
                            pruning_score=pruning_score,
                            importance_thresh=fastgs_importance_thresh,
                            grad_thresh=fastgs_grad_thresh,
                            grad_abs_thresh=fastgs_grad_abs_thresh,
                            dense=fastgs_dense_percent,
                        )
                if (
                    use_fastgs
                    and iteration % fastgs_prune_interval == 0
                    and (fastgs_prune_start < iteration < fastgs_prune_end)
                ):
                    camlist = self._sample_covisible_cameras(
                        num_cams=min(fastgs_num_sample_cams, len(self.viewpoints))
                    )
                    (_, pruning_score) = compute_gaussian_score_fastgs(
                        camlist,
                        self.gaussians,
                        self.background,
                        loss_thresh=fastgs_loss_thresh,
                        lambda_dssim=self.opt_params.lambda_dssim,
                        DENSIFY=False,
                    )
                    self.gaussians.final_prune_fastgs(
                        min_opacity=settings[
                            "GSBackEnd.color_refinement.gaussians.final_prune_fastgs.min_opacity"
                        ],
                        pruning_score=pruning_score,
                        prune_thresh=fastgs_prune_thresh,
                    )
                self.gaussians.optimizer.step()
                self.gaussians.optimizer.zero_grad(set_to_none=True)
                lr = self.gaussians.update_learning_rate(iteration)
                self.keyframe_optimizers.step()
                self.keyframe_optimizers.zero_grad(set_to_none=True)
                update_pose(viewpoint_cam)
