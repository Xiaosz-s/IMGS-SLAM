from parameters import module_settings

settings = module_settings("motion_filter")
from parameters import get_config
import cv2
import torch
import numpy as np
import torch.nn.functional as F
import geom.projective_ops as pops
from modules.corr import CorrBlock
from torchvision import transforms
from midas.omnidata import OmnidataModel


def sharpness(image):
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    return cv2.Laplacian(gray, cv2.CV_64F).var()


class MotionFilter:

    def __init__(
        self, net, video, config, device=settings["MotionFilter.__init__.device"]
    ):
        self.cnet = net.cnet
        self.fnet = net.fnet
        self.update = net.update
        self.video = video
        self.thresh = config["thresh"]
        self.init_thresh = config["init_thresh"]
        self.device = device
        self.count = 0
        self.omni_dep = None
        self.deltas = [0]
        self.current_kf_tstamp = None
        self.skip_blur = config["skip_blur"]
        self.cache = [None] * settings["blur_window"]
        self.shapeness = [0] * settings["blur_window"]
        self.MEAN = torch.as_tensor(settings["image_mean"], device=self.device)[
            :, None, None
        ]
        self.STDV = torch.as_tensor(settings["image_std"], device=self.device)[
            :, None, None
        ]
        self.prev_loss_ratio = None
        self.loss_increase_streak = 0
        self.prev_frame_cache = None
        self.loss_peak_min_increase = config["loss_peak_min_increase"]
        self.loss_ratio_epsilon = config["loss_ratio_epsilon"]

    @torch.cuda.amp.autocast(enabled=True)
    def context_encoder(self, image):
        (net, inp) = self.cnet(image).split([128, 128], dim=2)
        return (net.tanh().squeeze(0), inp.relu().squeeze(0))

    @torch.cuda.amp.autocast(enabled=True)
    def feature_encoder(self, image):
        return self.fnet(image).squeeze(0)

    @torch.cuda.amp.autocast(enabled=True)
    @torch.no_grad()
    def prior_extractor(self, im_tensor):
        input_size = im_tensor.shape[-2:]
        trans_totensor = transforms.Compose(
            [transforms.Resize(tuple(settings["prior_size"]), antialias=True)]
        )
        im_tensor = trans_totensor(im_tensor).cuda()
        if self.omni_dep is None:
            self.omni_dep = OmnidataModel(
                "depth", get_config()["runtime"]["depth_weights"], device="cuda"
            )
            self.omni_normal = OmnidataModel(
                "normal", get_config()["runtime"]["normal_weights"], device="cuda"
            )
        depth = self.omni_dep(im_tensor)[None] * settings["prior_depth_scale"]
        depth = F.interpolate(depth, input_size, mode="bicubic")
        depth = depth.float().squeeze()
        normal = self.omni_normal(im_tensor) * 2.0 - 1.0
        normal = F.interpolate(normal, input_size, mode="bicubic")
        normal = normal.float().squeeze()
        return (depth, normal)

    def compute_missing_ratio(
        self,
        delta,
        ht=None,
        wd=None,
        return_mask=settings["MotionFilter.compute_missing_ratio.return_mask"],
    ):
        if delta.dim() == 5 and delta.shape[2] == 2:
            delta = delta.permute(0, 3, 4, 2)
        elif delta.dim() == 4 and delta.shape[1] == 2:
            delta = delta.permute(0, 2, 3, 1)
        elif delta.dim() == 4 and delta.shape[-1] == 2:
            pass
        else:
            raise ValueError(
                f"Unexpected delta shape {delta.shape}, expected [1, H, W, 2] or [1, 2, H, W]."
            )
        if ht is None or wd is None:
            (ht, wd) = (delta.shape[1], delta.shape[2])
        coords0 = pops.coords_grid(ht, wd, device=self.device)[None]
        coords1 = coords0 + delta
        x_in = (coords1[..., 0] > settings["flow_border"]) & (
            coords1[..., 0] < wd - settings["flow_border"]
        )
        y_in = (coords1[..., 1] > settings["flow_border"]) & (
            coords1[..., 1] < ht - settings["flow_border"]
        )
        valid_mask = (x_in & y_in).squeeze(1)
        valid_ratio = valid_mask.float().mean().item()
        overlap_ratio = valid_ratio
        missing_ratio = 1.0 - overlap_ratio
        if return_mask:
            return (missing_ratio, overlap_ratio, valid_mask)
        return (missing_ratio, overlap_ratio)

    @torch.cuda.amp.autocast(enabled=True)
    @torch.no_grad()
    def track(
        self,
        tstamp,
        image,
        intrinsics=None,
        is_last=settings["MotionFilter.track.is_last"],
    ):
        ht = image.shape[-2] // 8
        wd = image.shape[-1] // 8
        intrinsics[:, :4] /= 8.0
        s = sharpness(image[0].permute(1, 2, 0).cpu().numpy())
        inputs = image[None].to(self.device) / 255.0
        inputs = inputs.sub_(self.MEAN).div_(self.STDV)
        gmap = self.feature_encoder(inputs)
        if self.video.counter.value == 0:
            (depth, normal) = self.prior_extractor(inputs[0])
            (net, inp) = self.context_encoder(inputs[:, [0]])
            (self.net, self.inp, self.fmap) = (net, inp, gmap)
            self.video.append(
                tstamp,
                image[0],
                None,
                1.0,
                depth,
                normal,
                intrinsics,
                gmap,
                net[0],
                inp[0],
            )
            self.current_kf_tstamp = tstamp
        else:
            coords0 = pops.coords_grid(ht, wd, device=self.device)[None, None]
            corr = CorrBlock(self.fmap[None, [0]], gmap[None, [0]])(coords0)
            (_, delta, weight) = self.update(self.net[None], self.inp[None], corr)
            self.deltas.append(delta.norm(dim=-1).mean().item())
            delta_for_metric = delta[0, 0].unsqueeze(0)
            (current_loss_ratio, _, _) = self.compute_missing_ratio(
                delta_for_metric, ht, wd, return_mask=True
            )
            if self.prev_loss_ratio is not None:
                if current_loss_ratio > self.prev_loss_ratio + self.loss_ratio_epsilon:
                    self.loss_increase_streak += 1
                elif (
                    current_loss_ratio < self.prev_loss_ratio - self.loss_ratio_epsilon
                ):
                    self.loss_increase_streak = 0
            thresh = self.init_thresh if not self.video.is_initialized else self.thresh
            delta_mean = delta.norm(dim=-1).mean().item()
            reset_loss_history = False
            if delta_mean > thresh:
                index_min = np.argmax(self.shapeness)
                if self.skip_blur and self.shapeness[index_min] > s:
                    (tstamp, image, intrinsics, gmap, inputs) = self.cache[index_min]
                self.shapeness = [0] * settings["blur_window"]
                self.cache = [None] * settings["blur_window"]
                (depth, normal) = self.prior_extractor(inputs[0])
                self.count = 0
                (net, inp) = self.context_encoder(inputs[:, [0]])
                (self.net, self.inp, self.fmap) = (net, inp, gmap)
                self.video.append(
                    tstamp,
                    image[0],
                    None,
                    None,
                    depth,
                    normal,
                    intrinsics,
                    gmap,
                    net[0],
                    inp[0],
                )
                self.current_kf_tstamp = tstamp
                reset_loss_history = True
            elif (
                self.prev_loss_ratio is not None
                and current_loss_ratio < self.prev_loss_ratio - self.loss_ratio_epsilon
                and (self.loss_increase_streak >= self.loss_peak_min_increase)
                and (self.prev_frame_cache is not None)
            ):
                (prev_tstamp, prev_image, prev_intrinsics, prev_gmap, prev_inputs) = (
                    self.prev_frame_cache
                )
                (depth, normal) = self.prior_extractor(prev_inputs[0])
                (net, inp) = self.context_encoder(prev_inputs[:, [0]])
                (self.net, self.inp, self.fmap) = (net, inp, prev_gmap)
                self.video.append(
                    prev_tstamp,
                    prev_image[0],
                    None,
                    None,
                    depth,
                    normal,
                    prev_intrinsics,
                    prev_gmap,
                    net[0],
                    inp[0],
                )
                self.current_kf_tstamp = prev_tstamp
                with self.video.get_lock():
                    self.video.mapping_only[self.video.counter.value - 1] = True
                self.count = 0
                self.shapeness = [0] * settings["blur_window"]
                self.cache = [None] * settings["blur_window"]
                reset_loss_history = True
            else:
                self.shapeness[tstamp % settings["blur_window"]] = s
                self.cache[tstamp % settings["blur_window"]] = [
                    tstamp,
                    image,
                    intrinsics,
                    gmap,
                    inputs,
                ]
                self.count += 1
            self.prev_frame_cache = [tstamp, image, intrinsics.clone(), gmap, inputs]
            if reset_loss_history:
                self.prev_loss_ratio = None
                self.loss_increase_streak = 0
            else:
                self.prev_loss_ratio = current_loss_ratio
