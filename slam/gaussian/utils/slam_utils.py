from parameters import module_settings

settings = module_settings("gaussian.utils.slam_utils")
import math
import numpy as np
import torch
from scipy.spatial.transform import Rotation as R


def to_se3_vec(pose_mat):
    quat = R.from_matrix(pose_mat[:3, :3]).as_quat()
    return np.hstack((pose_mat[:3, 3], quat))


def skew_sym_mat(x):
    device = x.device
    dtype = x.dtype
    ssm = torch.zeros(3, 3, device=device, dtype=dtype)
    ssm[0, 1] = -x[2]
    ssm[0, 2] = x[1]
    ssm[1, 0] = x[2]
    ssm[1, 2] = -x[0]
    ssm[2, 0] = -x[1]
    ssm[2, 1] = x[0]
    return ssm


def SO3_exp(theta):
    device = theta.device
    dtype = theta.dtype
    W = skew_sym_mat(theta)
    W2 = W @ W
    angle = torch.norm(theta)
    I = torch.eye(3, device=device, dtype=dtype)
    if angle < settings["SO3_exp.small_angle_threshold"]:
        return I + W + 0.5 * W2
    else:
        return I + torch.sin(angle) / angle * W + (1 - torch.cos(angle)) / angle**2 * W2


def V(theta):
    dtype = theta.dtype
    device = theta.device
    I = torch.eye(3, device=device, dtype=dtype)
    W = skew_sym_mat(theta)
    W2 = W @ W
    angle = torch.norm(theta)
    if angle < settings["V.small_angle_threshold"]:
        V = I + 0.5 * W + 1.0 / 6.0 * W2
    else:
        V = (
            I
            + W * ((1.0 - torch.cos(angle)) / angle**2)
            + W2 * ((angle - torch.sin(angle)) / angle**3)
        )
    return V


def SE3_exp(tau):
    dtype = tau.dtype
    device = tau.device
    rho = tau[:3]
    theta = tau[3:]
    R = SO3_exp(theta)
    t = V(theta) @ rho
    T = torch.eye(4, device=device, dtype=dtype)
    T[:3, :3] = R
    T[:3, 3] = t
    return T


def update_pose(camera):
    tau = torch.cat([camera.cam_trans_delta, camera.cam_rot_delta], axis=0)
    T_w2c = torch.eye(4, device=tau.device)
    T_w2c[0:3, 0:3] = camera.R
    T_w2c[0:3, 3] = camera.T
    new_w2c = SE3_exp(tau) @ T_w2c
    new_R = new_w2c[0:3, 0:3]
    new_T = new_w2c[0:3, 3]
    camera.update_RT(new_R, new_T)
    camera.cam_rot_delta.data.fill_(0)
    camera.cam_trans_delta.data.fill_(0)


def depths_to_points(view, depthmap, world_frame):
    (W, H) = (view.image_width, view.image_height)
    fx = W / (2 * math.tan(view.FoVx / 2.0))
    fy = H / (2 * math.tan(view.FoVy / 2.0))
    intrins = (
        torch.tensor([[fx, 0.0, W / 2.0], [0.0, fy, H / 2.0], [0.0, 0.0, 1.0]])
        .float()
        .cuda()
    )
    (grid_x, grid_y) = torch.meshgrid(
        torch.arange(W, device="cuda").float() + 0.5,
        torch.arange(H, device="cuda").float() + 0.5,
        indexing="xy",
    )
    points = torch.stack([grid_x, grid_y, torch.ones_like(grid_x)], dim=-1).reshape(
        -1, 3
    )
    if world_frame:
        c2w = view.world_view_transform.T.inverse()
        rays_d = points @ intrins.inverse().T @ c2w[:3, :3].T
        rays_o = c2w[:3, 3]
        points = depthmap.reshape(-1, 1) * rays_d + rays_o
    else:
        rays_d = points @ intrins.inverse().T
        points = depthmap.reshape(-1, 1) * rays_d
    return points


def depth_to_normal(view, depth, world_frame=settings["depth_to_normal.world_frame"]):
    points = depths_to_points(view, depth, world_frame).reshape(*depth.shape[1:], 3)
    normal_map = torch.zeros_like(points)
    dx = torch.cat([points[2:, 1:-1] - points[:-2, 1:-1]], dim=0)
    dy = torch.cat([points[1:-1, 2:] - points[1:-1, :-2]], dim=1)
    normal_map[1:-1, 1:-1, :] = torch.nn.functional.normalize(
        torch.cross(dx, dy, dim=-1), dim=-1
    )
    return (normal_map, points)


def image_gradient(image):
    c = image.shape[0]
    conv_y = torch.tensor(
        [[3, 0, -3], [10, 0, -10], [3, 0, -3]], dtype=torch.float32, device="cuda"
    )
    conv_x = torch.tensor(
        [[3, 10, 3], [0, 0, 0], [-3, -10, -3]], dtype=torch.float32, device="cuda"
    )
    normalizer = 1.0 / torch.abs(conv_y).sum()
    p_img = torch.nn.functional.pad(image, (1, 1, 1, 1), mode="reflect")[None]
    img_grad_v = normalizer * torch.nn.functional.conv2d(
        p_img, conv_x.view(1, 1, 3, 3).repeat(c, 1, 1, 1), groups=c
    )
    img_grad_h = normalizer * torch.nn.functional.conv2d(
        p_img, conv_y.view(1, 1, 3, 3).repeat(c, 1, 1, 1), groups=c
    )
    return (img_grad_v[0], img_grad_h[0])


def image_gradient_mask(image, eps=settings["image_gradient_mask.eps"]):
    c = image.shape[0]
    conv_y = torch.ones((1, 1, 3, 3), dtype=torch.float32, device="cuda")
    conv_x = torch.ones((1, 1, 3, 3), dtype=torch.float32, device="cuda")
    p_img = torch.nn.functional.pad(image, (1, 1, 1, 1), mode="reflect")[None]
    p_img = torch.abs(p_img) > eps
    img_grad_v = torch.nn.functional.conv2d(
        p_img.float(), conv_x.repeat(c, 1, 1, 1), groups=c
    )
    img_grad_h = torch.nn.functional.conv2d(
        p_img.float(), conv_y.repeat(c, 1, 1, 1), groups=c
    )
    return (img_grad_v[0] == torch.sum(conv_x), img_grad_h[0] == torch.sum(conv_y))


def depth_reg(depth, gt_image, huber_eps=settings["depth_reg.huber_eps"], mask=None):
    (mask_v, mask_h) = image_gradient_mask(depth)
    (gray_grad_v, gray_grad_h) = image_gradient(gt_image.mean(dim=0, keepdim=True))
    (depth_grad_v, depth_grad_h) = image_gradient(depth)
    (gray_grad_v, gray_grad_h) = (gray_grad_v[mask_v], gray_grad_h[mask_h])
    (depth_grad_v, depth_grad_h) = (depth_grad_v[mask_v], depth_grad_h[mask_h])
    w_h = torch.exp(-10 * gray_grad_h**2)
    w_v = torch.exp(-10 * gray_grad_v**2)
    err = (w_h * torch.abs(depth_grad_h)).mean() + (
        w_v * torch.abs(depth_grad_v)
    ).mean()
    return err


def get_loss_tracking(
    config,
    image,
    depth,
    opacity,
    viewpoint,
    initialization=settings["get_loss_tracking.initialization"],
):
    if config["Training"]["monocular"]:
        return get_loss_tracking_rgb(config, image, depth, opacity, viewpoint)
    return get_loss_tracking_rgbd(config, image, depth, opacity, viewpoint)


def get_loss_tracking_rgb(config, image, depth, opacity, viewpoint):
    gt_image = viewpoint.original_image.cuda()
    (_, h, w) = gt_image.shape
    mask_shape = (1, h, w)
    rgb_boundary_threshold = config["Training"]["rgb_boundary_threshold"]
    rgb_pixel_mask = (gt_image.sum(dim=0) > rgb_boundary_threshold).view(*mask_shape)
    rgb_pixel_mask = rgb_pixel_mask * viewpoint.grad_mask
    l1 = opacity * torch.abs(image * rgb_pixel_mask - gt_image * rgb_pixel_mask)
    return l1.mean()


def get_loss_tracking_rgbd(
    config,
    image,
    depth,
    opacity,
    viewpoint,
    initialization=settings["get_loss_tracking_rgbd.initialization"],
):
    alpha = config["Training"]["alpha"]
    gt_depth = torch.from_numpy(viewpoint.depth).to(
        dtype=torch.float32, device=image.device
    )[None]
    depth_pixel_mask = (
        gt_depth > settings["get_loss_tracking_rgbd.minimum_depth"]
    ).view(*depth.shape)
    opacity_mask = (
        opacity > settings["get_loss_tracking_rgbd.opacity_threshold"]
    ).view(*depth.shape)
    l1_rgb = get_loss_tracking_rgb(config, image, depth, opacity, viewpoint)
    depth_mask = depth_pixel_mask * opacity_mask
    l1_depth = torch.abs(depth * depth_mask - gt_depth * depth_mask)
    return alpha * l1_rgb + (1 - alpha) * l1_depth.mean()


def get_loss_normal(depth_mean, viewpoint):
    prior_normal = viewpoint.normal.cuda()
    prior_normal = prior_normal.reshape(3, *depth_mean.shape[-2:]).permute(1, 2, 0)
    prior_normal_normalized = torch.nn.functional.normalize(prior_normal, dim=-1)
    (normal_mean, _) = depth_to_normal(viewpoint, depth_mean, world_frame=False)
    normal_error = 1 - (prior_normal_normalized * normal_mean).sum(dim=-1)
    valid_mask = (
        prior_normal.norm(dim=-1)
        >= settings["get_loss_normal.normal_validity_threshold"]
    ) & ~torch.isnan(normal_error)
    normal_error = torch.where(valid_mask, normal_error, torch.zeros_like(normal_error))
    if valid_mask.sum() > 0:
        return (normal_error * valid_mask).sum() / valid_mask.sum()
    else:
        return torch.tensor(0.0, device=depth_mean.device)


def get_loss_mapping_rgb(config, image, depth, viewpoint):
    gt_image = viewpoint.original_image.cuda()
    (_, h, w) = gt_image.shape
    mask_shape = (1, h, w)
    rgb_boundary_threshold = config["Training"]["rgb_boundary_threshold"]
    detail_mask = getattr(viewpoint, "detail_mask", None)
    detail_weight = config["Training"]["detail_mask_weight"]
    rgb_pixel_mask = (gt_image.sum(dim=0) > rgb_boundary_threshold).view(*mask_shape)
    weight_map = torch.ones_like(rgb_pixel_mask, dtype=image.dtype, device=image.device)
    if detail_mask is not None:
        detail_mask = detail_mask.to(image.device)
        if detail_mask.dim() == 2:
            detail_mask = detail_mask.unsqueeze(0)
        weight_map = weight_map + (detail_weight - 1.0) * detail_mask
    diff = torch.abs(image - gt_image)
    weighted_diff = diff * rgb_pixel_mask * weight_map
    denom = (rgb_pixel_mask * weight_map).sum().clamp(min=1.0)
    return weighted_diff.sum() / denom


def get_loss_mapping_rgbd(config, image, depth, viewpoint):
    alpha = config["Training"]["alpha"]
    rgb_boundary_threshold = config["Training"]["rgb_boundary_threshold"]
    detail_mask = getattr(viewpoint, "detail_mask", None)
    detail_weight = config["Training"]["detail_mask_weight"]
    gt_image = viewpoint.original_image.cuda()
    gt_depth = viewpoint.depth.to(dtype=torch.float32, device=image.device)[None]
    rgb_pixel_mask = (gt_image.sum(dim=0) > rgb_boundary_threshold).view(*depth.shape)
    depth_pixel_mask = torch.logical_and(
        gt_depth > settings["get_loss_mapping_rgbd.minimum_depth"],
        depth > settings["get_loss_mapping_rgbd.minimum_depth"],
    ).view(*depth.shape)
    weight_map = torch.ones_like(rgb_pixel_mask, dtype=image.dtype, device=image.device)
    if detail_mask is not None:
        detail_mask = detail_mask.to(image.device)
        if detail_mask.dim() == 2:
            detail_mask = detail_mask.unsqueeze(0)
        weight_map = weight_map + (detail_weight - 1.0) * detail_mask
    diff_rgb = torch.abs(image - gt_image)
    weighted_rgb = diff_rgb * rgb_pixel_mask * weight_map
    rgb_denom = (rgb_pixel_mask * weight_map).sum().clamp(min=1.0)
    l1_rgb = weighted_rgb.sum() / rgb_denom
    depth_safe = torch.clamp(
        depth, min=settings["get_loss_mapping_rgbd.inverse_depth_epsilon"]
    )
    gt_depth_safe = torch.clamp(
        gt_depth, min=settings["get_loss_mapping_rgbd.inverse_depth_epsilon"]
    )
    inv_depth_diff = torch.abs(1.0 / depth_safe - 1.0 / gt_depth_safe)
    if depth_pixel_mask.sum() > 0:
        depth_weights = depth_pixel_mask * weight_map
        l1_depth = (inv_depth_diff * depth_weights).sum() / depth_weights.sum().clamp(
            min=1.0
        )
    else:
        l1_depth = torch.tensor(0.0, device=depth.device)
    return (
        alpha * l1_rgb
        + (1 - alpha)
        * l1_depth
        * settings["get_loss_mapping_rgbd.depth_loss_multiplier"]
    )


def get_median_depth(
    depth, opacity=None, mask=None, return_std=settings["get_median_depth.return_std"]
):
    depth = depth.detach().clone()
    opacity = opacity.detach()
    valid = depth > 0
    if opacity is not None:
        valid = torch.logical_and(
            valid, opacity > settings["get_median_depth.opacity_threshold"]
        )
    if mask is not None:
        valid = torch.logical_and(valid, mask)
    valid_depth = depth[valid]
    if return_std:
        return (valid_depth.median(), valid_depth.std(), valid)
    return valid_depth.median()
