from parameters import module_settings

settings = module_settings("gaussian.utils.fastgs_utils")
import torch
import random
from gaussian.renderer import render_fastgs
from gaussian.utils.loss_utils import l1_loss, ssim


def sampling_cameras(viewpoint_dict, num_cams=settings["sampling_cameras.num_cams"]):
    viewpoint_list = list(viewpoint_dict.values())
    num_cams = min(num_cams, len(viewpoint_list))
    camlist = []
    indices = list(range(len(viewpoint_list)))
    for _ in range(num_cams):
        loc = random.randint(0, len(indices) - 1)
        idx = indices.pop(loc)
        camlist.append(viewpoint_list[idx])
    return camlist


def sampling_cameras_continuous(
    viewpoint_dict, num_cams=settings["sampling_cameras_continuous.num_cams"]
):
    viewpoint_list = list(viewpoint_dict.values())
    total_cams = len(viewpoint_list)
    if num_cams >= total_cams:
        return viewpoint_list
    max_start_idx = total_cams - num_cams
    start_idx = random.randint(0, max_start_idx)
    camlist = viewpoint_list[start_idx : start_idx + num_cams]
    return camlist


def get_loss(reconstructed_image, original_image):
    l1_loss = torch.mean(torch.abs(reconstructed_image - original_image), 0).detach()
    l1_loss_norm = (l1_loss - torch.min(l1_loss)) / (
        torch.max(l1_loss)
        - torch.min(l1_loss)
        + settings["get_loss.normalization_epsilon"]
    )
    return l1_loss_norm


def compute_photometric_loss(
    viewpoint_cam, image, lambda_dssim=settings["compute_photometric_loss.lambda_dssim"]
):
    gt_image = viewpoint_cam.original_image.cuda()
    Ll1 = l1_loss(image, gt_image)
    ssim_value = ssim(image, gt_image)
    loss = (1.0 - lambda_dssim) * Ll1 + lambda_dssim * (1.0 - ssim_value)
    return loss


def compute_gaussian_score_fastgs(
    camlist,
    gaussians,
    background,
    loss_thresh=settings["compute_gaussian_score_fastgs.loss_thresh"],
    lambda_dssim=settings["compute_gaussian_score_fastgs.lambda_dssim"],
    DENSIFY=settings["compute_gaussian_score_fastgs.DENSIFY"],
):
    full_metric_counts = None
    full_metric_score = None
    for view in range(len(camlist)):
        my_viewpoint_cam = camlist[view]
        render_pkg = render_fastgs(
            my_viewpoint_cam, gaussians, background, get_flag=False
        )
        render_image = render_pkg["render"]
        photometric_loss = compute_photometric_loss(
            my_viewpoint_cam, render_image, lambda_dssim
        )
        gt_image = my_viewpoint_cam.original_image.cuda()
        l1_loss_norm = get_loss(render_image, gt_image)
        metric_map = (l1_loss_norm > loss_thresh).int()
        render_pkg = render_fastgs(
            my_viewpoint_cam,
            gaussians,
            background,
            get_flag=True,
            metric_map=metric_map,
        )
        accum_loss_counts = render_pkg["accum_metric_counts"]
        if DENSIFY:
            if full_metric_counts is None:
                full_metric_counts = accum_loss_counts.clone()
            else:
                full_metric_counts += accum_loss_counts
        if full_metric_score is None:
            full_metric_score = photometric_loss * accum_loss_counts.clone()
        else:
            full_metric_score += photometric_loss * accum_loss_counts
    pruning_score = (full_metric_score - torch.min(full_metric_score)) / (
        torch.max(full_metric_score)
        - torch.min(full_metric_score)
        + settings["compute_gaussian_score_fastgs.normalization_epsilon"]
    )
    if DENSIFY:
        importance_score = torch.div(
            full_metric_counts, len(camlist), rounding_mode="floor"
        )
    else:
        importance_score = None
    return (importance_score, pruning_score)
