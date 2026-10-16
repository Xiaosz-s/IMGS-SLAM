from typing import NamedTuple
import torch.nn as nn
import torch
from . import _C


def rasterize_gaussians(
    means3D,
    means2D,
    sh,
    colors_precomp,
    opacities,
    scales,
    rotations,
    cov3Ds_precomp,
    theta,
    rho,
    raster_settings,
):
    return _RasterizeGaussians.apply(
        means3D,
        means2D,
        sh,
        colors_precomp,
        opacities,
        scales,
        rotations,
        cov3Ds_precomp,
        theta,
        rho,
        raster_settings,
    )


class _RasterizeGaussians(torch.autograd.Function):

    @staticmethod
    def forward(
        ctx,
        means3D,
        means2D,
        sh,
        colors_precomp,
        opacities,
        scales,
        rotations,
        cov3Ds_precomp,
        theta,
        rho,
        raster_settings,
    ):
        args = (
            raster_settings.bg,
            means3D,
            colors_precomp,
            opacities,
            scales,
            rotations,
            raster_settings.scale_modifier,
            cov3Ds_precomp,
            raster_settings.viewmatrix,
            raster_settings.projmatrix,
            raster_settings.tanfovx,
            raster_settings.tanfovy,
            raster_settings.image_height,
            raster_settings.image_width,
            sh,
            raster_settings.sh_degree,
            raster_settings.campos,
        )
        (
            num_rendered,
            color,
            alpha,
            depth,
            radii,
            geomBuffer,
            binningBuffer,
            imgBuffer,
            n_touched,
        ) = _C.rasterize_gaussians(*args)
        ctx.means2d_columns = means2D.shape[-1]
        ctx.raster_settings = raster_settings
        ctx.num_rendered = num_rendered
        ctx.save_for_backward(
            colors_precomp,
            means3D,
            scales,
            rotations,
            cov3Ds_precomp,
            radii,
            sh,
            geomBuffer,
            binningBuffer,
            imgBuffer,
            alpha,
        )
        if raster_settings.get_flag and raster_settings.metric_map is not None:
            metric = (
                raster_settings.metric_map.to(device=means3D.device, dtype=torch.int32)
                .reshape(raster_settings.image_height, raster_settings.image_width)
                .contiguous()
            )
            metric_counts = _C.score_gaussians(
                geomBuffer,
                binningBuffer,
                imgBuffer,
                metric,
                means3D.shape[0],
                num_rendered,
            )
        else:
            metric_counts = torch.zeros(
                means3D.shape[0], dtype=torch.int32, device=means3D.device
            )
        return (color, radii, depth, n_touched, metric_counts)

    @staticmethod
    def backward(
        ctx, grad_color, grad_radii, grad_depth, grad_n_touched, grad_metric_counts
    ):
        num_rendered = ctx.num_rendered
        raster_settings = ctx.raster_settings
        (
            colors_precomp,
            means3D,
            scales,
            rotations,
            cov3Ds_precomp,
            radii,
            sh,
            geomBuffer,
            binningBuffer,
            imgBuffer,
            alpha,
        ) = ctx.saved_tensors
        args = (
            raster_settings.bg,
            means3D,
            radii,
            colors_precomp,
            scales,
            rotations,
            raster_settings.scale_modifier,
            cov3Ds_precomp,
            raster_settings.viewmatrix,
            raster_settings.projmatrix,
            raster_settings.projmatrix_raw,
            raster_settings.tanfovx,
            raster_settings.tanfovy,
            grad_color,
            grad_depth,
            sh,
            raster_settings.sh_degree,
            raster_settings.campos,
            geomBuffer,
            num_rendered,
            binningBuffer,
            imgBuffer,
            alpha,
        )
        (
            grad_means2D,
            grad_colors_precomp,
            grad_opacities,
            grad_means3D,
            grad_cov3Ds_precomp,
            grad_sh,
            grad_scales,
            grad_rotations,
            grad_tau,
        ) = _C.rasterize_gaussians_backward(*args)
        grad_tau = torch.sum(grad_tau.view(-1, 6), dim=0)
        grad_rho = grad_tau[:3].view(1, -1)
        grad_theta = grad_tau[3:].view(1, -1)
        if ctx.means2d_columns == 3:
            grad_means2D = torch.cat(
                (grad_means2D[:, :2], grad_means2D[:, 2:4].norm(dim=-1, keepdim=True)),
                dim=-1,
            )
        grads = (
            grad_means3D,
            grad_means2D,
            grad_sh,
            grad_colors_precomp,
            grad_opacities,
            grad_scales,
            grad_rotations,
            grad_cov3Ds_precomp,
            grad_theta,
            grad_rho,
            None,
        )
        return grads


class GaussianRasterizationSettings(NamedTuple):
    image_height: int
    image_width: int
    tanfovx: float
    tanfovy: float
    bg: torch.Tensor
    scale_modifier: float
    viewmatrix: torch.Tensor
    projmatrix: torch.Tensor
    projmatrix_raw: torch.Tensor
    sh_degree: int
    campos: torch.Tensor
    get_flag: bool = None
    metric_map: torch.Tensor = None


class GaussianRasterizer(nn.Module):

    def __init__(self, raster_settings):
        super().__init__()
        self.raster_settings = raster_settings

    def forward(
        self,
        means3D,
        means2D,
        opacities,
        shs=None,
        colors_precomp=None,
        scales=None,
        rotations=None,
        cov3D_precomp=None,
        theta=None,
        rho=None,
    ):
        raster_settings = self.raster_settings
        if (
            shs is None
            and colors_precomp is None
            or (shs is not None and colors_precomp is not None)
        ):
            raise Exception(
                "Please provide excatly one of either SHs or precomputed colors!"
            )
        if (
            (scales is None or rotations is None)
            and cov3D_precomp is None
            or (
                (scales is not None or rotations is not None)
                and cov3D_precomp is not None
            )
        ):
            raise Exception(
                "Please provide exactly one of either scale/rotation pair or precomputed 3D covariance!"
            )
        if shs is None:
            shs = torch.Tensor([])
        if colors_precomp is None:
            colors_precomp = torch.Tensor([])
        if scales is None:
            scales = torch.Tensor([])
        if rotations is None:
            rotations = torch.Tensor([])
        if cov3D_precomp is None:
            cov3D_precomp = torch.Tensor([])
        if theta is None:
            theta = torch.Tensor([])
        if rho is None:
            rho = torch.Tensor([])
        return rasterize_gaussians(
            means3D,
            means2D,
            shs,
            colors_precomp,
            opacities,
            scales,
            rotations,
            cov3D_precomp,
            theta,
            rho,
            raster_settings,
        )
