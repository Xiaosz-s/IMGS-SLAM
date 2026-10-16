#include <torch/extension.h>
#include "rasterize_points.h"

torch::Tensor ScoreGaussiansCUDA(torch::Tensor geom, torch::Tensor binning, torch::Tensor image, torch::Tensor metric, int points, int instances);

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {
  m.def("score_gaussians", &ScoreGaussiansCUDA);
  m.def("rasterize_gaussians", &RasterizeGaussiansCUDA);
  m.def("rasterize_gaussians_backward", &RasterizeGaussiansBackwardCUDA);
}