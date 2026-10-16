from pathlib import Path
import yaml
from setuptools import setup
from torch.utils.cpp_extension import BuildExtension, CUDAExtension

root = Path(__file__).resolve().parents[2]
config = yaml.safe_load((root / "config.yaml").read_text())["build"]
setup(
    name="gaussian_backend",
    packages=["gaussian_backend"],
    ext_modules=[
        CUDAExtension(
            "gaussian_backend._C",
            sources=[
                "cuda_rasterizer/rasterizer_impl.cu",
                "cuda_rasterizer/forward.cu",
                "cuda_rasterizer/backward.cu",
                "rasterize_points.cu",
                "ext.cpp",
                "score.cu",
            ],
            extra_compile_args={
                "cxx": config["cxx_flags"],
                "nvcc": config["cuda_flags"]
                + [
                    "-DCUB_WRAPPED_NAMESPACE=raster_cub",
                    "-I" + str(Path(__file__).resolve().parent / "third_party/glm"),
                ],
            },
        )
    ],
    cmdclass={"build_ext": BuildExtension},
)
