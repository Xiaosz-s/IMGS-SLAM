from pathlib import Path
import yaml
from setuptools import setup
from torch.utils.cpp_extension import BuildExtension, CUDAExtension

root = Path(__file__).resolve().parents[2]
config = yaml.safe_load((root / "config.yaml").read_text())["build"]
setup(
    name="droid_backends",
    packages=[],
    ext_modules=[
        CUDAExtension(
            "droid_backends",
            sources=[
                "src/droid.cpp",
                "src/droid_kernels.cu",
                "src/correlation_kernels.cu",
                "src/altcorr_kernel.cu",
            ],
            include_dirs=[str(root / "native/eigen")],
            extra_compile_args={
                "cxx": config["cxx_flags"],
                "nvcc": config["cuda_flags"],
            },
        )
    ],
    cmdclass={"build_ext": BuildExtension},
)
