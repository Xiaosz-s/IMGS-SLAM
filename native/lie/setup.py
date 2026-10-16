from pathlib import Path
import yaml
from setuptools import setup
from torch.utils.cpp_extension import BuildExtension, CUDAExtension

root = Path(__file__).resolve().parents[2]
config = yaml.safe_load((root / "config.yaml").read_text())["build"]
setup(
    name="lietorch",
    packages=["lietorch"],
    ext_modules=[
        CUDAExtension(
            "lietorch_backends",
            sources=[
                "lietorch/src/lietorch.cpp",
                "lietorch/src/lietorch_gpu.cu",
                "lietorch/src/lietorch_cpu.cpp",
            ],
            include_dirs=[
                str(Path(__file__).resolve().parent / "lietorch/include"),
                str(root / "native/eigen"),
            ],
            extra_compile_args={
                "cxx": config["cxx_flags"],
                "nvcc": config["cuda_flags"],
            },
        )
    ],
    cmdclass={"build_ext": BuildExtension},
)
