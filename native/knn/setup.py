from pathlib import Path
import yaml
from setuptools import setup
from torch.utils.cpp_extension import BuildExtension, CUDAExtension

root = Path(__file__).resolve().parents[2]
config = yaml.safe_load((root / "config.yaml").read_text())["build"]
setup(
    name="simple_knn",
    packages=["simple_knn"],
    ext_modules=[
        CUDAExtension(
            "simple_knn._C",
            sources=["spatial.cu", "simple_knn.cu", "ext.cpp"],
            extra_compile_args={
                "cxx": config["cxx_flags"],
                "nvcc": config["cuda_flags"],
            },
        )
    ],
    cmdclass={"build_ext": BuildExtension},
)
