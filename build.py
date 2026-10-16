import os
import subprocess
import sys
from pathlib import Path

import yaml


def main():
    root = Path(__file__).resolve().parent
    config = yaml.safe_load((root / "config.yaml").read_text())["build"]
    prefix = Path(sys.prefix)
    if config["cuda_home"] is not None:
        os.environ["CUDA_HOME"] = str(Path(config["cuda_home"]).expanduser().resolve())
    elif (prefix / "bin" / "nvcc").exists():
        os.environ["CUDA_HOME"] = str(prefix)
    os.environ["TORCH_CUDA_ARCH_LIST"] = ";".join(config["architectures"])
    os.environ["MAX_JOBS"] = str(config["jobs"])
    compiler = prefix / "bin" / config["compiler_prefix"]
    if Path(str(compiler) + "-cc").exists():
        os.environ["CC"] = str(compiler) + "-cc"
        os.environ["CXX"] = str(compiler) + "-c++"
        os.environ["CUDAHOSTCXX"] = str(compiler) + "-c++"
    for name in ["tracking", "lie", "knn", "rasterizer"]:
        subprocess.run(
            [
                sys.executable,
                "-m",
                "pip",
                "install",
                "--no-deps",
                "--no-build-isolation",
                str(root / "native" / name),
            ],
            check=True,
        )


if __name__ == "__main__":
    main()
