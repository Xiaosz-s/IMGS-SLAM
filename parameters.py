import os
from pathlib import Path
import yaml

_config = None


def configure(path):
    global _config
    path = Path(path).resolve()
    os.environ["SLAM_CONFIG"] = str(path)
    with path.open() as stream:
        _config = yaml.safe_load(stream)
    return _config


def get_config():
    if _config is None:
        configure(
            os.environ.get("SLAM_CONFIG", Path(__file__).with_name("config.yaml"))
        )
    return _config


def module_settings(name):
    return get_config()["parameters"][name]
