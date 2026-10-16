from parameters import module_settings

settings = module_settings("midas.omnidata")
import torch
from torchvision import transforms
from midas.dpt_depth import DPTDepthModel


class OmnidataModel:
    backbone = "vitb_rn50_384"
    patch_size = 32
    channel_dict = {"depth": 1, "normal": 3}

    def __init__(
        self,
        task=settings["OmnidataModel.__init__.task"],
        model_path=None,
        device=settings["OmnidataModel.__init__.device"],
    ):
        self.model_path = model_path
        self.task = task
        self.channel = self.channel_dict[task]
        self.device = device
        self.model = DPTDepthModel(backbone=self.backbone, num_channels=self.channel)
        checkpoint = torch.load(self.model_path, map_location=device)
        assert "state_dict" in checkpoint, "No state_dict found in checkpoint"
        state_dict = {}
        for k, v in checkpoint["state_dict"].items():
            state_dict[k[len("model.") :]] = v
        self.model.load_state_dict(state_dict)
        self.model.to(device)
        im2tensor = [transforms.ToTensor()]
        if task == "depth":
            im2tensor.append(transforms.Normalize(mean=0.5, std=0.5))
        self.im2tensor = transforms.Compose(im2tensor)

    def __call__(
        self,
        im_tensor,
        im_fname=None,
        down_factor=settings["OmnidataModel.__call__.down_factor"],
        target_size=None,
    ):
        with torch.no_grad():
            output = self.model(im_tensor)
        return output
