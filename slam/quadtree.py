from parameters import module_settings

settings = module_settings("quadtree")
import cv2
import numpy as np
import torch


class Node:

    def __init__(self, x0, y0, width, height):
        self.x0 = x0
        self.y0 = y0
        self.width = width
        self.height = height
        self.children = []

    def get_pixels(self, img):
        return img[self.y0 : self.y0 + self.height, self.x0 : self.x0 + self.width]

    def compute_error(self, img):
        pixels = self.get_pixels(img)
        avg_color = np.mean(pixels, axis=(0, 1))
        r_channel = pixels[:, :, 0].astype(float)
        g_channel = pixels[:, :, 1].astype(float)
        b_channel = pixels[:, :, 2].astype(float)
        r_mse = np.mean((r_channel - avg_color[0]) ** 2)
        g_mse = np.mean((g_channel - avg_color[1]) ** 2)
        b_mse = np.mean((b_channel - avg_color[2]) ** 2)
        error = (
            r_mse * settings["Node.compute_error.luminance_red"]
            + g_mse * settings["Node.compute_error.luminance_green"]
            + b_mse * settings["Node.compute_error.luminance_blue"]
        )
        normalized_error = (
            error
            * img.shape[0]
            * img.shape[1]
            / settings["Node.compute_error.error_normalization"]
        )
        return normalized_error

    def max_pixel_deviation(
        self,
        img,
        use_lab=settings["Node.max_pixel_deviation.use_lab"],
        ref=settings["Node.max_pixel_deviation.ref"],
        metric=settings["Node.max_pixel_deviation.metric"],
    ):
        region = self.get_pixels(img)
        if region.size == 0:
            return 0.0
        if use_lab:
            region_proc = cv2.cvtColor(region, cv2.COLOR_RGB2LAB)
        else:
            region_proc = region
        (h, w, c) = region_proc.shape
        flat = region_proc.reshape(-1, c).astype(np.float32)
        if ref == "median":
            ref_color = np.median(flat, axis=0)
        else:
            ref_color = np.mean(flat, axis=0)
        diffs = flat - ref_color[None, :]
        if metric == "l2":
            dist = np.linalg.norm(diffs, axis=1)
        else:
            dist = np.linalg.norm(diffs, axis=1)
        return float(dist.max())


class QTree:

    def __init__(
        self,
        img,
        threshold=settings["QTree.__init__.threshold"],
        min_pixel_size=settings["QTree.__init__.min_pixel_size"],
        pixel_dev_threshold=None,
        use_lab=settings["QTree.__init__.use_lab"],
        max_area=None,
    ):
        self.img = img
        self.threshold = threshold
        self.min_pixel_size = min_pixel_size
        self.pixel_dev_threshold = pixel_dev_threshold
        self.use_lab = use_lab
        self.max_area = max_area
        (height, width) = img.shape[:2]
        self.root = Node(0, 0, width, height)
        self.all_children = []

    def subdivide(self):
        self._recursive_subdivide(self.root)
        self.all_children = self._find_children(self.root)

    def _recursive_subdivide(self, node):
        if node.width <= self.min_pixel_size or node.height <= self.min_pixel_size:
            return
        need_split = False
        if node.compute_error(self.img) > self.threshold:
            need_split = True
        if self.pixel_dev_threshold is not None:
            max_dev = node.max_pixel_deviation(self.img, use_lab=self.use_lab)
            if max_dev > self.pixel_dev_threshold:
                need_split = True
        if self.max_area is not None:
            block_area = node.width * node.height
            if block_area > self.max_area:
                need_split = True
        if not need_split:
            return
        w1 = int(np.floor(node.width / 2.0))
        w2 = node.width - w1
        h1 = int(np.floor(node.height / 2.0))
        h2 = node.height - h1
        if w1 == 0 or w2 == 0 or h1 == 0 or (h2 == 0):
            return
        n1 = Node(node.x0, node.y0, w1, h1)
        n2 = Node(node.x0, node.y0 + h1, w1, h2)
        n3 = Node(node.x0 + w1, node.y0, w2, h1)
        n4 = Node(node.x0 + w1, node.y0 + h1, w2, h2)
        self._recursive_subdivide(n1)
        self._recursive_subdivide(n2)
        self._recursive_subdivide(n3)
        self._recursive_subdivide(n4)
        node.children = [n1, n2, n3, n4]

    def _find_children(self, node):
        if not node.children:
            return [node]
        else:
            all_children = []
            for child in node.children:
                grandchildren = self._find_children(child)
                all_children.extend(grandchildren)
            return all_children


def generate_quadtree_depth_mask(
    image_tensor,
    depth_tensor,
    threshold=settings["generate_quadtree_depth_mask.threshold"],
    min_pixel_size=settings["generate_quadtree_depth_mask.min_pixel_size"],
    device=None,
    return_qtree=settings["generate_quadtree_depth_mask.return_qtree"],
    pixel_dev_threshold=None,
    use_lab=settings["generate_quadtree_depth_mask.use_lab"],
    max_area=None,
):
    if device is None:
        device = image_tensor.device
    if len(image_tensor.shape) != 3:
        raise ValueError(
            f"期望image_tensor的shape为(C, H, W)，实际为{image_tensor.shape}"
        )
    (C, H, W) = image_tensor.shape
    if len(depth_tensor.shape) != 2:
        raise ValueError(f"期望depth_tensor的shape为(H, W)，实际为{depth_tensor.shape}")
    if depth_tensor.shape[0] != H or depth_tensor.shape[1] != W:
        raise ValueError(
            f"depth_tensor的尺寸{depth_tensor.shape}与image_tensor的尺寸({H}, {W})不匹配"
        )
    if image_tensor.is_cuda:
        img_np = image_tensor.cpu().permute(1, 2, 0).numpy()
    else:
        img_np = image_tensor.permute(1, 2, 0).numpy()
    if img_np.dtype == np.float32 or img_np.dtype == np.float64:
        if img_np.min() < 0:
            img_np = (img_np + 1.0) / 2.0
        img_np = (img_np * 255).astype(np.uint8)
    if depth_tensor.is_cuda:
        depth_np = depth_tensor.cpu().numpy()
    else:
        depth_np = depth_tensor.numpy()
    qtree = QTree(
        img_np,
        threshold=threshold,
        min_pixel_size=min_pixel_size,
        pixel_dev_threshold=pixel_dev_threshold,
        use_lab=use_lab,
        max_area=max_area,
    )
    qtree.subdivide()
    mask_np = np.zeros((H, W), dtype=bool)
    for child in qtree.all_children:
        center_x = child.x0 + child.width // 2
        center_y = child.y0 + child.height // 2
        block_depth = depth_np[
            child.y0 : child.y0 + child.height, child.x0 : child.x0 + child.width
        ]
        valid_mask = block_depth > 0
        if not np.any(valid_mask):
            continue
        (valid_y, valid_x) = np.where(valid_mask)
        global_y = valid_y + child.y0
        global_x = valid_x + child.x0
        distances = np.sqrt((global_x - center_x) ** 2 + (global_y - center_y) ** 2)
        nearest_idx = np.argmin(distances)
        nearest_y = global_y[nearest_idx]
        nearest_x = global_x[nearest_idx]
        mask_np[nearest_y, nearest_x] = True
    mask_tensor = torch.from_numpy(mask_np).to(device)
    if return_qtree:
        return (mask_tensor, qtree, img_np)
    else:
        return mask_tensor
