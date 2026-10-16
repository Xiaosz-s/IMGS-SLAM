import argparse
import os
import re
import time
from pathlib import Path
from types import SimpleNamespace

from parameters import configure


def stream(queue, options):
    import cv2
    import numpy as np
    import torch

    files = sorted(Path(options["images"]).iterdir())[options["start"] :]
    if options["length"] is not None:
        files = files[: options["length"]]
    calibration = np.asarray(options["calibration"], dtype=float)
    matrix = np.array(
        [
            [calibration[0], 0, calibration[2]],
            [0, calibration[1], calibration[3]],
            [0, 0, 1],
        ]
    )
    for index, path in enumerate(files):
        image = cv2.cvtColor(cv2.imread(str(path)), cv2.COLOR_BGR2RGB)
        intrinsics = torch.tensor(calibration[:4])
        if options["undistort"] and len(calibration) > 4:
            image = cv2.undistort(image, matrix, calibration[4:])
        border = options["crop_border"]
        if border:
            image = image[border:-border, border:-border]
            intrinsics[2:] -= border
        height, width = image.shape[:2]
        scale = np.sqrt(options["target_area"] / (height * width))
        target_h, target_w = int(height * scale), int(width * scale)
        target_h -= target_h % options["image_multiple"]
        target_w -= target_w % options["image_multiple"]
        image = cv2.resize(image, (target_w, target_h))
        image = torch.as_tensor(image).permute(2, 0, 1)
        intrinsics[[0, 2]] *= target_w / width
        intrinsics[[1, 3]] *= target_h / height
        queue.put((index, image[None], intrinsics[None], index == len(files) - 1))
    time.sleep(options["reader_finish_delay"])


def save(system, trajectory, options):
    import numpy as np
    from lietorch import SE3

    output = Path(options["output"])
    count = system.video.counter.value
    indices = system.video.tstamp[:count].cpu().numpy().astype(int)
    poses = SE3(system.video.poses[:count]).inv().data.cpu().numpy()
    files = sorted(Path(options["images"]).iterdir())[options["start"] :]
    if options["length"] is not None:
        files = files[: options["length"]]
    timestamps = np.array(
        [float(re.findall(r"[+]?(?:\d*\.\d+|\d+)", p.name)[-1]) for p in files]
    )[:, None]
    np.save(
        output / options["intrinsics_name"],
        system.video.intrinsics[0].cpu().numpy() * 8,
    )
    np.savetxt(
        output / options["keyframe_name"],
        np.concatenate((timestamps[indices], poses), axis=1),
    )
    np.savetxt(
        output / options["trajectory_name"],
        np.concatenate((timestamps, trajectory), axis=1),
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", type=Path, default=Path(__file__).with_name("config.yaml")
    )
    args = parser.parse_args()
    config = configure(args.config)
    options = dict(config["runtime"])
    for key in [
        "images",
        "output",
        "tracking_weights",
        "depth_weights",
        "normal_weights",
    ]:
        path = Path(options[key])
        options[key] = str(
            path.resolve()
            if path.is_absolute()
            else (args.config.resolve().parent / path).resolve()
        )
    config["runtime"].update(options)
    os.environ["OMP_NUM_THREADS"] = str(options["threads"])
    os.environ["OPENBLAS_NUM_THREADS"] = str(options["threads"])
    import torch
    from torch.multiprocessing import Process, Queue
    from slam import Slam

    torch.set_num_threads(options["threads"])
    torch.cuda.set_device(options["device"])
    torch.multiprocessing.set_start_method("spawn")
    Path(options["output"]).mkdir(parents=True, exist_ok=True)
    files = sorted(Path(options["images"]).iterdir())[options["start"] :]
    if options["length"] is not None:
        files = files[: options["length"]]
    if not files:
        raise ValueError(options["images"])
    buffer = options["buffer"]
    if buffer < 0:
        buffer = min(
            options["buffer_limit"],
            len(files) // options["buffer_divisor"] + options["buffer_margin"],
        )
    queue = Queue(maxsize=options["queue_size"])
    reader = Process(target=stream, args=(queue, options))
    reader.start()
    system = None
    try:
        while True:
            index, image, intrinsics, last = queue.get()
            if system is None:
                system = Slam(
                    SimpleNamespace(
                        weights=options["tracking_weights"],
                        config=str(args.config),
                        output=options["output"],
                        buffer=buffer,
                        image_size=list(image.shape[-2:]),
                    )
                )
            system.track(index, image, intrinsics=intrinsics, is_last=last)
            if last:
                break
        reader.join()
        trajectory = system.terminate()
        save(system, trajectory, options)
    finally:
        if reader.is_alive():
            reader.terminate()
        reader.join()
        queue.close()


if __name__ == "__main__":
    main()
