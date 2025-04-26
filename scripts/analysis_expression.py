import numpy as np
import pandas as pd

from protein_enrichment.helpers import str2dict

from skimage.io import imread
from skimage.io import imsave
from skimage.filters import median
from skimage.filters import gaussian
from skimage.filters.thresholding import threshold_niblack
from skimage.morphology import disk
from skimage.morphology import binary_erosion
from skimage.util import img_as_ubyte
from pathlib import Path


# Relativ path to the input directory and output directory.
IN: Path = Path("../data", "expression")
OUT: Path = Path("../data", "expression", "expression.csv")

if __name__ == "__main__":
    data = []
    mask = np.zeros((12, 960, 960), dtype=bool)

    for embryo in IN.iterdir():
        if not embryo.is_dir() or embryo.name.startswith("_"):
            continue

        meta = str2dict(embryo.name)

        try:
            image = imread(Path(embryo, "image.tif"))
            slices, _, _ = image.shape
        except FileNotFoundError:
            continue
        except PermissionError:
            continue

        for z in range(slices):
            processed = median(image[z], footprint=disk(1))
            processed = gaussian(processed, sigma=1)
            processed = img_as_ubyte(processed)
            threshold = threshold_niblack(processed, window_size=51)
            np.clip(threshold, a_min=2, a_max=None)
            mask[z] = binary_erosion(processed > threshold, disk(1))

        imsave(Path(embryo, "mask.tif"), img_as_ubyte(mask))
        roi = image[mask]
        mean = roi.mean()
        med = np.percentile(roi, 50)
        p99 = np.percentile(roi, 99)

        data += [pd.DataFrame({**meta,
                               "Mean": [mean],
                               "Median": [med],
                               "P99": [p99]})]

    pd.concat(data).to_csv(OUT, index=False)
