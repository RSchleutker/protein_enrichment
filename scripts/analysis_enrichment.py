import numpy as np
import pandas as pd

from protein_enrichment.python.cell import Cell
from protein_enrichment.helpers import str2dict

from skimage.io import imsave
from tifffile import TiffFile
from pathlib import Path


# GFP and YFP/Venus have different emission wavelengths and accordingly the
# optimal resolution calculated by the Leica LAS X software differs for both
# proteins. Accordingly, different values have to be used for fusions with
# either GFP or YFP/Venus.
BUFFER_GFP = np.zeros((768, 768), dtype="uint8")
BUFFER_YFP = np.zeros((728, 728), dtype="uint8")


if __name__ == "__main__":

    for embryo in Path("data", "enrichment").iterdir():
        if not embryo.is_dir() or embryo.name.startswith("_"):
            continue
        # elif embryo.joinpath("_measurements.csv").exists():
        #     continue

        print(f"Analyzing: {embryo.name}")

        data: list[dict] = []  # Will be used to construct a DataFrame.

        # Read metadata from the folder name.
        meta = str2dict(embryo.name)

        # Image will not be read completely. Instead, a handler is used and
        # only the z-slice of interest is read into memory to save some
        # time.
        with TiffFile(Path(embryo, embryo.name + ".tif")) as image:
            for file in (x for x in embryo.iterdir() if x.suffix == ".csv"):
                if file.name.startswith("_"):
                    continue

                if meta["Prot"] == "Gli":
                    cell = Cell.from_csv(file, image, BUFFER_YFP)
                else:
                    cell = Cell.from_csv(file, image, BUFFER_GFP)

                try:
                    cell.find_anchors(iterations=2)
                    intensities = cell.measure()
                except Exception as error:
                    print(error)
                    continue
                else:
                    data.append(meta | str2dict(file.stem) | intensities)

                    # Export images only if not existing already.
                    if not (path := file.with_suffix(".tif")).exists():
                        imsave(path, cell.image)
                    if not (path := file.with_suffix(".mask.tif")).exists():
                        imsave(path, cell.show_mask())

        if data:
            pd.DataFrame.\
                from_records(data).\
                to_csv(embryo.joinpath("_measurements.csv"), index=False)
