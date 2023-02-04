# BUILT INS
import sys
import os.path
import csv

# MODULES
from .main import MFD
from .gtif import writef


def main(
    dtm_path: str,
    mannings_path: str,
    hydrogram: list,
    lng: float,
    lat: float,
) -> None:
    """
    Runs a floods distribution modelation. 

    Parameters:
    area <str>: Name of the area
    lng <float>: Longitude of the brak point
    lat <float>: Latitude of the break point
    hydrogram <list[typle[float, float]]>: A list of pair values with time and flow representing the break hydrogram

    Returns:
    None: The script will write three raster files on the data directory
    """

    floods, drafts, speeds = None, None, None
    try:
        model = MFD(dtm_path, mannings_path, radius=3000, mute=False)
        floods, drafts, speeds = model.drainpaths((lng, lat), hydrogram)
    except KeyboardInterrupt as e:
        print(e)
        print("Keyboard Interruption")
    finally:
        if not (floods is None or drafts is None or speeds is None):
            data_dir = os.path.dirname(dtm_path)
            writef(os.path.join(data_dir, "floods.tif"), floods, dtm_path)
            writef(os.path.join(data_dir, "drafts.tif"), drafts, dtm_path)
            writef(os.path.join(data_dir, "speeds.tif"), speeds, dtm_path)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python -m mfdfloods <path:data_dir> <float:lng> <float:lat> [float:radius]")
        exit()

    kwargs = dict()
    data_dir = os.path.abspath(sys.argv[1])
    kwargs["lng"] = float(sys.argv[2])
    kwargs["lat"] = float(sys.argv[3])

    dtm_path = os.path.join(data_dir, "dtm.tif")
    if not os.path.isfile(dtm_path):
        raise FileNotFoundError(dtm_path + " does not exists")
    else:
        kwargs["dtm_path"] = dtm_path

    mannings_path = os.path.join(data_dir, "mannings.tif")
    if not os.path.isfile(mannings_path):
        raise FileNotFoundError(mannings_path + " does not exists")
    else:
        kwargs["mannings_path"] = mannings_path

    hydrogram_name = os.path.join(data_dir, "hydrogram.csv")
    if not os.path.isfile(hydrogram_name):
        raise FileNotFoundError(hydrogram_name + " does not exists")
    else:
        with open(hydrogram_name) as f:
            reader = csv.reader(f, delimiter=",", quotechar='"')
            kwargs["hydrogram"] = [row for row in reader]

    if len(sys.argv) == 5:
        kwargs["radius"] = float(sys.argv[4])

    main(**kwargs)
