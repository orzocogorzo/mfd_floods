# BUILT INS
import sys
import os.path
import csv

# MODULES
from .main import MFD
from .gtif import writef


def main(area: str, lng: float, lat: float, hydrogram: list[tuple[float, float]]) -> None:
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
        model = MFD(
            dtm_path="data/%s_dtm.tif" % area,
            manning_path="data/%s_mannings.tif" % area,
            mute=False
        )
        floods, drafts, speeds = model.drainpaths((lng, lat), hydrogram)
    except KeyboardInterrupt as e:
        print(e)
        print("Keyboard Interruption")
    finally:
        if not (floods is None or drafts is None or speeds is None):
            writef("data/%s_floods_%s-%s.tif" % (area, lng, lat), floods, "data/%s_dtm.tif" % (area))
            writef("data/%s_drafts_%s-%s.tif" % (area, lng, lat), drafts, "data/%s_dtm.tif" % (area))
            writef("data/%s_speeds_%s-%s.tif" % (area, lng, lat), speeds, "data/%s_dtm.tif" % (area))


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python -m mfd_floods <str:area_name> <float:lng> <float:lat> <path:hydrogram_csv>")
        exit()

    kwargs = dict()
    kwargs["area"] = str(sys.argv[1])
    kwargs["lng"] = float(sys.argv[2])
    kwargs["lat"] = float(sys.argv[3])
    if os.path.isfile(sys.argv[4]):
        with open(sys.argv[4], "r") as f:
            reader = csv.reader(f, delimiter=",", quotechar='"')
            kwargs["hydrogram"] = [row for row in reader]
    else:
        raise TypeError("Hydrogram file doesn't exists")

    main(**kwargs)
