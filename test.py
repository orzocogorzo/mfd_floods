# BUILT INS
from ctypes import ArgumentError
import sys
import os.path
import csv

# MODULES
from mfdfloods import MFD, gtif


def test (area: str, lng: float, lat: float, hydrogram: list[tuple[float, float]]):
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
        gtif.writef("data/%s_floods_%s-%s.tif" % (area, lng, lat), floods, "data/%s_dtm.tif" % (area))
        gtif.writef("data/%s_drafts_%s-%s.tif" % (area, lng, lat), drafts, "data/%s_dtm.tif" % (area))
        gtif.writef("data/%s_speeds_%s-%s.tif" % (area, lng, lat), speeds, "data/%s_dtm.tif" % (area))


if __name__ == "__main__":
    if len(sys.argv) < 5:
        raise TypeError("Invalid arguments")

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

    test(**kwargs)
