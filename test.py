# BUILT INS
import sys

# MODULES
from mfdfloods import MFD, gtif


def test (area, lng, lat, break_flow=120, base_flow=50, break_time=100, cellsize=5, radius=2000):
    print(area, lng, lat, break_flow, base_flow, break_time, cellsize, radius)
    try:
        model = MFD(
            dtm_path="data/%s_dtm.tif" % area,
            manning_path="data/%s_mannings.tif" % area,
            cellsize=cellsize,
            radius=radius,
            mute=False
        )
        floods, drafts, speeds = model.drainpaths((lng, lat), break_flow, base_flow, break_time)
    except KeyboardInterrupt as e:
        print(e)
        print("Keyboard Interruption")
    finally:
        gtif.writef("data/%s_floods_%s-%s.tif" % (area, lng, lat), floods, "data/%s_dtm.tif" % (area))
        gtif.writef("data/%s_drafts_%s-%s.tif" % (area, lng, lat), drafts, "data/%s_dtm.tif" % (area))
        gtif.writef("data/%s_speeds_%s-%s.tif" % (area, lng, lat), speeds, "data/%s_dtm.tif" % (area))


if __name__ == "__main__":
    kwargs = dict()
    kwargs["area"] = str(sys.argv[1])
    kwargs["lng"] = float(sys.argv[2])
    kwargs["lat"] = float(sys.argv[3])
    if len(sys.argv) >= 4: kwargs["break_flow"] = int(sys.argv[4])
    if len(sys.argv) >= 5: kwargs["base_flow"] = int(sys.argv[5])
    if len(sys.argv) >= 6: kwargs["break_time"] = int(sys.argv[6])
    if len(sys.argv) >= 7: kwargs["cellsize"] = int(sys.argv[7])
    if len(sys.argv) >= 8: kwargs["radius"] = int(sys.argv[8])
    test(**kwargs)
