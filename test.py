import os.path
import sys
from osgeo import gdal
from mfdfloods.geotransform import GeoTransformFit
from random import randint
from numpy.typing import NDArray
from numpy import zeros
from typing import Generator


def adrift(array: NDArray, steps: int = 100) -> Generator:
    y, x = array.shape

    for i in range(steps):
        yield randint(0, y), randint(0, x)


def main(dtm_filename: str, man_filename: str, vp_filename: str) -> None:
    dtm = gdal.Open(dtm_filename)
    man = gdal.Open(man_filename)
    vp = gdal.Open(vp_filename)

    dtm_gt = dtm.GetGeoTransform()
    man_gt = man.GetGeoTransform()
    vp_gt = vp.GetGeoTransform()
    print(dtm_gt, man_gt, vp_gt)
    # import pdb; pdb.set_trace()

    dtm_arr = dtm.GetRasterBand(1).ReadAsArray()
    man_arr = man.GetRasterBand(1).ReadAsArray()
    vp_arr = vp.GetRasterBand(1).ReadAsArray()
    fit = GeoTransformFit(vp_arr, vp_gt, man_gt)
    man_res = zeros(dtm_arr.shape)
    fit_res = zeros(vp_arr.shape)

    path = adrift(man_arr)

    v, n = 0, 0
    try:
        while step := next(path):
            if fit[step] is None:
                continue

            if fit[step] == man_arr[step] or fit[step] < 0 and man_arr[step] < 0:
                v += 1
                print('[X]', fit[step], man_arr[step])
            else:
                print('[ ]', fit[step], man_arr[step])

            n += 1

            print("x: ", (dtm_gt[0] + step[1] * dtm_gt[1]) - (vp_gt[0] + fit.proxy(step)[1] * vp_gt[1]), "y: ", (dtm_gt[3] + step[0] * dtm_gt[3]) - (vp_gt[3] + fit.proxy(step)[0] * vp_gt[3]))
    except StopIteration:
        pass

    print(f"{v}/{n} = {v / n * 100}%")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("USAGE: python test.py dmt_filename mannings_filename")
        exit()

    if not os.path.isfile(sys.argv[1]):
        raise FileNotFoundError(sys.argv[1])
    if not os.path.isfile(sys.argv[2]):
        raise FileNotFoundError(sys.argv[2])
    if not os.path.isfile(sys.argv[3]):
        raise FileNotFoundError(sys.argv[3])

    main(*sys.argv[1:])
