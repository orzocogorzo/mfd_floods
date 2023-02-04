from numpy.typing import NDArray

class GeoTransformFit:
    def __init__(self, array: NDArray, gt: tuple, gt_ref: tuple, nodata=None):
        self._array = array
        self.rx = gt_ref[1] / gt[1]
        self.ry = gt_ref[5] / gt[5]
        self.dx = round((gt_ref[0] - gt[0]) / gt_ref[1])
        self.dy = round((gt_ref[3] - gt[3]) / gt_ref[5])
        self._nodata = nodata

    def proxy(self, rc: tuple) -> tuple:
        return (round(rc[0] * self.ry) + self.dy, round(rc[1] * self.rx) + self.dx)

    def __getitem__(self, rc: tuple) -> None:
        try:
            return self._array[self.proxy(rc)]
        except IndexError:
            return self._nodata

    def __setitem__(self, rc: tuple, value: float) -> None:
        try:
            self._array[self.proxy(rc)] = value
        except IndexError:
            pass
