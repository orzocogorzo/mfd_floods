from typing import Optional
from numpy.typing import NDArray

class GeoTransformFit:
    def __init__(self, array: NDArray, gt: tuple, gt_ref: tuple, nodata=None):
        self._array = array

        # Cellsizes
        self.sx = gt[1]
        self.sy = gt[5]

        # Relation between cellsizes
        self.rx = gt_ref[1] / gt[1]
        self.ry = gt_ref[5] / gt[5]

        # Origin coordinate deltas
        self.dx = gt_ref[0] - gt[0]
        self.dy = gt_ref[3] - gt[3]
        self._nodata = nodata

    def proxy(self, rc: tuple) -> tuple:
        # A * R + D / S = B
        return (round(rc[0] * self.ry + self.dy / self.sy), round(rc[1] * self.rx + self.dx / self.sx))

    def __getitem__(self, rc: tuple) -> Optional[float]:
        try:
            return self._array[self.proxy(rc)]
        except ValueError:
            return self._nodata
        except IndexError:
            return self._nodata

    def __setitem__(self, rc: tuple, value: float) -> None:
        try:
            self._array[self.proxy(rc)] = value
        except ValueError:
            pass
        except IndexError:
            pass
