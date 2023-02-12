# import sys
# sys.setrecursionlimit(3000)

# MODULES
import numpy as np
np.seterr(all="warn")


class Matrix (object):

    def __init__ (self, dtm_array):
        self.deltas = np.array([
            (-1, -1), (-1, 0), (-1, +1),
            (0, -1), (0, +1),
            (+1, -1), (+1, 0), (+1, +1)
        ])
        
        self.dtm = self.array(dtm_array)

    def displace (self, mat, direction, delta=1):
        direction = direction.lower()
        tmp = self.zeros(mat.shape)
        if direction == "l" or direction == "left":
            tmp[:, -1*delta:] = float("nan")
            tmp[:, :-1*delta] = mat[:, delta:]
        elif direction == "r" or direction == "right":
            tmp[:, :delta] = float("nan")
            tmp[:, delta:] = mat[:, :-1*delta]
        elif direction == "b" or direction == "bottom":
            tmp[:delta, :] = float("nan")
            tmp[delta:, :] = mat[:-1*delta, :]
        elif direction == "t" or direction == "top":
            tmp[-1*delta:, :] = float("nan")
            tmp[:-1*delta, :] = mat[delta:, :]

        return tmp

    def zeros (self, shape):
        return np.zeros(shape)

    def array (self, data):
        if not isinstance(data, np.ndarray):
            if isinstance(data, list):
                return np.array(data)
            else:
                raise Exception("Invalid dtm data")
        else:
            return data

    def where (self, index, a, b):
        return np.where(index, a, b)

    def argwhere (self, index):
        return np.argwhere(index)

    def log_and (self, a, b):
        return np.logical_and(a, b)
