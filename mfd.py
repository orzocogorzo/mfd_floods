# import sys
# sys.setrecursionlimit(3000)

# SYS
import math
from functools import reduce

# MODULES
import numpy as np
from osgeo import gdal
import richdem as rd


def _mat (mat, direction):
    direction = direction.lower()
    tmp = np.zeros(mat.shape)
    if direction == "l" or direction == "left":
        tmp[:, 0] = [float("inf") * mat.shape[1]]
    elif direction == "r" or direction == "right":
        tmp[:, -1] = [float("inf") * mat.shape[1]]
    elif direction == "t" or direction == "top":
        tmp[0, :] = [float("inf") * mat.shape[0]]
    elif direction == "b" or direction == "bottom":
        tmp[-1, :] = [float("inf") * mat.shape[0]]
        
    return tmp
    
    
def displace2top (mat, delta=1):
    tmp = _mat(mat, "T")
    tmp[delta:, :] = mat[:-1*delta, :]
    return tmp
    
    
def displace2bottom (mat, delta=1):
    tmp = _mat(mat, "B")
    tmp[:-1*delta, :] = mat[delta:, :]
    return tmp
    

def displace2right (mat, delta=1):
    tmp = _mat(mat, "R")
    tmp[:, :-1*delta] = mat[:, delta:]
    return tmp
    

def displace2left (mat, delta=1):
    tmp = _mat(mat, "L")
    tmp[:, delta:] = mat[:, :-1*delta]
    return tmp


def hydrogram (f):
    f = float(f)
    
    def _range (start, stop, step):
        x = start
        while x <= stop:
             yield x
             x += step
    
    time = _range(0.0, 1.0, 0.01)
    
    def _hydrogram ():
        for t in time:
            yield f**-t * f

    return _hydrogram()


class MFD (object):

    def __init__ (self, dtm_array, pxsize):
        self.deltas = [
            (-1, -1), (-1, 0), (-1, +1),
            (0, -1), (0, 0), (0, +1),
            (+1, -1), (+1, 0), (+1, +1)
        ]

        self.cardinalities = [
            "LT", "T", "RT",
            "L", None, "R",
            "LB", "B", "RB"
        ]
        
        if not isinstance(dtm_array, np.ndarray):
            if isinstance(dtm_array, list):
                self.dtm = np.array(dtm_array)
            else:
                raise Exception("Invalid dtm data")
        else:
            self.dtm = dtm_array
            
        self.dtm = rd.rdarray(self.dtm, no_data=-99)
        rd.FillDepressions(self.dtm, in_place=True)
            
        self.pxsize = pxsize
        
        self.slopes = {card: None for card in self.cardinalities if card}
        
        for course in self.slopes:
            if course == "T":
                self.slopes["T"] = displace2top(self.dtm) - self.dtm
            elif course == "L":
                self.slopes["L"] = displace2left(self.dtm) - self.dtm
            elif course == "B":
                self.slopes["B"] = displace2bottom(self.dtm) - self.dtm
            elif course == "R":
                self.slopes["R"] = displace2right(self.dtm) - self.dtm
            elif course == "RT":
                self.slopes["RT"] = displace2right(displace2top(self.dtm)) - self.dtm
            elif course == "RB":
                self.slopes["RB"] = displace2right(displace2bottom(self.dtm)) - self.dtm
            elif course == "LT":
                self.slopes["LT"] = displace2left(displace2top(self.dtm)) - self.dtm
            elif course == "LB":
                self.slopes["LB"] = displace2left(displace2bottom(self.dtm)) - self.dtm


    def get_downslopes (self, perimetter):
        return {
            self.deltas[i]: val for i, val in enumerate(perimetter)
            if val < 0
        }

    def get_perimetters (self, row, col):
        return np.array([
            self.slopes.get(cardinality)[row, col] if cardinality else float("inf")
            for cardinality in self.cardinalities
        ])
                
    def get_volumetry (self, perimetters, size=5):
        less = perimetters.min()
        return reduce(lambda a, d: a + (self.pxsize**2/2) * (self.pxsize * less)/2, perimetters, 0)

    def check_volumetry (self, src_alt, tgt_alt, volume):
        return (self.pxsize**2)/2 * (tgt_alt/2)

    # def get_current_flood (self, rc, last_flood):
    #     return reduce(
    #         lambda a, i: a - self.drainages[(rc[0] - self.deltas[i][0], rc[1] - self.deltas[i][1])],
    #         [i for i in range(0, 9, 1) if i != 4],
    #         last_flood
    #     )

    def drainpaths (self, start, flow):
        self.drainages = np.zeros(self.dtm.shape)
        queue = list()
        visited = dict()
        
        def _drainpaths (rcs, step_drainages):
            try:
                next_step = list()

                for rc in rcs:
                    if rc in visited: continue                        
                    visited[rc] = True
                        
                    src_flood = step_drainages[rc]
                    if src_flood < 1e-1: continue
                    
                    perimetters = self.get_perimetters(*rc)
                    downslopes = self.get_downslopes(perimetters)
                    if len(downslopes) == 0:
                        src_flood = self.drainages[rc] - self.drainages[(
                            self.deltas[perimetters.argmin()][0],
                            self.deltas[perimetters.argmin()][1]
                        )] + src_flood
                        overflow = src_flood - self.get_volumetry(perimetters)
                        if overflow > 0:
                            delta = self.deltas[perimetters.argmin()]
                            new_rc = (rc[0] + delta[0], rc[1] + delta[1])
                            step_drainages[new_rc] += overflow
                            if new_rc not in visited: next_step.append(new_rc)
                    else:

                        t_downslopes = reduce(lambda a, d: a + min(0, d), downslopes.values())
                        
                        for delta in downslopes:
                            new_rc = (rc[0] + delta[0], rc[1] + delta[1])
                            catchment_factor = downslopes[delta] / t_downslopes
                            new_catchment = catchment_factor * src_flood
                            step_drainages[new_rc] += new_catchment
                            if new_rc not in visited: next_step.append(new_rc)
                    
                if len(next_step): queue.append(next_step)
                if len(queue) > 0: _drainpaths(queue.pop(), step_drainages)
                return step_drainages
                        
            except KeyboardInterrupt as e:
                print("KeyboardInterrupt!!")
                return step_drainages
            except ValueError as e:
                print("ValueError!!")
                print(e)
                return step_drainages
            except RecursionError as e:
                print("RecursionError!!")
                print(e)
                return step_drainages
            except Exception as e:
                print("Exception!!")
                print(e)
                return step_drainages
                
        
        hyd = hydrogram(flow)
        for flood in hyd:
            step_drainages = np.zeros(self.drainages.shape)
            step_drainages[start] = float(flood)
            step_drainages = _drainpaths([start], step_drainages)
            self.drainages = np.where(self.drainages > step_drainages, self.drainages, step_drainages)

        return self.drainages


