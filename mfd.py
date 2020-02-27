# SYS
import sys
import linecache
import math
from functools import reduce

# MODULES
import numpy as np
from osgeo import gdal
import richdem as rd


# sys.setrecursionlimit(3000)
# np.seterr(all="raise")

def print_exception ():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))


def _mat (mat, direction, delta):
    direction = direction.lower()
    tmp = np.zeros(mat.shape)
    if direction == "l" or direction == "left":
        tmp[:, 0:delta] = [float("inf") * mat.shape[1]]*delta
    elif direction == "r" or direction == "right":
        tmp[:, -1:-1*delta] = [float("inf") * mat.shape[1]]*delta
    elif direction == "t" or direction == "top":
        tmp[0:delta, :] = [float("inf") * mat.shape[0]]*delta
    elif direction == "b" or direction == "bottom":
        tmp[-1:-1*delta, :] = [float("inf") * mat.shape[0]]*delta
        
    return tmp
    
    
def displace2top (mat, delta=1):
    tmp = _mat(mat, "T", delta)
    tmp[delta:, :] = mat[:-1*delta, :]
    return tmp
    
    
def displace2bottom (mat, delta=1):
    tmp = _mat(mat, "B", delta)
    tmp[:-1*delta, :] = mat[delta:, :]
    return tmp
    

def displace2right (mat, delta=1):
    tmp = _mat(mat, "R", delta)
    tmp[:, :-1*delta] = mat[:, delta:]
    return tmp
    

def displace2left (mat, delta=1):
    tmp = _mat(mat, "L", delta)
    tmp[:, delta:] = mat[:, :-1*delta]
    return tmp


def hydrogram (f):
    f = float(f)
    
    def _range (start, stop, step):
        x = start
        while x <= stop:
             yield x
             x += step
    
    time = _range(0.0, 1.0, 0.005)
    
    def _hydrogram ():
        for t in time:
            yield f**-t * f + 60

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
            
        self.dtm = rd.rdarray(self.dtm, no_data=float("nan"))
        rd.FillDepressions(self.dtm, in_place=True)
            
        self.pxsize = pxsize
        
        self.slopes = {card: None for card in self.cardinalities if card}
        self.volumetry = np.zeros(self.dtm.shape)
        
        for course in self.slopes:
            if course == "T":
                self.slopes["T"] = displace2top(self.dtm) - self.dtm
                self.volumetry = self.pxsize**2/2 * np.where(self.slopes["T"] > 0, self.slopes["T"], 0)/2
            elif course == "L":
                self.slopes["L"] = displace2left(self.dtm) - self.dtm
                self.volumetry += self.pxsize**2/2 * np.where(self.slopes["L"] > 0, self.slopes["L"], 0)/2
            elif course == "B":
                self.slopes["B"] = displace2bottom(self.dtm) - self.dtm
                self.volumetry += self.pxsize**2/2 * np.where(self.slopes["B"] > 0, self.slopes["B"], 0)/2
            elif course == "R":
                self.slopes["R"] = displace2right(self.dtm) - self.dtm
                self.volumetry += self.pxsize**2/2 * np.where(self.slopes["R"] > 0, self.slopes["R"], 0)/2
            elif course == "RT":
                self.slopes["RT"] = displace2right(displace2top(self.dtm)) - self.dtm
                self.volumetry += math.sqrt(self.pxsize**2+self.pxsize**2)**2/2 * np.where(self.slopes["RT"] > 0, self.slopes["RT"], 0)/2
            elif course == "RB":
                self.slopes["RB"] = displace2right(displace2bottom(self.dtm)) - self.dtm
                self.volumetry += math.sqrt(self.pxsize**2+self.pxsize**2)**2/2 * np.where(self.slopes["RB"] > 0, self.slopes["RB"], 0)/2
            elif course == "LT":
                self.slopes["LT"] = displace2left(displace2top(self.dtm)) - self.dtm
                self.volumetry += math.sqrt(self.pxsize**2+self.pxsize**2)**2/2 * np.where(self.slopes["LT"] > 0, self.slopes["LT"], 0)/2
            elif course == "LB":
                self.slopes["LB"] = displace2left(displace2bottom(self.dtm)) - self.dtm
                self.volumetry += math.sqrt(self.pxsize**2+self.pxsize**2)**2/2 * np.where(self.slopes["LB"] > 0, self.slopes["LB"], 0)/2

    def start_point (self, rc):
        perimetters = self.get_perimetters(rc)
        gateway = perimetters.argmin()
        return (rc[0] + self.deltas[gateway][0], rc[1] + self.deltas[gateway][1]), {rc: True, **{
            (rc[0] + delta[0], rc[1] + delta[1]): True
            for i, delta in enumerate(self.deltas)
            if i != gateway
        }}

    def get_perimetters (self, rc):
        return np.array([
            self.slopes.get(cardinality)[rc] if cardinality else float("inf")
            for cardinality in self.cardinalities
        ])

    def get_overflows (self, flood, perimetters, not_visiteds):
        volumetries = np.array([
            (i%2 == 0 and (self.pxsize**2/2) or (math.sqrt(self.pxsize**2+self.pxsize**2)**2)) * (p * self.pxsize)/2 * 8
            for i, p in enumerate(perimetters)
        ])
        return np.where(np.logical_and(not_visiteds, volumetries < flood), flood - volumetries, 0)

    def get_downslopes (self, perimetters):
        return np.where(perimetters < 0, perimetters*-1, 0)

    def get_plains (self, perimetters):
        return np.where(perimetters == 0, 1, 0)

    def drainpaths (self, src, flow):
        self.drainages = np.zeros(self.dtm.shape)
        self.last_step = self.drainages.copy()

        def _drainpaths (rcs, step_drainages, queue=list(), visited=dict()):
            try:
                next_step = list()
                for rc in rcs:
                    if rc in visited: continue                     
                    visited[rc] = True
                    
                    src_flood = step_drainages[rc]

                    if src_flood < 1e-2:
                        src_flood += min(self.last_step[rc], self.volumetry[rc])
                        if src_flood < 1e-2: continue
                       
                    perimetters = self.get_perimetters(rc)
                    downslopes = self.get_downslopes(perimetters)
                    plains = self.get_plains(perimetters)

                    if downslopes.sum() + plains.sum() == 0:
                        print("LOCAL DEPRESSION!!")
                        src_flood += min(self.last_step[rc], self.volumetry[rc])

                    overstream = src_flood - self.volumetry[rc]
                    downstream = src_flood - overstream
                    print("\n")
                    print("RC: ", rc, "SRCFLOOD: ", src_flood)
                    print("PERIMETTERS: ", perimetters)
                    print("VOLUME: ", self.volumetry[rc], "OVERSTREAM: ", overstream, "DOWNSTREAM: ", downstream)

                    overflows = self.get_overflows(src_flood, perimetters, np.array([
                        (rc[0] + delta[0], rc[1] + delta[1]) not in visited
                        for i, delta in enumerate(self.deltas)
                    ]))

                    total_overflows = overflows.sum()
                    if total_overflows > 0:
                        overstreams = overflows/total_overflows * overstream
                    else:
                        overstreams = np.zeros(overflows.shape)

                    total_downslopes = downslopes.sum()
                    if total_downslopes > 0:
                        downstreams = downslopes/total_downslopes * downstream
                    else:
                        downstreams = np.zeros(downslopes.shape)

                    total_plains = plains.sum()
                    if total_plains > 0:
                        plainstreams = plains/total_plains * downstream
                    else:
                        plainstreams = np.zeros(plains.shape)

                    print("OVERFLOWS: ", len(overflows[overflows > 0]), "DOWNSLOPES: ", len(downslopes[downslopes > 0]))
                    print("OVERSTREAMS: ", overstreams)
                    print("DOWNSTREAMS: ", downstreams)
                    print("TOTAL_OVERSTREAMS: ", overstreams.sum(), "TOTAL_DOWNSTREAMS: ", downstreams.sum())

                    for i, delta in enumerate(zip(overstreams, downstreams, plainstreams)):
                        new_rc = (rc[0] + self.deltas[i][0], rc[1] + self.deltas[i][1])
                        print("NEWRC: ", new_rc, "OVERSTREAM: ", overstream, "DOWNSTREAM: ", downstream)
                        if new_rc in visited: continue
                        step_drainages[new_rc] += delta[0] + delta[1] + delta[2]
                        next_step.append(new_rc)

                if len(next_step) > 0: queue.append(next_step)
                if len(queue) > 0: _drainpaths([rc for rc in queue.pop() if rc not in visited], step_drainages, queue, visited)

            except ValueError:
                print("ValueError!!")
                print_exception()
            except RecursionError:
                print("RecursionError!!")
                print_exception()
            except Exception:
                print("Exception!!")
                print_exception()
            finally:
                self.last_step = step_drainages.copy()
                return step_drainages

        try:
            start, visited = self.start_point(src)
            # self.drainages[start] = float(flow)
            # self.drainages = _drainpaths([start], self.drainages, list(), visited)
            hyd = hydrogram(flow)
            for flood in hyd:
                step_drainages = np.zeros(self.drainages.shape)
                step_drainages[start] = float(flood)
                step_drainages = _drainpaths(
                    [start],
                    step_drainages,
                    list(),
                    dict(visited)
                )
                self.drainages = np.where(self.drainages >= step_drainages, self.drainages, step_drainages)
        except KeyboardInterrupt:
            print("KeyboardInterruption!!")
        finally:
            return self.drainages


