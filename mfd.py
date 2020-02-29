
# SYS
from time import time, sleep
import sys
import linecache
import math

# MODULES
import numpy as np
import richdem as rd


# sys.setrecursionlimit(3000)
# np.seterr(all="raise")
# np.set_printoptions(precision=2)
# def print (*args, **kwargs):
#     pass
# def input (*args, **kwargs):
#     pass
# def debugger (msg):
#     def wrapper (fn):
#         def _wrapper (*args, **kwargs):
#             res = fn(*args, **kwargs)
#             print("[[DEBUG]]")
#             print(" DESCRIPTION: ", msg)
#             print(" ARGS: ", args)
#             print(" KWARGS: ", kwargs)
#             print(" RES: ", res)
#             input("Press Enter to continue...")
#             return res
#         return _wrapper
#     return wrapper

count = 0
def counter (fn):
    # count = [0]
    global count
    def wrapper (*args, **kwargs):
        global count
        start = time()
        count += 1
        res = fn(*args, **kwargs)
        print("COUNT: ", count)
        print(time() - start)
        return res
    return wrapper


def print_exception ():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))


def hydrogram (f):
    f = float(f)
    
    def _range (start, stop, step):
        x = start
        while x <= stop:
             yield x
             x += step
    
    timerange = _range(0.0, 1.0, 1e-2)
    
    def _hydrogram ():
        for t in timerange:
            yield f**-t * f + 60

    return _hydrogram()


class MFD (object):

    def __init__ (self, dtm_array, pxsize):
        self.deltas = np.array([
            (-1, -1), (-1, 0), (-1, +1),
            (0, -1), (0, 0), (0, +1),
            (+1, -1), (+1, 0), (+1, +1)
        ])
        
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

    def start_point (self, rc):
        perimetters = self.get_perimetters(rc, [True]*9)
        gateway = perimetters.argmin()
        return tuple(rc + self.deltas[gateway]), {rc: True, **{
            tuple(delta): True
            for i, delta in enumerate(rc + self.deltas)
            if i != gateway
        }}

    def get_perimetters (self, rc, not_visiteds):
        return np.where(not_visiteds, [self.dtm[tuple(delta)] - self.dtm[rc] for delta in rc + self.deltas], float("inf"))

    def get_volumetries (self, flood, perimetters, not_visiteds):
        return np.where(not_visiteds, [
            (i%2 == 0 and self.pxsize or math.sqrt(self.pxsize**2+self.pxsize**2))**2 * max(0, p) * .5 * (1/3) * 4
            for i, p in enumerate(perimetters)
        ], 0)

    def get_downslopes (self, perimetters, not_visiteds):
        return np.where(np.logical_and(not_visiteds, perimetters < 0), perimetters*-1, 0)

    def get_plains (self, perimetters, not_visiteds):
        return np.where(np.logical_and(not_visiteds, perimetters == 0), 1, 0)

    def drainpaths (self, src, flow):
        self.drainages = np.zeros(self.dtm.shape)
        self.last_step = np.zeros(self.dtm.shape)

        # @debugger("drainpaths")
        # @counter
        def _drainpaths (rcs, step_drainages, queue, visited):
            new_step = list()
            for rc in rcs:
                if rc in visited: continue

                visited[rc] = True
                src_flood = step_drainages[rc]

                if src_flood < 1e-2:
                    src_flood += self.last_step[rc]
                    if src_flood < 1e-2:
                        step_drainages[rc] = src_flood
                        continue
                    
                not_visiteds = list(map(lambda d: tuple(d) not in visited, rc + self.deltas)) # [True]*9
                perimetters = self.get_perimetters(rc, not_visiteds)
                downslopes = self.get_downslopes(perimetters, not_visiteds)
                plains = self.get_plains(perimetters, not_visiteds)
                volumetries = self.get_volumetries(src_flood, perimetters, not_visiteds)

                # print("\nRC: ", rc)
                # print("SRCFLOOD: ", src_flood)
                # print("!VISITEDS: ", not_visiteds)
                # print("PERIMETTERS: ", perimetters)
                # print("DOWNSLOPES: ", downslopes)
                # print("PLAINS: ", plains)
                # print("VOLUMETRIES: ", volumetries)

                if downslopes.sum() + plains.sum() == 0:
                    src_flood += min(self.last_step[rc], max(0, volumetries.min()))
                    over_flood = min(0, src_flood - max(0, volumetries.min()))
                    step_drainages[rc] = src_flood - over_flood
                    drived_flood = 0
                else:
                    over_flood = min(0, src_flood - max(0, volumetries.min()))
                    drived_flood = src_flood - over_flood

                overflows = np.where(src_flood > volumetries, src_flood - volumetries, 0)
                overstreams = overflows/overflows.sum() * over_flood if overflows.sum() else np.zeros(volumetries.shape)
                downstreams = downslopes/downslopes.sum() * drived_flood if downslopes.sum() else np.zeros(volumetries.shape)
                plainstreams = plains/plains.sum() * drived_flood if downslopes.sum() == 0 and plains.sum() > 0 else np.zeros(plains.shape)

                # print("OVERFLOWS: ", overflows)
                # print("DOWNSTREAMS: ", downstreams)
                # print("PLAINSTREAMS: ", plainstreams)

                for i, delta in enumerate(zip(overstreams, downstreams, plainstreams)):
                    new_rc = tuple(rc + self.deltas[i])
                    catchment = sum(delta) # delta[0] + delta[1] + delta[2]
                    if catchment == 0 or new_rc in visited: continue
                    # print("NEWRC: ", new_rc, "CATCHMENT: ", catchment)
                    step_drainages[new_rc] += catchment
                    new_step.append(new_rc)

            # input("Press Enter to continue...")
            if len(new_step) > 0: queue.append(new_step)
            if len(queue) > 0: _drainpaths(queue.pop(0), step_drainages, queue, visited)

            return step_drainages

        try:
            start, visited = self.start_point(src)
            # self.drainages[start] = float(flow)
            # self.drainages = _drainpaths(start, self.drainages, list(), visited)
            hyd = hydrogram(flow)
            step_drainages = np.zeros(self.drainages.shape)
            total_flood = 0
            for flood in hyd:
                # print(flood)
                total_flood += flood
                step_drainages[:,:] = 0
                step_drainages[start] = float(flood)
                step_drainages = _drainpaths(
                    [start],
                    step_drainages,
                    list(),
                    dict(visited)
                )
                self.last_step = step_drainages.copy()
                self.drainages = np.where(self.drainages >= step_drainages, self.drainages, step_drainages)
                
        except KeyboardInterrupt:
            print("KeyboardInterruption!!")
        except Exception:
            print("Exception!")
            print_exception()
        finally:
            print("INT: ", total_flood)
            print("OUT: ", self.drainages.sum())
            return self.drainages


