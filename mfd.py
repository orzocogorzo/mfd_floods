
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

    def __init__ (self, dtm_array, cellsize):
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
            
        self.cellsize = cellsize

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
            (i%2 == 0 and self.cellsize or math.sqrt(self.cellsize**2+self.cellsize**2))**2 * p * .5 * (1/3) * 4
            for i, p in enumerate(perimetters)
        ], 0)

    def get_downslopes (self, perimetters, not_visiteds):
        return np.where(np.logical_and(not_visiteds, perimetters < 0), perimetters*-1, 0)

    def drainpaths (self, src, flow):
        self.drainages = np.zeros(self.dtm.shape)
        self.last_step = np.zeros(self.dtm.shape)

        def _drainpaths (rcs, step_drainages, queue, visited):
            new_step = list()
            for rc in rcs:
                if rc in visited: continue

                visited[rc] = True
                src_flood = step_drainages[rc]

                if src_flood < 1e-2:
                    src_flood += self.last_step[rc]
                    step_drainages[rc] = src_flood
                    if src_flood < 1e-2: continue
                    
                not_visiteds = list(map(lambda d: tuple(d) not in visited, rc + self.deltas)) # [True]*9
                perimetters = self.get_perimetters(rc, not_visiteds)
                downslopes = self.get_downslopes(perimetters, not_visiteds)
                volumetries = self.get_volumetries(src_flood, perimetters, not_visiteds)
                lessvol = max(0, volumetries.min())

                if downslopes.sum() == 0:
                    src_flood += min(self.last_step[rc], lessvol)
                    over_flood = max(0, src_flood - lessvol)
                    step_drainages[rc] = src_flood - over_flood
                    drived_flood = 0
                else:
                    over_flood = max(0, src_flood - lessvol)
                    drived_flood = src_flood - over_flood

                overflows = np.where(src_flood > volumetries, src_flood - volumetries, 0)
                overstreams = overflows/overflows.sum() * over_flood if overflows.sum() else np.zeros(volumetries.shape)
                drivedstreams = downslopes/downslopes.sum() * drived_flood if downslopes.sum() else np.zeros(volumetries.shape)

                # print("\nRC: ", rc)
                # print("SRCFLOOD: ", src_flood)
                # print("!VISITEDS: ", not_visiteds)
                # print("PERIMETTERS: ", perimetters)
                # print("DOWNSLOPES: ", downslopes)
                # print("VOLUMETRIES: ", volumetries)
                # print("OVERFLOOD: ", over_flood)
                # print("OVERFLOWS: ", overstreams)
                # print("DOWNSTREAMS: ", drivedstreams)

                for i, delta in enumerate(zip(overstreams, drivedstreams)):
                    new_rc = tuple(rc + self.deltas[i])
                    catchment = sum(delta)
                    if catchment == 0 or new_rc in visited: continue
                    # print("NEWRC: ", new_rc, "CATCHMENT: ", catchment)
                    step_drainages[new_rc] += catchment
                    new_step.append(new_rc)

            # input("Press Enter to continue...")
            # if len(new_step) > 0: _drainpaths(new_step, step_drainages, queue, visited)
            if len(new_step) > 0: queue.append(new_step)
            if len(queue) > 0: _drainpaths(queue.pop(0), step_drainages, queue, visited)

            return step_drainages

        try:
            start, visited = self.start_point(src)
            # self.drainages[start] = float(flow)
            # self.drainages = _drainpaths(
            #     [start],
            #     self.drainages,
            #     list(),
            #     dict(visited)
            # )
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




#       def drainpaths (self, start, flow):
#         self.drainages = np.zeros(self.dtm.shape)
#         # queue = list()
#         # visited = dict()
        
#         def _drainpaths (rcs, step_drainages, queue=list(), visited=dict()):
#             try:
#                 next_step = list()

#                 for rc in rcs:
#                     if rc in visited: continue                        
#                     visited[rc] = True
                        
#                     src_flood = step_drainages[rc]
#                     if src_flood < 1e-1: continue
                    
#                     perimetters = self.get_perimetters(rc)
#                     downslopes = self.get_downslopes(perimetters)
#                     if len(downslopes) == 0:
#                         src_flood = self.drainages[rc] - self.drainages[(
#                             self.deltas[perimetters.argmin()][0],
#                             self.deltas[perimetters.argmin()][1]
#                         )] + src_flood
#                         overflow = src_flood - self.get_volumetry(perimetters)
#                         if overflow > 0:
#                             delta = self.deltas[perimetters.argmin()]
#                             new_rc = (rc[0] + delta[0], rc[1] + delta[1])
#                             step_drainages[new_rc] += overflow
#                             if new_rc not in visited: next_step.append(new_rc)
#                     else:

#                         t_downslopes = min(sum(downslopes.values()), 0)
                        
#                         for delta in downslopes:
#                             new_rc = (rc[0] + delta[0], rc[1] + delta[1])
#                             catchment_factor = downslopes[delta] / t_downslopes
#                             new_catchment = catchment_factor * src_flood
#                             step_drainages[new_rc] += new_catchment
#                             if new_rc not in visited: next_step.append(new_rc)
                    
#                 if len(next_step): queue.append(next_step)
#                 if len(queue): _drainpaths(queue.pop(), step_drainages, queue, visited)
#                 return step_drainages
                        
#             except KeyboardInterrupt as e:
#                 print("KeyboardInterrupt!!")
#                 return step_drainages
#             except ValueError as e:
#                 print("ValueError!!")
#                 print(e)
#                 return step_drainages
#             except RecursionError as e:
#                 print("RecursionError!!")
#                 print(e)
#                 return step_drainages
#             except Exception as e:
#                 print("Exception!!")
#                 print(e)
#                 return step_drainages
                
        
#         hyd = hydrogram(flow)
#         for flood in hyd:
#             step_drainages = np.zeros(self.drainages.shape)
#             step_drainages[start] = float(flood)
#             step_drainages = _drainpaths([start], step_drainages, list(), dict())
#             self.drainages = np.where(self.drainages > step_drainages, self.drainages, step_drainages)

#         return self.drainages