
# SYS
import math

# VENDOR
import richdem as rd

# MODULES
from matrix import Matrix
from hydrogram import hydrogram
from debug import print_exception


class MFD (Matrix):

    def __init__ (self, dtm_array, cellsize):
        Matrix.__init__(self, dtm_array)
            
        self.cellsize = cellsize

        self.dtm = rd.rdarray(self.dtm, no_data=float("nan"))
        rd.FillDepressions(self.dtm, in_place=True)

    def start_point (self, rc):
        slopes = self.get_slopes(rc, [True]*9)
        gateway = slopes.argmin()
        return tuple(rc + self.deltas[gateway]), {rc: True, **{
            tuple(delta): True
            for i, delta in enumerate(rc + self.deltas)
            if i != gateway
        }}

    def get_slopes (self, rc, not_visiteds):
        return self.where(not_visiteds, [self.dtm[tuple(delta)] - self.dtm[rc] for delta in rc + self.deltas], float("inf"))

    def get_volumetries (self, flood, slopes, not_visiteds):
        return self.where(not_visiteds, [
            (i%2 == 0 and self.cellsize or math.sqrt(self.cellsize**2+self.cellsize**2))**2 * p * .5 * (1/3) * 4
            for i, p in enumerate(slopes)
        ], 0)

    def get_downslopes (self, slopes, not_visiteds):
        return self.where(self.log_and(not_visiteds, slopes < 0), slopes*-1, 0)

    def drainpaths (self, src, flow):
        self.drainages = self.zeros(self.dtm.shape)
        self.last_step = self.zeros(self.dtm.shape)
        self.worsts = dict()

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
                slopes = self.get_slopes(rc, not_visiteds)
                downslopes = self.get_downslopes(slopes, not_visiteds)
                volumetries = self.get_volumetries(src_flood, slopes, not_visiteds)
                lessvol = max(0, volumetries.min())

                if downslopes.sum() == 0:
                    src_flood += min(self.last_step[rc], lessvol)
                    over_flood = max(0, src_flood - lessvol)
                    step_drainages[rc] = src_flood - over_flood
                    drived_flood = 0
                else:
                    over_flood = max(0, src_flood - lessvol)
                    drived_flood = src_flood - over_flood

                overflows = self.where(src_flood > volumetries, src_flood - volumetries, 0)
                overstreams = overflows/overflows.sum() * over_flood if overflows.sum() else self.zeros(volumetries.shape)
                drivedstreams = downslopes/downslopes.sum() * drived_flood if downslopes.sum() else self.zeros(volumetries.shape)

                # print("\nRC: ", rc)
                # print("SRCFLOOD: ", src_flood)
                # print("!VISITEDS: ", not_visiteds)
                # print("slopes: ", slopes)
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
            step_drainages = self.zeros(self.drainages.shape)
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
                self.drainages = self.where(self.drainages >= step_drainages, self.drainages, step_drainages)
                
        except KeyboardInterrupt:
            print("KeyboardInterruption!!")
        except Exception:
            print("Exception!")
            print_exception()
        finally:
            print("INT: ", total_flood)
            print("OUT: ", self.drainages.sum())
            return self.drainages