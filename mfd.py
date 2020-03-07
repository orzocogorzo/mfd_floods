
# SYS
import math

# VENDOR
import richdem as rd

# MODULES
from matrix import Matrix
from hydrogram import hydrogram
from debug import print_exception, crono


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

    @crono
    def drainpaths (self, src, break_flow, base_flow, break_time):
        self.drainages = self.zeros(self.dtm.shape)

        def _drainpaths (rcs, queue=list()):
            try:
                next_level = list()
                for rc in rcs:
                    if rc in self.visited: continue

                    src_flood = self.drainages[rc] + self.drainages[rc] * self.flood_factor

                    if src_flood < 1e-2:
                        self.drainages[rc] += src_flood
                        src_flood = self.drainages[rc]
                        if src_flood < 1e-2:
                            next_step[rc] = True
                            continue
                        
                    not_visiteds = list(map(lambda d: tuple(d) not in self.visited and tuple(d) != rc, rc + self.deltas)) # [True]*9
                    slopes = self.get_slopes(rc, not_visiteds)
                    downslopes = self.get_downslopes(slopes, not_visiteds)
                    volumetries = self.get_volumetries(src_flood, slopes, not_visiteds)
                    lessvol = max(0, volumetries.min())

                    if downslopes.sum() == 0:
                        src_flood += min(self.drainages[rc], lessvol)
                        over_flood = max(0, src_flood - lessvol)
                        self.drainages[rc] = src_flood - over_flood
                        drived_flood = 0
                        if over_flood == 0:
                            next_step[rc] = True
                            continue
                    else:
                        over_flood = max(0, src_flood - lessvol)
                        drived_flood = src_flood - over_flood

                    overflows = self.where(src_flood > volumetries, src_flood - volumetries, 0)
                    overstreams = overflows/overflows.sum() * over_flood if overflows.sum() else self.zeros(volumetries.shape)
                    drivedstreams = downslopes/downslopes.sum() * drived_flood if downslopes.sum() else self.zeros(volumetries.shape)

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
                        if catchment == 0 or new_rc in self.visited: continue
                        # print("NEWRC: ", new_rc, "CATCHMENT: ", catchment)
                        self.drainages[new_rc] += catchment
                        next_level.append(new_rc)

                    self.visited[rc] = True
                    if rc in next_step: del next_step[rc]

                if len(next_level) > 0: queue.append(next_level)
                if len(queue) > 0: _drainpaths(queue.pop(0), queue)
            except Exception:
                print_exception()
            finally:
                # print("NEXTSTEP: ", next_step)
                return next_step

        try:
            start, self.visited = self.start_point(src)
            hyd = hydrogram(break_flow, base_flow, break_time)
            last_flood = None
            self.drainages[start] = break_flow
            next_step = {start: True}
            for flood in hyd:
                print(flood)
                self.flood_factor = flood / last_flood if last_flood else 0
                next_step = _drainpaths(
                    dict(next_step)
                )
                last_flood = flood
                # exit(1)
                
        except KeyboardInterrupt:
            print("KeyboardInterruption!!")
        except Exception:
            print("Exception!")
            print_exception()
        finally:
            return self.drainages