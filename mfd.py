
# SYS
import math
from collections import defaultdict

# VENDOR
import richdem as rd

# MODULES
from matrix import Matrix
from hydrogram import hydrogram
from debug import print_exception, crono


class MFD (Matrix):

    def __init__ (self, dtm_array, manning_array, cellsize):
        Matrix.__init__(self, dtm_array)
            
        self.cellsize = cellsize

        self.dtm = rd.rdarray(self.dtm, no_data=float("nan"))
        self.mannings = Matrix.array(self, manning_array)
        rd.FillDepressions(self.dtm, in_place=True)

    def start_point (self, rc):
        slopes = self.get_slopes(rc, [True]*9)
        gateway = slopes.argmin()
        return tuple(rc + self.deltas[gateway]), {rc: True, **{
            tuple(delta): True
            for i, delta in enumerate(rc + self.deltas)
            if i != gateway
        }}, slopes[gateway]

    def get_slopes (self, rc, not_visiteds):
        return self.array([self.dtm[tuple(delta)] - self.dtm[rc] for delta in rc + self.deltas])

    def get_incoming_speeds (self, rc, not_visiteds):
        return self.where(not_visiteds != True, [self.speeds[tuple(delta)] for delta in rc + self.deltas], 0)

    def get_volumetries (self, slopes, not_visiteds):
        return self.where(not_visiteds, [
            math.pow(self.cellsize*0.5, 2.0) * 0.5 * p * (1/3) * 4
            for i, p in enumerate(slopes)
        ], 0)
        
    def get_downslopes (self, slopes, not_visiteds):
        return self.where(self.log_and(not_visiteds, slopes < 0), slopes*-1, 0)

    def get_draft (self, rc):
        return self.floods[rc]/math.pow(self.cellsize, 2.0)

    def get_speeds (self, slopes, draft, manning, incoming_speed, not_visiteds):
        return self.where(not_visiteds, list(map(lambda slope: (self.get_speed(draft, manning, slope) + incoming_speed)/2, slopes)), 0)

    def get_speed (self, draft, manning, slope):
        return max(1e-001, (1.0/manning) * math.pow(self.cellsize + 2*draft, 2.0/3.0) * math.pow(max(0, (-1*slope))/5.0, 0.5))

    @crono
    def drainpaths (self, src, break_flow, base_flow, break_time):
        self.floods = self.zeros(self.dtm.shape)
        self.drafts = self.zeros(self.dtm.shape)
        self.speeds = self.zeros(self.dtm.shape)

        def _drainpaths (rcs, queue=list(), level=0):
            try:
                next_level = list()
                for rc in rcs:
                    if rc in self.visited: continue
                    if not self.mannings[rc] or not self.dtm[rc]: continue
                        
                    src_flood = self.floods[rc] + self.floods[rc] * self.flood_factor

                    if src_flood < 1e-2:
                        self.floods[rc] = src_flood
                        self.drafts[rc] = self.get_draft(rc)
                        next_step[rc] = True
                        continue
                        
                    if self.speeds[rc] < 5.0:
                        if level > 0:
                            next_step[rc] = True
                            continue
                        
                        self.floods[rc] = src_flood * (self.speeds[rc]/5)
                        src_flood = self.floods[rc]
                        self.drafts[rc] = self.get_draft(rc)

                    not_visiteds = self.array(list(map(lambda d: tuple(d) not in self.visited and tuple(d) != rc, rc + self.deltas))) # [True]*9
                    slopes = self.get_slopes(rc, not_visiteds)
                    downslopes = self.get_downslopes(slopes, not_visiteds)
                    volumetries = self.get_volumetries(slopes, not_visiteds)
                    less_volumetry = max(0, volumetries.min())
                    incoming_speeds = self.get_incoming_speeds(rc, not_visiteds)
                    incoming_speed =  incoming_speeds.sum()/len(incoming_speeds > 0)

                    if downslopes.sum() == 0:
                        over_flood = max(0, src_flood - less_volumetry)
                        drived_flood = 0
                        if over_flood == 0:
                            self.floods[rc] = src_flood
                            self.drafts[rc] = self.get_draft(rc)
                            self.speeds[rc] = 0
                            next_step[rc] = True
                            continue
                    else:
                        over_flood = max(0, src_flood - less_volumetry)
                        drived_flood = src_flood - over_flood

                    overflows = self.where(src_flood > volumetries, src_flood - volumetries, 0)
                    overstreams = overflows/overflows.sum() * over_flood if overflows.sum() else self.zeros(volumetries.shape)
                    drivedstreams = downslopes/downslopes.sum() * drived_flood if downslopes.sum() else self.zeros(volumetries.shape)
                    flows = overstreams + drivedstreams
                    speeds = self.get_speeds(slopes, self.get_draft(rc), self.mannings[rc], incoming_speed, not_visiteds)

                    for i, (flow, speed) in enumerate(zip(flows, speeds)):
                        new_rc = tuple(rc + self.deltas[i])
                        if flow == 0 or speed == 0 or new_rc in self.visited or not self.mannings[new_rc] or not self.dtm[new_rc]: continue
                        self.floods[new_rc] += flow * (flow / flows.sum() + speed / speeds[flows > 0].sum())/2
                        self.drafts[new_rc] = self.get_draft(new_rc)
                        self.speeds[new_rc] = (self.speeds[new_rc] or speed + speed)/2
                        next_level.append(new_rc)

                    self.visited[rc] = True
                    if rc in next_step: del next_step[rc]

                if len(next_level) > 0: queue.append(next_level)
                if len(queue) > 0: _drainpaths(queue.pop(0), queue, level+1)
            except Exception:
                print_exception()
            finally:
                return next_step

        try:
            start, self.visited, slope = self.start_point(src)
            hyd = hydrogram(break_flow, base_flow, break_time)
            last_flood = None
            self.floods[start] = break_flow
            self.drafts[start] = self.get_draft(start)
            self.speeds[start] = self.get_speed(break_flow/self.cellsize**2, self.mannings[start], slope)
            next_step = {start: True}
            for flood in hyd:
                # print(flood)
                self.flood_factor = (flood / last_flood) if last_flood else 0
                next_step = _drainpaths(
                    dict(next_step)
                )
                last_flood = flood
                
        except KeyboardInterrupt:
            print("KeyboardInterruption!")
        except Exception:
            print("Exception!")
            print_exception()
        finally:
            return self.floods, self.drafts, self.speeds