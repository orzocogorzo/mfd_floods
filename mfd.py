
# SYS
import math
from collections import defaultdict

# VENDOR
import richdem as rd

# MODULES
from matrix import Matrix
from hydrogram import hydrogram
from debug import print_exception, crono, truncate


class MFD (Matrix):

    def __init__ (self, dtm_array, manning_array, cellsize):
        Matrix.__init__(self, dtm_array)
            
        self.cellsize = cellsize

        self.dtm = rd.rdarray(self.dtm, no_data=float("nan"))
        # rd.FillDepressions(self.dtm, in_place=True)
        self.mannings = self.array(manning_array)

    def start_point (self, rc):
        slopes = self.get_slopes(rc, [True]*9)
        gateway = slopes.argmin()
        return tuple(rc + self.deltas[gateway]), {
            rc: True
        }, slopes[gateway]

    def get_slopes (self, rc, not_visiteds):
        return self.array([self.dtm[tuple(delta)] - self.dtm[rc] for delta in rc + self.deltas])

    def get_volumetries (self, slopes, not_visiteds):
        return self.where(not_visiteds, math.pow(self.cellsize*0.5, 2.0) * 0.5 * slopes/2 * (1/3), 0)
        
    def get_downslopes (self, slopes, not_visiteds):
        return self.where(self.log_and(not_visiteds, slopes < 0), slopes*-1, 0)

    def get_upslopes (self, slopes, not_visiteds):
        return self.where(self.log_and(not_visiteds, slopes >= 0), slopes, 0)

    def get_draft (self, rc, flood):
        return flood/math.pow(self.cellsize, 2.0)

    def get_speeds (self, slopes, draft, manning, incoming_speed, not_visiteds):
        return self.where(not_visiteds, list(map(lambda slope: (self.get_speed(draft, manning, slope) + incoming_speed)/2, slopes)), 0)

    def get_speed (self, draft, manning, slope):
        return max(1e-2, (1.0/manning) * math.pow((self.cellsize*draft)/(self.cellsize+2*draft), 2.0/3.0) * math.pow(max(0, (-1*slope))/5.0, 0.5))

    @crono
    def drainpaths (self, src, break_flow, base_flow, break_time):
        floods = self.zeros(self.dtm.shape)
        drafts = self.zeros(self.dtm.shape)
        speeds = self.zeros(self.dtm.shape)
        flood_factor = 0
        visited = dict()

        def _drainpaths (rcs, next_step=dict(), queue=list(), level=0):
            try:
                next_level = list()
                for rc in rcs:
                    if rc in visited or not self.mannings[rc] or not self.dtm[rc]: continue
                    
                    src_flood = floods[rc]

                    if src_flood < self.cellsize*0.02:
                        # print("NOT ENOUGHT")
                        if level == 0:
                            floods[rc] += src_flood * flood_factor
                            drafts[rc] = self.get_draft(rc, floods[rc])
                        next_step[rc] = True
                        continue

                    not_visiteds = self.array(list(map(lambda d: tuple(d) not in visited and tuple(d) != rc, rc + self.deltas)))
                    slopes = self.get_slopes(rc, not_visiteds)
                    downslopes = self.get_downslopes(slopes, not_visiteds)
                    upslopes = self.get_upslopes(slopes, not_visiteds)
                    under_volume = self.get_volumetries(downslopes, not_visiteds)
                    over_volume = self.get_volumetries(upslopes, not_visiteds)
                    incoming_speed = self.where(not_visiteds != True, [speeds[tuple(delta)] for delta in rc + self.deltas], 0).mean()

                    if downslopes.sum() == 0:
                        over_flood = max(0, src_flood - self.where(not_visiteds, over_volume, float("inf")).min()*8)
                        drived_flood = 0
                        # traped_flood = src_flood - over_flood
                        if over_flood == 0:
                            if level == 0:
                                floods[rc] += src_flood * flood_factor
                                drafts[rc] = self.get_draft(rc, floods[rc])
                                speeds[rc] = 0
                            next_step[rc] = True
                            continue
                    else:
                        drived_flood = min(src_flood, under_volume.sum())
                        over_flood = src_flood - drived_flood
                        # traped_flood = 0

                    over_cacthments = self.where(self.log_and(not_visiteds, src_flood > over_volume*8), src_flood - over_volume*8, 0)
                    overfloods = over_cacthments**1.1/(over_cacthments**1.1).sum() * over_flood if over_cacthments.sum() else self.zeros(slopes.shape)
                    drivedfloods = downslopes**1.1/(downslopes**1.1).sum() * drived_flood if downslopes.sum() else self.zeros(slopes.shape)
                    rc_floods = overfloods + drivedfloods
                    rc_speeds = self.get_speeds(slopes, self.get_draft(rc, src_flood), self.mannings[rc], incoming_speed, not_visiteds)

                    # evacuation = 0
                    for i, (flood, speed) in enumerate(zip(rc_floods, rc_speeds)):
                        new_rc = tuple(rc + self.deltas[i])
                        if flood == 0 or speed == 0 or new_rc in visited or not self.mannings[new_rc] or not self.dtm[new_rc]: continue 
                        drafts[new_rc] = self.get_draft(new_rc, flood)
                        speeds[new_rc] = (speeds[new_rc] or speed + speed)/2 
                        floods[new_rc] += (flood/rc_floods.sum()+speed/rc_speeds[rc_floods > 0].sum())/2 * rc_floods.sum() * max(self.cellsize, speeds[new_rc])/self.cellsize
                        # evacuation += floods[new_rc]
                        next_level.append(new_rc)

                    # print(truncate(src_flood), truncate(evacuation))
                    visited[rc] = True
                    if rc in next_step: del next_step[rc]

                if len(next_level) > 0: queue.append(next_level)
                if len(queue) > 0: _drainpaths(queue.pop(0), next_step, queue, level+1)
            except Exception:
                print_exception()
            finally:
                return next_step

        try:
            start, visited, slope = self.start_point(src)
            hyd = hydrogram(break_flow, base_flow, break_time)
            last_flood = None
            floods[start] = break_flow
            drafts[start] = self.get_draft(start, break_flow)
            speeds[start] = self.get_speed(break_flow/self.cellsize**2, self.mannings[start], slope)
            next_step = {start: True}
            for flood in hyd:
                print(truncate(flood))
                flood_factor = (flood/last_flood) if last_flood else 0
                next_step = _drainpaths(
                    next_step,
                    dict()
                )
                last_flood = flood
                
        except KeyboardInterrupt:
            print("KeyboardInterruption!")
        except Exception:
            print("Exception!")
            print_exception()
        finally:
            return floods, drafts, speeds