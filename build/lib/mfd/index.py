
# BULTINS
import sys
import math
from collections import defaultdict

# VENDOR
import richdem as rd
# from matplotlib import pyplot

# MODULES
from .matrix import Matrix
from .hydrogram import hydrogram
from .debug import print_exception, crono, truncate, progress_bar, progress_counter


class MFD (Matrix):

    def __init__ (self, dtm, mannings, cellsize, radius):
        Matrix.__init__(self, dtm)
            
        self.cellsize = cellsize
        self.radius = radius

        self.dtm = rd.rdarray(self.dtm, no_data=float("nan"))
        rd.FillDepressions(self.dtm, in_place=True)
        self.mannings = self.array(mannings)

    def start_point (self, rc, drafts):
        slopes = self.get_slopes(rc, drafts)
        gateway = slopes.argmin()
        return tuple(rc + self.deltas[gateway]), {
            rc: True
        }, slopes[gateway]

    def get_slopes (self, rc, drafts):
        return self.array([(self.dtm[tuple(delta)] + drafts[tuple(delta)]) - (self.dtm[rc] + drafts[rc]) for delta in rc + self.deltas])

    def get_volumetries (self, slopes, not_visiteds):
        return self.where(not_visiteds, math.pow(self.cellsize, 2.0) * 0.5 * slopes * (1/3), 0)
        
    def get_downslopes (self, slopes, not_visiteds):
        return self.where(self.log_and(not_visiteds, slopes < 0), slopes*-1, 0)

    def get_upslopes (self, slopes, not_visiteds):
        return self.where(self.log_and(not_visiteds, slopes >= 0), slopes, 0)

    def get_draft (self, rc, flood):
        return flood/math.pow(self.cellsize, 2.0)

    def get_speeds (self, slopes, draft, manning, not_visiteds):
        return self.where(not_visiteds, list(map(lambda slope: self.get_speed(draft, manning, slope), slopes)), 0)

    def get_speed (self, draft, manning, slope):
        return max(1e-3, (1.0/manning) * math.pow(self.cellsize+2*draft, 2.0/3.0) * math.pow(max(0, (-1*slope))/5.0, 0.5))
        # return max(1e-3, (1.0/manning) * math.pow((self.cellsize*draft)/(self.cellsize+2*draft), 2.0/3.0) * math.pow(max(0, (-1*slope))/5.0, 0.5))

    @crono
    def drainpaths (self, src, break_flow, base_flow, break_time):
        floods = self.zeros(self.dtm.shape)
        drafts = self.zeros(self.dtm.shape)
        speeds = self.zeros(self.dtm.shape)
        slopes = self.zeros(self.dtm.shape)
        drainages = self.zeros(self.dtm.shape)
        flood_factor = 0
        visited = dict()
        self.is_over = False

        def _drainpaths (rcs, next_step=dict(), queue=list(), level=0):
            try:
                if self.is_over: return
                next_level = list()
                reacheds = list()
                for rc in rcs:
                    if rc in visited: continue
                    
                    src_flood = floods[rc]
                    src_speed = speeds[rc]
                    # src_draft = drafts[rc]
                    src_slope = slopes[rc]

                    if src_speed / self.cellsize < 1 and src_flood / self.cellsize < 0.5:
    
                        if level == 0:
                            floods[rc] += src_flood * flood_factor * min(1, src_speed / self.cellsize)
                            drafts[rc] = self.get_draft(rc, floods[rc])
                            speeds[rc] = self.get_speed(drafts[rc], self.mannings[rc], src_slope)
                        
                        next_step[rc] = True
                        continue

                    rc_slopes = self.get_slopes(rc, drafts)
                    not_visiteds = [True] * 9  # self.array(list(map(lambda d: tuple(d) not in visited and tuple(d) != rc, rc + self.deltas)))
                    downslopes = self.get_downslopes(rc_slopes, not_visiteds)
                    upslopes = self.get_upslopes(rc_slopes, not_visiteds)
                    under_volume = self.get_volumetries(downslopes, not_visiteds)
                    over_volume = self.get_volumetries(upslopes, not_visiteds)

                    if sum(downslopes) == 0:
                        over_flood = max(0, src_flood - over_volume.min() * 8)
                        drived_flood = 0
                        if over_flood == 0:
                            if level == 0:
                                floods[rc] += src_flood * flood_factor * min(1, src_speed / self.cellsize)
                                drafts[rc] = self.get_draft(rc, floods[rc])
                                speeds[rc] = 0
                            next_step[rc] = True
                            continue
                    else:
                        drived_flood = min(src_flood, sum(under_volume))
                        over_flood = src_flood - drived_flood

                    visited[rc] = True
                    if rc in next_step: del next_step[rc]

                    over_cacthments = self.where(src_flood > over_volume * 8, src_flood - over_volume * 8, 0)
                    # CATCHMENT DISTRIBUTION. Powers of catchment and slopes defined as the level of concentration/dispersion
                    # of the floods drived by the slopes.
                    overfloods = over_cacthments ** 1 / sum(over_cacthments ** 1) * over_flood if sum(over_cacthments) else self.zeros((9,))
                    drivedfloods = downslopes ** 2 / sum(downslopes ** 2) * drived_flood if sum(downslopes) else self.zeros((9,))
                    rc_floods = overfloods + drivedfloods
                    rc_speeds = self.get_speeds(rc_slopes, drafts[rc], self.mannings[rc], self.log_and(rc_floods > 0, not_visiteds))

                    rc_acum_flood = sum(rc_floods)
                    rc_acum_flood2 = sum(rc_floods ** 1)
                    rc_acum_speed2 = sum(rc_speeds ** 2)
                    for i, (flood, speed) in enumerate(zip(rc_floods, rc_speeds)):
                        new_rc = tuple(rc + self.deltas[i])
                        if not self.mannings[new_rc] or not self.dtm[new_rc]:
                            self.is_over = True
                            return
                        slopes[new_rc] = slopes[new_rc] or rc_slopes[i] + rc_slopes[i] / 2
                        speeds[new_rc] = (speeds[new_rc] or speed + speed) / 2

                        # CATCHMENT ASSIGNATION. Based on a ponderation of flood by the speed and powered as the level of 
                        # concentration/dispersion drived by the speed.
                        floods[new_rc] += (flood ** 1 / rc_acum_flood2 + speed ** 2 / rc_acum_speed2) / 2 * rc_acum_flood 
                        drafts[new_rc] = self.get_draft(new_rc, floods[new_rc])

                        # DRAINAGE: Define the critical level of flood when the terrain can drain all the
                        # water and it's impossible the accumulate flood.
                        drainages[new_rc] += 1
                        # if floods[new_rc] / self.cellsize < 1e-4: continue
                        if (floods[new_rc] / self.cellsize < 1e-4 and drainages[new_rc] > 10) or floods[new_rc] / self.cellsize < 1e-5:
                            continue

                        if speed / self.cellsize > 1:
                            reacheds.append(new_rc)
                        else:
                            next_level.append(new_rc)
                
                if len(reacheds) > 0: queue.insert(0, reacheds)
                if len(next_level) > 0: queue.append(next_level)
                if len(queue) > 0: _drainpaths(queue.pop(0), next_step, queue, level + 1)
            except Exception:
                print_exception()
            finally:
                return next_step

        try:
            start, visited, slope = self.start_point(src, drafts)
            hyd = hydrogram(break_flow, base_flow, break_time)
            last_flood = None
            floods[start] = break_flow
            drafts[start] = self.get_draft(start, break_flow)
            speeds[start] = self.get_speed(break_flow / self.cellsize ** 2, self.mannings[start], slope)
            next_step = {start: True}
            i = 0
            # steps, news, outs, lens = list(), list(), list(), list()
            # progress = progress_bar(break_time)
            progress = progress_counter("FLOOD")
            for flood in hyd:
                # print(flood)
                progress(i, flood)
                flood_factor = (flood / last_flood) if last_flood else 0
                # last_step = next_step
                next_step = _drainpaths(
                    next_step,
                    dict()
                )
                # outs.append(flood)
                # steps.append(len(next_step))
                # news.append(len(list(filter(lambda k: k not in last_step, next_step))))
                # lens.append(self.array([math.sqrt(sum(coord**2)) for coord in abs(self.argwhere(floods > 0) - start) * self.cellsize]).max())
                distance = self.array([math.sqrt(sum(coord**2)) for coord in abs(self.argwhere(floods > 0) - start) * self.cellsize]).max()
                last_flood = flood
                i += 1
                if self.is_over or distance >= self.radius: break
                
        except KeyboardInterrupt:
            print("KeyboardInterruption!")
        except Exception:
            print("Exception!")
            print_exception()
        finally:
            # pyplot.plot(list(range(0, len(steps), 10)), [math.log(d) if d else d for d in [sum(outs[i:i+10]) for i in range(0, len(steps), 10)]], "b-")
            # pyplot.plot(list(range(0, len(steps), 10)), [math.log(d) if d else d for d in [sum(steps[i:i+10]) for i in range(0, len(steps), 10)]], "r-")
            # pyplot.plot(list(range(0, len(steps), 10)), [math.log(d) if d else d for d in [sum(news[i:i+10]) for i in range(0, len(steps), 10)]], "g--")
            # pyplot.plot(list(range(0, len(steps), 10)), [math.log(d) if d else d for d in [sum(lens[i:i+10]) for i in range(0, len(lens), 10)]], "y--")
            # pyplot.show()
            return floods, drafts, speeds
