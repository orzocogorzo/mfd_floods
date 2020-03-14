
# BULTINS
import sys
import math
from collections import defaultdict

# VENDOR
import richdem as rd

# MODULES
from matrix import Matrix
from hydrogram import hydrogram
from debug import print_exception, crono, truncate, progress_bar


class MFD (Matrix):

    def __init__ (self, dtm, mannings, cellsize):
        Matrix.__init__(self, dtm)
            
        self.cellsize = cellsize

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
        flood_factor = 0
        visited = dict()
        last_alt = []

        def _drainpaths (rcs, next_step=dict(), queue=list(), level=0):
            try:
                next_level = list()
                for rc in rcs:
                    if rc in visited or not self.mannings[rc] or not self.dtm[rc]: continue                    
                    src_flood = floods[rc]

                    if src_flood / math.pow(self.cellsize, 2.0) < self.cellsize * 0.01:
                        # print("NOT ENOUGHT")
                        ####
                        # La distancia entre el inicio y el punto más bajo alcanzado crece a medida
                        # que aumenta el área calculada. El criterio actual es más exigente al principio
                        # que al final y no distingue entre montaña y plano. La distancia entre el punto
                        # actual y el más bajo alcanzado ha de medirse como excesiva en base a un criterio
                        # que
                        ####
                        if self.dtm[rc] - last_alt[-1] > (self.dtm[start] - last_alt[-1]) * 0.5:
                            visited[rc] = True
                            continue

                        if level == 0:
                            floods[rc] += src_flood * flood_factor * min(1, speeds[rc] / self.cellsize)
                            drafts[rc] = self.get_draft(rc, floods[rc])
                        
                        next_step[rc] = True
                        continue

                    slopes = self.get_slopes(rc, drafts)
                    not_visiteds = [True] * 9  # self.array(list(map(lambda d: tuple(d) not in visited and tuple(d) != rc, rc + self.deltas)))
                    downslopes = self.get_downslopes(slopes, not_visiteds)
                    upslopes = self.get_upslopes(slopes, not_visiteds)
                    under_volume = self.get_volumetries(downslopes, not_visiteds)
                    over_volume = self.get_volumetries(upslopes, not_visiteds)

                    if downslopes.sum() == 0:
                        over_flood = max(0, src_flood - over_volume.min() * 8)
                        drived_flood = 0
                        if over_flood == 0:
                            if level == 0:
                                floods[rc] += src_flood * flood_factor * min(1, speeds[rc] / self.cellsize)
                                drafts[rc] = self.get_draft(rc, floods[rc])
                                speeds[rc] = 0
                            next_step[rc] = True
                            continue
                    else:
                        drived_flood = min(src_flood, under_volume.sum())
                        over_flood = src_flood - drived_flood

                    visited[rc] = True
                    if rc in next_step: del next_step[rc]

                    over_cacthments = self.where(src_flood > over_volume * 8, src_flood - over_volume * 8, 0)
                    overfloods = over_cacthments ** 2 / (over_cacthments ** 2).sum() * over_flood if over_cacthments.sum() else self.zeros(slopes.shape)
                    drivedfloods = downslopes ** 2 / (downslopes ** 2).sum() * drived_flood if downslopes.sum() else self.zeros(slopes.shape)
                    rc_floods = overfloods + drivedfloods
                    rc_speeds = self.get_speeds(slopes, self.get_draft(rc, src_flood), self.mannings[rc], self.log_and(rc_floods > 0, not_visiteds))
                    # print(rc_speeds)
                    # rc_floods = sum(rc_floods) * self.where(rc_floods > 0, rc_floods ** (rc_speeds/self.cellsize), 0) / \
                    #     max(1e-2, sum(self.where(rc_floods > 0, rc_floods ** (rc_speeds/self.cellsize), 0)))

                    reacheds = list()
                    for i, (flood, speed) in enumerate(zip(rc_floods, rc_speeds)):
                        new_rc = tuple(rc + self.deltas[i])
                        if flood == 0 or speed == 0 or not self.mannings[new_rc] or not self.dtm[new_rc]: continue
                        speeds[new_rc] = (speeds[new_rc] or speed + speed) / 2
                        floods[new_rc] += (flood ** 2 / sum(rc_floods ** 2) + speed ** 2 / sum(rc_speeds ** 2)) / 2 * rc_floods.sum() * min(1, speed / self.cellsize)
                        drafts[new_rc] = self.get_draft(new_rc, floods[new_rc])
                        if new_rc not in visited:
                            if speed / self.cellsize > 1:
                                reacheds.append(new_rc)
                            else:
                                next_level.append(new_rc)

                    if self.dtm[rc] < last_alt[-1]:
                        last_alt.append(self.dtm[rc])
                    if len(reacheds) > 0: _drainpaths(reacheds, next_step, queue, level + 1)

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
            last_alt.append(self.dtm[start])
            next_step = {start: True}
            i = 0
            progress = progress_bar(break_time)
            for flood in hyd:
                progress(i)
                flood_factor = (flood / last_flood) if last_flood else 0
                next_step = _drainpaths(
                    next_step,
                    dict()
                )
                last_flood = flood
                i += 1
                
        except KeyboardInterrupt:
            print("KeyboardInterruption!")
        except Exception:
            print("Exception!")
            print_exception()
        finally:
            return floods, drafts, speeds