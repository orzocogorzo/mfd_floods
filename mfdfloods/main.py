# BULTINS
import math
from numpy.typing import NDArray

# VENDOR
import richdem as rd

# MODULES
from .matrix import Matrix
from .hydrogram import gen_hydrogram
from .gtif import openf, as_array, get_rowcol
from .geotransform import GeoTransformFit
from .debug import print_exception, progress_counter  # , crono, truncate, progress_bar
# from .geotransform import GeoTransformFit


class MFD(Matrix):

    def __init__(self, dtm_path: str, manning_path: str, nodata: float = -999, mute: bool = True):
        self.dtm_ds = openf(dtm_path)
        self.dtm_gt = self.dtm_ds.GetGeoTransform()
        self.mannings_ds = openf(manning_path)
        self.mannings_gt = self.mannings_ds.GetGeoTransform()

        Matrix.__init__(self, as_array(self.dtm_ds))

        self.cellsize = self.dtm_gt[1]
        self.cellarea = math.pow(self.cellsize, 2.0)
        self.nodata = nodata

        self.dtm = rd.rdarray(self.dtm, no_data=nodata)
        rd.FillDepressions(self.dtm, in_place=True)

        self.mannings = self.array(as_array(self.mannings_ds))
        self.mannings = GeoTransformFit(
            self.mannings,
            self.mannings_gt,
            self.dtm_gt
        )

        self.mute = mute

    def __del__(self):
        del self.dtm_ds
        del self.mannings_ds

    def start_point(self, rc, drafts):
        slopes = self.get_slopes(rc, drafts)
        gateway = slopes.argmin()
        return tuple(rc + self.deltas[gateway]), {
            rc: True,
            **{
                tuple(delta): True
                for delta in (rc + self.deltas[gateway] * -1) + self.deltas
            }
        }, slopes[gateway]

    def get_slopes(self, rc, drafts) -> NDArray:
        # Get periferic alt deltas
        return self.array([(self.dtm[tuple(delta)] + drafts[tuple(delta)]) - (self.dtm[rc] + drafts[rc]) for delta in rc + self.deltas])

    def get_volumetries(self, slopes) -> NDArray:
        # Volumetrie of the pyramide from the center to the edge (half of the cell)
        return self.cellarea * 0.5 * slopes * (1 / 3)

    def get_downslopes(self, slopes) -> NDArray:
        # Negativa alt deltas
        return self.where(slopes < 0, slopes * -1, 0)

    def get_upslopes(self, slopes) -> NDArray:
        # Positive alt deltas
        return self.where(slopes >= 0, slopes, 0)

    def get_draft(self, flood) -> float:
        return flood / self.cellarea

    def get_speeds(self, slopes, draft, manning) -> NDArray:
        return self.array(list(map(lambda slope: self.get_speed(draft, manning, slope), slopes)))

    def get_speed(self, draft, manning, slope) -> float:
        return max(1e-3, (1. / manning) * math.pow(self.cellsize + 2 * draft, 2. / 3.) * math.pow(max(0, (-1 * slope)) / 5.0, 0.5))
        # return max(1e-3, (1.0/manning) * math.pow((self.cellsize*draft)/(self.cellsize+2*draft), 2.0/3.0) * math.pow(max(0, (-1*slope))/5.0, 0.5))

    # @crono
    def drainpaths (self, src: tuple[float, float], hydrogram_curve: list[tuple[float, float]]):
        floods = self.zeros(self.dtm.shape)
        drafts = self.zeros(self.dtm.shape)
        speeds = self.zeros(self.dtm.shape)
        slopes = self.zeros(self.dtm.shape)
        drainages = self.zeros(self.dtm.shape)
        flood_factor = 0
        visited = dict()
        self.is_over = False

        def _drainpaths(rcs, next_step=dict(), queue=list(), level=0):
            try:
                if self.is_over: return
                next_level = dict()
                reacheds = dict()

                for rc, src_flood in rcs.items():
                    if rc in visited: continue

                    # src_flood = floods[rc]
                    src_speed = speeds[rc]
                    src_slope = slopes[rc]

                    if src_flood < self.cellsize and src_speed / self.cellsize < .5:
                        if level == 0:
                            # floods[rc] += src_flood * flood_factor * min(1, src_speed / self.cellsize)
                            floods[rc] = src_flood + src_flood * flood_factor * min(1, src_speed / self.cellsize)
                            drafts[rc] = self.get_draft(floods[rc])
                            speeds[rc] = self.get_speed(drafts[rc], self.mannings[rc], src_slope)

                        drainages[rc] += 1
                        next_step[rc] = floods[rc]
                        continue

                    rc_slopes = self.get_slopes(rc, drafts)
                    downslopes = self.get_downslopes(rc_slopes )
                    upslopes = self.get_upslopes(rc_slopes)
                    under_volume = self.get_volumetries(downslopes)
                    over_volume = self.get_volumetries(upslopes)

                    if sum(downslopes) == 0:
                        over_flood = max(0, src_flood - over_volume.min() * 8)
                        drived_flood = 0
                        if over_flood == 0:
                            if level == 0:
                                floods[rc] = src_flood + src_flood * flood_factor * min(1, src_speed / self.cellsize)
                                drafts[rc] = self.get_draft(floods[rc])
                                speeds[rc] = 0

                            drainages[rc] += 1
                            next_step[rc] = True
                            continue
                    else:
                        src_flood = floods[rc]
                        drived_flood = min(src_flood, sum(under_volume))
                        over_flood = src_flood - drived_flood

                    visited[rc] = True
                    if rc in next_step: del next_step[rc]

                    over_cacthments = self.where(src_flood > over_volume * 8, src_flood - over_volume * 8, 0)
                    # CATCHMENT DISTRIBUTION. Powers of catchment and slopes defined as the level of concentration/dispersion
                    # of the floods drived by the slopes.
                    overfloods = over_cacthments ** 1 / sum(over_cacthments ** 1) * over_flood if sum(over_cacthments) else self.zeros((9,))
                    drivedfloods = downslopes ** 1 / sum(downslopes ** 1) * drived_flood if sum(downslopes) else self.zeros((9,))
                    rc_floods = overfloods + drivedfloods
                    rc_speeds = self.get_speeds(rc_slopes, drafts[rc], self.mannings[rc])

                    rc_acum_flood = sum(rc_floods)
                    rc_acum_flood3 = sum(rc_floods ** 1)
                    rc_acum_speed3 = sum(rc_speeds ** 3)
                    for i, (flood, speed) in enumerate(zip(rc_floods, rc_speeds)):
                        new_rc = tuple(rc + self.deltas[i])
                        if not self.mannings[new_rc] or not self.dtm[new_rc] or self.mannings[new_rc] == self.nodata or self.dtm[new_rc] == self.nodata:
                            self.is_over = True
                            return

                        slopes[new_rc] = slopes[new_rc] or rc_slopes[i] + rc_slopes[i] / 2
                        speeds[new_rc] = (speeds[new_rc] or speed + speed) / 2

                        # CATCHMENT ASSIGNATION. Based on a ponderation of flood by the speed and powered as the level of 
                        # concentration/dispersion drived by the speed.
                        floods[new_rc] += (flood ** 1 / rc_acum_flood3 + speed ** 3 / rc_acum_speed3) / 2 * rc_acum_flood
                        drafts[new_rc] = self.get_draft(floods[new_rc])

                        # DRAINAGE: Define the critical level of flood when the terrain can drain all the
                        # water and it's impossible the accumulate flood.
                        drainages[new_rc] += 1
                        if (floods[new_rc] / self.cellsize < 1e-4 and drainages[new_rc] > 10) or floods[new_rc] / self.cellsize < 1e-5:
                            if new_rc in next_level: del next_level[new_rc]
                            if new_rc in next_step: del next_step[new_rc]
                            continue

                        if speed / self.cellsize > 1:
                            reacheds[new_rc] = floods[new_rc]
                        else:
                            next_level[new_rc] = floods[new_rc]

                if len(reacheds) > 0: queue.insert(0, reacheds)
                if len(next_level) > 0: queue.append(next_level)
                if len(queue) > 0: _drainpaths(queue.pop(0), next_step, queue, level + 1)
            except Exception:
                print_exception()
            finally:
                return next_step

        try:
            src = get_rowcol(*src, ds=self.dtm_ds)
            start, visited, slope = self.start_point(src, drafts)
            hyd = gen_hydrogram(hydrogram_curve)
            break_flow = 0
            while break_flow == 0:
                break_flow = next(hyd)
            last_flood = None
            trapped = 0
            floods[start] = break_flow
            drafts[start] = self.get_draft(break_flow)
            speeds[start] = self.get_speed(break_flow / self.cellsize ** 2, self.mannings[start], slope)
            next_step = {start: break_flow}
            i = 0
            # steps, news, outs, lens = list(), list(), list(), list()
            if self.mute is False:
                progress = progress_counter("FLOOD")
            else:
                progress = lambda i, f: f

            while True:
                try:
                    flood = next(hyd)
                except StopIteration:
                    break

                progress(i, flood)
                flood_factor = (flood / last_flood) if last_flood else 0
                floods = floods * max(1, flood_factor)
                last_step = next_step
                next_step = _drainpaths(
                    next_step,
                    dict()
                )
                # outs.append(flood)
                # steps.append(len(next_step))
                # news.append(len(list(filter(lambda k: k not in last_step, next_step))))
                # lens.append(self.array([math.sqrt(sum(coord**2)) for coord in abs(self.argwhere(floods > 0) - start) * self.cellsize]).max())
                if len(list(filter(lambda k: k not in last_step, next_step))):
                    trapped = 0
                else:
                    trapped += 1
                # distance = self.array([math.sqrt(sum(coord**2)) for coord in abs(self.argwhere(floods > 0) - start) * self.cellsize]).max()
                last_flood = flood
                i += 1
                if self.is_over or i > 1e+5 or len(next_step) == 0 or trapped >= 1000: break

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
