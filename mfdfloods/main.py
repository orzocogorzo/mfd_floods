# BULTINS
import math
from typing import Optional

# VENDOR
import richdem as rd
import numpy as np
from numpy.typing import NDArray

# MODULES
from .matrix import Matrix
from .hydrogram import gen_hydrogram, hydrogram_statistics
from .gtif import openf, as_array, get_rowcol
from .geotransform import GeoTransformFit
from .debug import print_exception, progress_counter


def del_key(key, handle) -> None:
    try:
        del handle[key]
    except KeyError:
        pass

class MFD(Matrix):
    def __init__(
        self,
        dtm_path: str,
        mannings_path: str,
        nodata: float = -99,
        radius: float = 2000,
        convergence_factor: float = 2.,
        slope_trawl: float = 2.,
        mute: bool = True
    ) -> None:
        self.dtm_ds = openf(dtm_path)
        self.dtm_gt = self.dtm_ds.GetGeoTransform()
        self.mannings_ds = openf(mannings_path)
        self.mannings_gt = self.mannings_ds.GetGeoTransform()

        Matrix.__init__(self, as_array(self.dtm_ds))

        self.cellsize = (self.dtm_gt[1] + abs(self.dtm_gt[5])) / 2.
        self.cellarea = math.pow(self.cellsize, 2.)
        self.nodata = nodata

        self.dtm = rd.rdarray(self.dtm, no_data=nodata)
        rd.FillDepressions(self.dtm, in_place=True)

        self.mannings = self.array(as_array(self.mannings_ds))
        self.mannings = GeoTransformFit(
            self.mannings,
            self.mannings_gt,
            self.dtm_gt,
        )

        self.radius = radius
        self.convergence_factor = convergence_factor
        self.slope_trawl = slope_trawl
        self.max_drain = 1e+1
        self.mute = mute

    def __del__(self) -> None:
        try:
            del self.dtm_ds
            del self.mannings_ds
        except AttributeError:
            pass

    def start_point(self, rc: tuple, drafts: NDArray) -> tuple:
        slopes = self.get_slopes(rc, drafts)
        direction = slopes.argmin()
        deltas = self.get_deltas(rc)
        gateway = deltas[direction]
        
        self.overcomes[rc] = True
        for delta in deltas:
            if not (delta[0] == gateway[0] and delta[1] == gateway[1]):
                self.overcomes[tuple(delta)] = True

        return tuple(gateway), slopes[direction]

    def get_deltas(self, rc: tuple) -> NDArray:
        # Non visited deltas
        return self.array([
            delta for delta in rc + self.deltas
            # if tuple(delta) not in self.overcomes
        ])

    def get_slopes(self, rc: tuple, drafts: NDArray, self_draft: Optional[float] = None) -> NDArray:
        # Get peripheric alt deltas
        if self_draft is None:
            self_draft = float(drafts[rc])

        return self.array([
            (float(self.dtm[tuple(delta)]) + float(drafts[tuple(delta)])) - (self.dtm[rc] + self_draft)
            for delta in self.get_deltas(rc)
        ])

    def get_slope(self, slopes: NDArray) -> float:
        # Max cell traversal slope
        try:
            slopes = np.append(slopes, .0)
            return slopes.min() - slopes.max()
        except Exception:
            return .0

    def get_volumetries(self, slopes: NDArray) -> NDArray:
        # Volumetrie of the pyramide from the center to the edge (half of the cell)
        return self.cellarea * .25 * slopes * (1. / 3.)

    def get_downslopes(self, slopes: NDArray) -> NDArray:
        # Negativa alt deltas
        return self.where(slopes < .0, slopes * -1, .0)

    def get_upslopes(self, slopes: NDArray) -> NDArray:
        # Positive alt deltas
        return self.where(slopes >= .0, slopes, .0)

    def get_draft(self, flood: float) -> float:
        # return (flood + self.get_volumetries(slopes * .5).sum()) / self.cellarea
        return flood / self.cellarea

    def get_speeds(self, slopes: NDArray, draft: float, manning) -> NDArray:
        return self.array([self.get_speed(draft, manning, slope) for slope in slopes])

    def get_speed(self, draft: float, manning, slope: float) -> float:
        # Manning formula
        return max(.0, (1. / manning) * math.pow(self.cellsize + 2. * draft, 2. / 3.) * math.pow(max(.0, abs(slope)) / self.cellsize, .5))

    def drainpaths(self, source: tuple, hydrogram_curve: list) -> tuple[NDArray, NDArray, NDArray]:
        floods = self.zeros(self.dtm.shape)
        drafts = self.zeros(self.dtm.shape)
        speeds = self.zeros(self.dtm.shape)
        drainages = self.zeros(self.dtm.shape)
        flood_factor = .0
        self.is_over = False
        self.overcomes = {}

        def _drainpaths(
            rcs: dict,
            next_step: dict,
            catchments: dict,
            level: int = 1,
            queue: list = [],
            visited: dict = {},
        ) -> tuple[dict, dict]:
            try:
                reacheds = {}
                if self.is_over:
                    return {}, {}

                for rc in rcs:
                    if rc in self.overcomes:
                        continue

                    src_deltas = self.get_deltas(rc)
                    src_flood = max(.0, float(floods[rc]) + catchments.get(rc, .0))
                    src_draft = self.get_draft(src_flood)
                    src_slopes = self.get_slopes(rc, drafts, src_draft)
                    src_slope = self.get_slope(src_slopes)
                    src_speed = self.get_speed(src_draft, self.mannings[rc], src_slope)

                    if src_speed / level / self.cellsize < 1.:
                        if drainages[rc] <= self.max_drain:
                            next_step[rc] = True
                        else:
                            self.overcomes[rc] = True
                            del_key(rc, next_step)
                        continue

                    downslopes = self.get_downslopes(src_slopes)
                    upslopes = self.get_upslopes(src_slopes)
                    under_volume = self.get_volumetries(downslopes)
                    over_volume = self.get_volumetries(upslopes)

                    if downslopes.sum() == 0:
                        over_flood = max(.0, src_flood - over_volume.min() * 8)
                        drived_flood = .0
                        if over_flood == 0:
                            if drainages[rc] <= self.max_drain:
                                next_step[rc] = True
                            else:
                                self.overcomes[rc] = True
                                del_key(rc, next_step)
                            continue
                    else:
                        drived_flood = min(src_flood, under_volume.sum())
                        over_flood = src_flood - drived_flood

                    over_catchments = self.where(src_flood > over_volume * 8, src_flood - over_volume * 8, 0)
                    over_floods = over_catchments / over_catchments.sum() * over_flood if over_catchments.sum() else self.zeros((len(src_deltas),))
                    over_floods = self.where(over_floods > 1e-3, over_floods, 0)
                    drived_floods = downslopes / downslopes.sum() * drived_flood if downslopes.sum() else self.zeros((len(src_deltas),))
                    drived_floods = self.where(drived_floods > 1e-3, drived_floods, 0)
                    src_floods = over_floods + drived_floods
                    src_speeds = self.get_speeds(downslopes, float(drafts[rc]), self.mannings[rc])

                    if src_floods.sum() == 0:
                        if drainages[rc] <= self.max_drain:
                            next_step[rc] = True
                        else:
                            self.overcomes[rc] = True
                            del_key(rc, next_step)
                        continue
                    
                    src_acum_flood = src_floods.sum()
                    powered_flood = (src_floods ** self.convergence_factor).sum()
                    powered_speed = (src_speeds ** self.slope_trawl).sum()
                    cell_reacheds = 0
                    for i, (flood, speed) in enumerate(zip(src_floods, src_speeds)):
                        new_rc = tuple(src_deltas[i])
                        try:
                            if self.mannings[new_rc] == self.nodata or self.dtm[new_rc] == self.nodata:
                                raise IndexError
                        except IndexError:
                            self.is_over = True
                            return {}, {}

                        if flood > 0 and speed > 0:
                            speed = max(speeds[new_rc], speed)
                            flood = ((flood ** self.convergence_factor / powered_flood + speed ** self.slope_trawl / powered_speed) / 2 * src_acum_flood) / level # * (speed / level / self.cellsize)
                            catchments[new_rc] = catchments.get(new_rc, .0) + flood
                            catchments[rc] = catchments.get(rc, .0) - flood

                            if new_rc in visited:
                                continue

                            if speed / level / self.cellsize < 1.:
                                if drainages[new_rc] <= self.max_drain and not new_rc in reacheds:
                                    next_step[new_rc] = True
                                continue

                            reacheds[new_rc] = True
                            cell_reacheds += 1

                    if cell_reacheds == 0:
                        next_step[rc] = True
                    else:
                        visited[rc] = True

                if len(reacheds) > 0:
                    queue.append((reacheds, level + 1))

                if len(queue) > 0:
                    reacheds, level = queue.pop()
                    next_step, catchments = _drainpaths(
                        reacheds,
                        next_step={**next_step, **reacheds},
                        catchments=catchments,
                        level=level,
                        queue=queue,
                        visited=visited,
                    )
                elif len(visited) > 0:
                    for rc in visited:
                        self.overcomes[rc] = True

            except KeyboardInterrupt:
                self.is_over = True
                return {}, {}
            except Exception as e:
                raise e
            finally:
                return next_step, catchments

        try:
            source = get_rowcol(*source, ds=self.dtm_ds)
            self.overcomes[source] = True
            start, slope = self.start_point(source, drafts)

            hyd = gen_hydrogram(hydrogram_curve)
            hyd_statistics = hydrogram_statistics(hydrogram_curve)
            break_flood = 0
            while break_flood == 0:
                break_flood = next(hyd)
                floods[start] = break_flood
                drafts[start] = self.get_draft(break_flood)
                speeds[start] = self.get_speed(
                    float(drafts[start]),
                    self.mannings[start],
                    self.get_slope(self.get_slopes(start, drafts)),
                )

            if self.mute is False:
                progress = progress_counter("FLOOD")
            else:
                progress = lambda i, f: f

            i = 0
            last_flood = break_flood
            flood = break_flood
            distance = .0
            trapped = 0
            peak = 0
            next_step = {start: True}

            while True:
                progress(i, flood)
                next_step, catchments = _drainpaths(next_step, {}, {}, 1, [], {})

                try:
                    flood = next(hyd)
                    peak = max(flood, peak)
                    flood_factor = flood / last_flood
                except (ZeroDivisionError, StopIteration):
                    print("\nExit condition: Hydrogram drained")
                    break

                for rc in catchments:
                    catchment = max(.0, catchments[rc]) * flood_factor
                    if catchment <= 0:
                        continue

                    floods[rc] += catchment
                    drafts[rc] = self.get_draft(catchment)
                    slope = self.get_slope(self.get_slopes(rc, drafts))
                    speeds[rc] = self.get_speed(float(drafts[rc]), self.mannings[rc], slope)
                    drainages[rc] += 1

                edge = np.sqrt(np.power(abs(self.argwhere(floods > 0) - start) * self.cellsize, 2.).sum(1)).max()
                if distance == int(edge) and peak == hyd_statistics["peak"]:
                    trapped += 1
                else:
                    trapped = 0

                distance = int(edge)
                i += 1
                if self.is_over:
                    print("\nExit condition: Flood is over dtm boundaries")
                    break
                elif i > 5e+3:
                    print("\nExit condition: Max recursion limit")
                    break
                elif trapped >= 2e+2:
                    print("\nExit condition: Flood's stability reached")
                    break
                elif distance >= self.radius:
                    print("\nExit condition: Distance limit reached")
                    break

                last_flood = flood

        except KeyboardInterrupt:
            self.is_over = True
            print("KeyboardInterruption!")
        except Exception:
            print_exception()
        finally:
            return floods, drafts, speeds
