# BULTINS
import math
from typing import Optional
from osgeo.gdal import Dataset

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
    max_drain: float = 5e1
    dtm_ds: Dataset
    mannings_ds: Dataset
    water_ds: Dataset | None = None
    dtm_gt: tuple
    mannings_gt: tuple
    water_gt: tuple
    dtm: rd.rdarray
    mannings: GeoTransformFit
    water: GeoTransformFit
    cellsize: float
    cellarea: float
    nodata: float
    radius: float
    convergence_factor: float
    speed_trawl: float
    mute: bool

    def __init__(
        self,
        dtm_path: str,
        mannings_path: str,
        water_path: str | None = None,
        nodata: float = -99,
        radius: float = 2000,
        convergence_factor: float = 1.5,
        speed_trawl: float = 1.5,
        mute: bool = True,
    ) -> None:
        self.dtm_ds = openf(dtm_path)
        self.dtm_gt = self.dtm_ds.GetGeoTransform()
        self.mannings_ds = openf(mannings_path)
        self.mannings_gt = self.mannings_ds.GetGeoTransform()
        if water_path is not None:
            self.water_ds = openf(water_path)
            self.water_gt = self.water_ds.GetGeoTransform()

        Matrix.__init__(self, as_array(self.dtm_ds))

        self.cellsize = (self.dtm_gt[1] + abs(self.dtm_gt[5])) / 2.0
        self.cellarea = math.pow(self.cellsize, 2.0)
        self.nodata = nodata

        self.dtm = rd.rdarray(self.dtm, no_data=nodata)
        rd.FillDepressions(self.dtm, in_place=True)

        self.mannings = GeoTransformFit(
            self.array(as_array(self.mannings_ds)),
            self.mannings_gt,
            self.dtm_gt,
        )

        if self.water_ds is not None:
            self.water = GeoTransformFit(
                self.array(as_array(self.water_ds)),
                self.water_gt,
                self.dtm_gt,
            )
        else:
            self.water = GeoTransformFit(
                self.dtm * 0,
                self.dtm_gt,
                self.dtm_gt,
            )

        self.radius = radius
        self.convergence_factor = convergence_factor
        self.speed_trawl = speed_trawl
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
        return self.array([delta for delta in rc + self.deltas])

    def get_slopes(
        self, rc: tuple, drafts: NDArray, self_draft: Optional[float] = None
    ) -> NDArray:
        # Get peripheric alt deltas
        if self_draft is None:
            self_draft = float(drafts[rc])

        return self.array(
            [
                (float(self.dtm[tuple(delta)]) + float(drafts[tuple(delta)]))
                - (self.dtm[rc] + self_draft)
                for delta in self.get_deltas(rc)
            ]
        )

    def get_slope(self, slopes: NDArray) -> float:
        # Max cell traversal slope
        try:
            slopes = np.append(slopes, 0.0)
            return slopes.min() - slopes.max()
        except Exception:
            return 0.0

    def get_volumetries(self, slopes: NDArray) -> NDArray:
        # Volumetrie of the pyramide from the center to the edge (half of the cell)
        return self.cellarea * 0.25 * slopes * (1.0 / 3.0)

    def get_downslopes(self, slopes: NDArray) -> NDArray:
        # Negativa alt deltas
        return self.where(slopes < 0.0, slopes * -1, 0.0)

    def get_upslopes(self, slopes: NDArray) -> NDArray:
        # Positive alt deltas
        return self.where(slopes >= 0.0, slopes, 0.0)

    def get_draft(self, flood: float) -> float:
        # return (flood + self.get_volumetries(slopes * .5).sum()) / self.cellarea
        return flood / self.cellarea

    def get_speeds(self, slopes: NDArray, draft: float, manning) -> NDArray:
        return self.array([self.get_speed(draft, manning, slope) for slope in slopes])

    def get_speed(self, draft: float, manning, slope: float) -> float:
        # Manning formula
        return max(
            0.0,
            (1.0 / manning)
            * math.pow(self.cellsize + 2.0 * draft, 2.0 / 3.0)
            * math.pow(max(0.0, abs(slope)) / self.cellsize, 0.5),
        )

    def drainpaths(
        self, source: tuple, hydrogram_curve: list
    ) -> tuple[NDArray, NDArray, NDArray]:
        floods = self.zeros(self.dtm.shape)
        drafts = self.zeros(self.dtm.shape)
        speeds = self.zeros(self.dtm.shape)
        drainages = self.zeros(self.dtm.shape)
        flood_factor = 0.0
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
                    if self.water[rc] == 1:
                        self.overcomes[rc] = True
                        del_key(rc, next_step)
                        del_key(rc, reacheds)
                        del_key(rc, visited)

                    if self.overcomes.get(rc) is not None:
                        continue

                    if type(rcs[rc]) is int and rcs[rc] > 1:
                        next_step[rc] = rcs[rc] - 1
                        continue

                    src_deltas = self.get_deltas(rc)
                    src_flood = max(0.0, float(floods[rc]) + catchments.get(rc, 0.0))
                    src_draft = self.get_draft(src_flood)
                    src_slopes = self.get_slopes(rc, drafts, src_draft)
                    src_slope = self.get_slope(src_slopes)
                    src_speed = self.get_speed(src_draft, self.mannings[rc], src_slope)

                    if src_speed / level / self.cellsize < 1:
                        if drainages[rc] <= self.max_drain:
                            next_step[rc] = True
                        else:
                            self.overcomes[rc] = True
                            del_key(rc, next_step)
                            del_key(rc, reacheds)
                            del_key(rc, visited)
                        continue

                    downslopes = self.get_downslopes(src_slopes)
                    upslopes = self.get_upslopes(src_slopes)
                    under_volume = self.get_volumetries(downslopes)
                    over_volume = self.get_volumetries(upslopes)

                    if downslopes.sum() == 0:
                        over_flood = max(0.0, src_flood - over_volume.min() * 8)
                        drived_flood = 0.0
                        if over_flood == 0:
                            if drainages[rc] <= self.max_drain:
                                next_step[rc] = True
                            else:
                                self.overcomes[rc] = True
                                del_key(rc, next_step)
                                del_key(rc, visited)
                                del_key(rc, reacheds)
                            continue
                    else:
                        drived_flood = min(src_flood, under_volume.sum())
                        over_flood = src_flood - drived_flood

                    over_catchments = self.where(
                        src_flood > over_volume * 8, src_flood - over_volume * 8, 0
                    )
                    over_floods = (
                        over_catchments / over_catchments.sum() * over_flood
                        if over_catchments.sum()
                        else self.zeros((len(src_deltas),))
                    )
                    over_floods = self.where(over_floods > 1e-3, over_floods, 0)
                    drived_floods = (
                        downslopes / downslopes.sum() * drived_flood
                        if downslopes.sum()
                        else self.zeros((len(src_deltas),))
                    )
                    drived_floods = self.where(drived_floods > 1e-3, drived_floods, 0)
                    src_floods = over_floods + drived_floods
                    src_speeds = self.get_speeds(
                        downslopes, float(drafts[rc]), self.mannings[rc]
                    )

                    if src_floods.sum() == 0:
                        if drainages[rc] <= self.max_drain:
                            next_step[rc] = True
                        else:
                            self.overcomes[rc] = True
                            del_key(rc, next_step)
                            del_key(rc, reacheds)
                            del_key(rc, visited)
                        continue

                    src_acum_flood = src_floods.sum()
                    powered_flood = (
                        np.power(src_floods, self.convergence_factor)
                    ).sum()
                    powered_speed = (np.power(src_speeds, self.speed_trawl)).sum()
                    cell_reacheds = 0
                    for i, (flood, speed) in enumerate(zip(src_floods, src_speeds)):
                        new_rc = tuple(src_deltas[i])
                        try:
                            if (
                                self.mannings[new_rc] == self.nodata
                                or self.dtm[new_rc] == self.nodata
                            ):
                                raise IndexError
                        except IndexError:
                            self.is_over = True
                            return {}, {}

                        if flood > 0 and speed > 0:
                            speed = max(speeds[new_rc], speed)
                            flood = (
                                (
                                    np.power(flood, self.convergence_factor)
                                    / powered_flood
                                    + np.power(speed, self.speed_trawl) / powered_speed
                                )
                                / 2
                                * src_acum_flood
                            ) / level
                            flood = flood * min(1, speed / level / self.cellsize)
                            catchments[new_rc] = catchments.get(new_rc, 0.0) + flood
                            catchments[rc] = catchments.get(rc, 0.0) - flood

                            if visited.get(rc) is not None:
                                continue

                            if (iters := speed / level / self.cellsize) < 1.0:
                                if (
                                    drainages[new_rc] <= self.max_drain
                                    and reacheds.get(new_rc) is None
                                ):
                                    next_step[new_rc] = round(1 / iters)
                                continue

                            reacheds[new_rc] = True
                            cell_reacheds += 1

                    if cell_reacheds == 0:
                        next_step[rc] = True
                        del_key(rc, visited)
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
            distance = 0.0
            trapped = 0
            peak = 0
            next_step = {start: True}
            catchments = {}

            while True:
                progress(i, flood)
                prev_catchments = catchments
                next_step, catchments = _drainpaths(next_step, {}, {}, 1, [], {})

                try:
                    flood = next(hyd)
                    peak = max(flood, peak)
                    flood_factor = flood / last_flood
                except (ZeroDivisionError, StopIteration):
                    print("\nExit condition: Hydrogram drained")
                    break

                # for rc in catchments:
                for rc in next_step:
                    catchments[rc] = catchments.get(rc)
                    if catchments[rc] is None:
                        if self.overcomes.get(rc) is not None:
                            del_key(rc, next_step)
                            continue

                        catchments[rc] = prev_catchments.get(rc, 0)
                    catchment = catchments[rc] * flood_factor
                    if catchment <= 0:
                        continue

                    floods[rc] += catchment
                    drafts[rc] = self.get_draft(catchment)
                    slope = self.get_slope(self.get_slopes(rc, drafts))
                    speeds[rc] = self.get_speed(
                        float(drafts[rc]), self.mannings[rc], slope
                    )
                    drainages[rc] += 1

                prev_catchments = catchments

                edge = np.sqrt(
                    np.power(
                        abs(self.argwhere(floods > 0) - start) * self.cellsize, 2.0
                    ).sum(1)
                ).max()
                if distance == int(edge) and peak == hyd_statistics["peak"]:
                    trapped += 1
                else:
                    trapped = 0

                distance = int(edge)
                i += 1
                if self.is_over:
                    print("\nExit condition: Flood is over dtm boundaries")
                    break
                elif i > 5e3:
                    print("\nExit condition: Max recursion limit")
                    break
                elif trapped >= 2e2:
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
