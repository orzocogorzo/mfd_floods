# BUILT INS
import sys
import os
import csv
import math
from typing import Optional
from numpy.typing import NDArray

# VENDOR
import richdem as rd
import numpy as np
# import matplotlib.pyplot as plt

# MODULES
from mfdfloods.matrix import Matrix
from mfdfloods.hydrogram import gen_hydrogram
from mfdfloods.gtif import openf, as_array, get_rowcol, writef
from mfdfloods.geotransform import GeoTransformFit
from mfdfloods.debug import print_exception, progress_counter  # , crono, truncate, progress_bar


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
        convergence_factor: float = 2,
        slope_trawl: float = 2,
        mute: bool = True
    ) -> None:
        self.dtm_ds = openf(dtm_path)
        self.dtm_gt = self.dtm_ds.GetGeoTransform()
        self.mannings_ds = openf(mannings_path)
        self.mannings_gt = self.mannings_ds.GetGeoTransform()

        Matrix.__init__(self, as_array(self.dtm_ds))

        self.cellsize = (self.dtm_gt[1] + abs(self.dtm_gt[5])) / 2
        self.cellarea = math.pow(self.cellsize, 2.0)
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
        self.max_drain = 1e+8
        self.mute = mute

    def __del__(self) -> None:
        del self.dtm_ds
        del self.mannings_ds

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
        # Get periferic alt deltas
        if self_draft is None:
            self_draft = float(drafts[rc])

        return self.array([
            (float(self.dtm[tuple(delta)]) + float(drafts[tuple(delta)])) - (self.dtm[rc] + self_draft)
            for delta in self.get_deltas(rc)
        ])

    def get_slope(self, slopes: NDArray) -> float:
        # Max cell traversal slope
        try:
            slopes = np.append(slopes, 0)
            return slopes.min() - slopes.max()
        except Exception:
            return 0

    def get_volumetries(self, slopes: NDArray) -> NDArray:
        # Volumetrie of the pyramide from the center to the edge (half of the cell)
        return self.cellarea * 0.25 * slopes * (1 / 3)

    def get_downslopes(self, slopes: NDArray) -> NDArray:
        # Negativa alt deltas
        return self.where(slopes < 0, slopes * -1, 0)

    def get_upslopes(self, slopes: NDArray) -> NDArray:
        # Positive alt deltas
        return self.where(slopes >= 0, slopes, 0)

    def get_draft(self, flood: float) -> float:
        # return (flood + self.get_volumetries(slopes * .5).sum()) / self.cellarea
        return flood / self.cellarea

    def get_speeds(self, slopes: NDArray, draft: float, manning) -> NDArray:
        return self.array([self.get_speed(draft, manning, slope) for slope in slopes])

    def get_speed(self, draft: float, manning, slope: float) -> float:
        # Manning formula
        return max(0, (1. / manning) * math.pow(self.cellsize + 2. * draft, 2. / 3.) * math.pow(max(.0, abs(slope)) / self.cellsize, .5))

    # @crono
    def drainpaths(self, source: tuple, hydrogram_curve: list) -> tuple[NDArray, NDArray, NDArray]:
        floods = self.zeros(self.dtm.shape)
        drafts = self.zeros(self.dtm.shape)
        speeds = self.zeros(self.dtm.shape)
        drainages = self.zeros(self.dtm.shape)
        flood_factor = 0
        self.is_over = False
        self.overcomes = {}

        def _drainpaths(
            rcs: dict,
        ) -> tuple[dict, dict]:
            next_step = {}
            catchments = {}
            try:
                if self.is_over:
                    return {}, {}

                for rc in rcs:
                    if rc in self.overcomes:
                        continue

                    src_deltas = self.get_deltas(rc)
                    src_flood = floods[rc]  # max(0, float(floods[rc]) + catchments.get(rc, 0))
                    src_draft = self.get_draft(float(src_flood))
                    src_slopes = self.get_slopes(rc, drafts, src_draft)
                    src_slope = self.get_slope(src_slopes)
                    src_speed = self.get_speed(src_draft, self.mannings[rc], src_slope)

                    if src_speed / self.cellsize < 1:
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
                        over_flood = max(0, src_flood - over_volume.min() * 8)
                        drived_flood = 0
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
                    over_floods = self.where(over_floods > 1e-2, over_floods, 0)
                    drived_floods = downslopes / downslopes.sum() * drived_flood if downslopes.sum() else self.zeros((len(src_deltas),))
                    drived_floods = self.where(drived_floods > 1e-2, drived_floods, 0)
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
                    for i, (flood, speed) in enumerate(zip(src_floods, src_speeds)):
                        new_rc = tuple(src_deltas[i])
                        try:
                            if self.mannings[new_rc] == self.nodata or self.dtm[new_rc] == self.nodata:
                                raise IndexError
                        except IndexError:
                            self.is_over = True
                            return {}, {}

                        if flood > 0 and speed > 0:
                            # speed = max(speeds[new_rc], speed)
                            flood = ((flood ** self.convergence_factor / powered_flood + speed ** self.slope_trawl / powered_speed) / 2 * src_acum_flood)
                            # draft = flood / self.cellarea
                            # catchments[new_rc] = catchments.get(new_rc, 0) + flood
                            # catchments[rc] = catchments.get(rc, 0) - flood
                            catchments[new_rc] = flood
                            # catchments[rc] = flood

                            if speed / self.cellsize < 1:
                                if drainages[new_rc] > self.max_drain:
                                    self.overcomes[new_rc] = True
                                    del_key(new_rc, next_step)
                            else:
                                next_step[new_rc] = True
                                self.overcomes[rc] = True


            except KeyboardInterrupt:
                self.is_over = True
                return {}, {}
            except Exception as e:
                raise e
            finally:
                return next_step, catchments

        # DEBUG
        # steps, outs, lens = list(), list(), list()
        # END

        try:
            source = get_rowcol(*source, ds=self.dtm_ds)
            self.overcomes[source] = True
            start, slope = self.start_point(source, drafts)

            hyd = gen_hydrogram(hydrogram_curve)
            break_flood = 0
            while break_flood == 0:
                # while break_flood < self.cellsize:
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
            distance = 0
            trapped = 0
            next_step = {start: True}

            while True:
                # import pdb; pdb.set_trace()
                progress(i, flood)
                next_step, catchments = _drainpaths(next_step)

                try:
                    flood = next(hyd)
                    flood_factor = flood / last_flood
                except (ZeroDivisionError, StopIteration):
                    print("\nExit condition: Hydrogram drained")
                    break

                for rc in catchments:
                    catchment = max(0, catchments[rc]) * flood_factor
                    if catchment <= 0:
                        continue

                    floods[rc] += catchment
                    drafts[rc] = self.get_draft(catchment)
                    slope = self.get_slope(self.get_slopes(rc, drafts))
                    speeds[rc] = self.get_speed(float(drafts[rc]), self.mannings[rc], slope)
                    drainages[rc] += 1

                # DEBUG
                # outs.append(flood)
                # steps.append(len(next_step))
                # lens.append(self.array([math.sqrt(sum(coord**2)) for coord in abs(self.argwhere(floods > 0) - start) * self.cellsize]).max())
                # END

                edge = np.sqrt(np.power(abs(self.argwhere(floods > 0) - start) * self.cellsize, 2).sum(1)).max()
                if distance == int(edge):
                    trapped += 1
                else:
                    trapped = 0

                distance = int(edge)
                i += 1
                if self.is_over:
                    print("\nExit condition: Flood is over dtm boundaries")
                    break
                elif i > 1e+4:
                    print("\nExit condition: Max recursion limit")
                    break
                elif trapped >= 5e+3:
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
            # plt.plot(list(range(0, len(steps), 10)), [math.log(d) if d else d for d in [sum(outs[i:i+10]) for i in range(0, len(steps), 10)]], "b-")
            # plt.plot(list(range(0, len(steps), 10)), [math.log(d) if d else d for d in [sum(steps[i:i+10]) for i in range(0, len(steps), 10)]], "r-")
            # plt.plot(list(range(0, len(steps), 10)), [math.log(d) if d else d for d in [sum(lens[i:i+10]) for i in range(0, len(lens), 10)]], "y--")
            # plt.show()
            return floods, drafts, speeds


def main(
    dtm_path: str,
    mannings_path: str,
    hydrogram: list,
    lng: float,
    lat: float,
) -> None:
    """
    Runs a floods distribution modelation. 

    Parameters:
    area <str>: Name of the area
    lng <float>: Longitude of the brak point
    lat <float>: Latitude of the break point
    hydrogram <list[typle[float, float]]>: A list of pair values with time and flow representing the break hydrogram

    Returns:
    None: The script will write three raster files on the data directory
    """

    floods, drafts, speeds = None, None, None
    try:
        model = MFD(dtm_path, mannings_path, radius=3000, mute=False)
        floods, drafts, speeds = model.drainpaths((lng, lat), hydrogram)
    except KeyboardInterrupt as e:
        print(e)
        print("Keyboard Interruption")
    finally:
        if not (floods is None or drafts is None or speeds is None):
            data_dir = os.path.dirname(dtm_path)
            writef(os.path.join(data_dir, "floods.tif"), floods, dtm_path)
            writef(os.path.join(data_dir, "drafts.tif"), drafts, dtm_path)
            writef(os.path.join(data_dir, "speeds.tif"), speeds, dtm_path)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python -m mfdfloods <path:data_dir> <float:lng> <float:lat> [float:radius]")
        exit()

    kwargs = dict()
    data_dir = os.path.abspath(sys.argv[1])
    kwargs["lng"] = float(sys.argv[2])
    kwargs["lat"] = float(sys.argv[3])

    dtm_path = os.path.join(data_dir, "dtm.tif")
    if not os.path.isfile(dtm_path):
        raise FileNotFoundError(dtm_path + " does not exists")
    else:
        kwargs["dtm_path"] = dtm_path

    mannings_path = os.path.join(data_dir, "mannings.tif")
    if not os.path.isfile(mannings_path):
        raise FileNotFoundError(mannings_path + " does not exists")
    else:
        kwargs["mannings_path"] = mannings_path

    hydrogram_name = os.path.join(data_dir, "hydrogram.csv")
    if not os.path.isfile(hydrogram_name):
        raise FileNotFoundError(hydrogram_name + " does not exists")
    else:
        with open(hydrogram_name) as f:
            reader = csv.reader(f, delimiter=",", quotechar='"')
            kwargs["hydrogram"] = [row for row in reader]

    if len(sys.argv) == 5:
        kwargs["radius"] = float(sys.argv[4])

    main(**kwargs)
