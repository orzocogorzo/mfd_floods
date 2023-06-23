# BUILT INS
import sys
import os.path
import csv

# MODULES
from .main import MFD
from .gtif import writef


def main(
    dtm_path: str,
    mannings_path: str,
    water_path: str | None,
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
        model = MFD(dtm_path, mannings_path, water_path, radius=3000, mute=False)
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
    if len(sys.argv) < 3:
        print("Usage: python -m mfdfloods <path:data_dir> <int:pk>")
        exit(1)

    kwargs = dict()
    data_dir = os.path.abspath(sys.argv[1])
    try:
        pk = sys.argv[2]
    except Exception as e:
        print("Invalid pk")
        exit(1)

    data = None
    with open(os.path.join(data_dir, "pks.csv")) as f:
        reader = csv.reader(f)
        headers = None
        for row in reader:
            if data is not None:
                break

            if headers is None:
                headers = row
                continue

            if row[headers.index("pk")] == pk:
                data = {headers[i]: row[i] for i in range(len(headers))}

    if data is None:
        print("Pk not found")
        exit(1)

    dtm_path = os.path.join(data_dir, "dtm.tif")
    if not os.path.isfile(dtm_path):
        raise FileNotFoundError(dtm_path + " does not exists")
    else:
        data["dtm_path"] = dtm_path

    mannings_path = os.path.join(data_dir, "mannings.tif")
    if not os.path.isfile(mannings_path):
        raise FileNotFoundError(mannings_path + " does not exists")
    else:
        data["mannings_path"] = mannings_path

    water_path = os.path.join(data_dir, "water.tif")
    if os.path.isfile(water_path):
        data["water_path"] = water_path
    else:
        data["water_path"] = None

    hydrogram_name = os.path.join(data_dir, "hydrogram.csv")
    if not os.path.isfile(hydrogram_name):
        raise FileNotFoundError(hydrogram_name + " does not exists")
    else:
        with open(hydrogram_name) as f:
            reader = csv.reader(f, delimiter=",", quotechar='"')
            data["hydrogram"] = [row for row in reader]

    if len(sys.argv) == 4:
        data["radius"] = float(sys.argv[3])

    main(
        data["dtm_path"],
        data["mannings_path"],
        data["water_path"],
        data["hydrogram"],
        float(data["lng"]),
        float(data["lat"]),
    )
