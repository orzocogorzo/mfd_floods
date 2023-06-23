# MFD Floods

A python script to model hidrologic behavior of downstream drainpaths. With a DTM cover of your study area you can define longitude and latitude (in the DTM distance unit) as a start point and an income flow to see how the water will flood the territory.

## Installation

With pip `pip install mfdfloods`

From source `python -m pip install -r /path/to/mfdfloods/directory`

## Dependencies

The script requires GDAL installed on your system and python-gdal as a python dependency.

To install GDAL execute `apt install gdal-bin libgdal-dev`.

## Usage as command line

MFD Floods requires a two raster layer, in a [GDAL raster format](https://gdal.org/user/raster_data_model.html), like [GTiff](https://gdal.org/drivers/raster/gtiff.html). The first required layer is the DTM covering the area of study. The second, is a raster layer containing manning values of the terrain. Place this two layers on a directory like this:

```
my_area/
├── dtm.tif
├── hydrogram.csv
├── mannings.tif
└── pks.csv
```

The other two required files are the **hidrogram.csv** and the **pks.csv**.

The first file should contain data about the exit hidrogram as two set of values structured as two columns

| Seconds | Cubic meters |
| ------- | ------------ |
| 10      | 3.45         |
| 20      | 6.21         |
| 30      | 7.34         |

The second file should contain pair of coordinates as start points for the modile identified by an ID like this

| id  | lng      | lat       |
| --- | -------- | --------- |
| 1   | 335713.6 | 4706584.6 |
| 2   | 335731.2 | 4706545.4 |
| 3   | 336621.9 | 4704577.0 |

Once you have your are of study ready, then run `python -m mfdfloods ./my_area <start_point_id>` to perform a the hidrological modeling. The script will output three new files:

1. drafts.tif
2. speeds.tif
3. floods.tif

## Usage as python dependency

Include mfdfloods as a module on your scripts with `from mfdfloods import MFD` then instantiate the class MFD to execute its drainpaths method.
