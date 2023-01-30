# MFD Floods

A python script to model hidrologic behavior of downstream drainpaths. With a DTM cover of your study area you can define longitude and latitude (in the DTM distance unit) as a start point and an income flow to see how the water will flood the territory.

## Installation

With pip `pip install mfdfloods`

From source `python -m pip install -r /path/to/mfdfloods/directory`

## Dependencies

The script requires GDAL installed on your system and python-gdal as a python dependency.

To install GDAL execute `apt install gdal-bin libgdal-dev`.

## Module

Execute .py from inside the folder to test the algorithm.

The test.py is a script that call the class MFD and execute its modelization with
the datasource from the `data/` folder. There you have to place your GeoJSON files
with the modelized line geometry.

`python test.py <path:data> <lng:float> <lat:float>`

Arguments:

1. **data** is the pathname of the data directory. There you have to place a
   `dtm.tif`, a `mannings.tif` and a `hydrogram.csv` files.
2. **lng** is the longitude in your reference dtm distance units.
3. **lat** is the latitude in your reference dtm distance units.

The output will be placed in your **data** floder as three files with the name _(draft|flood|speed)\_{lng}-{lat}.tif_.

## Use

Include mfdfloods as a module on your scripts with `from mfdfloods import MFD` then instantiate the class MFD to execute its drainpaths method.
