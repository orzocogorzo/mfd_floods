# MFD Floods
A python script to model hidrologic behavior of downstream drainpaths. With a DTM cover of your study area you can define longitude and latitude (in the DTM distance unit) as a start point and an income flow to see how the water will flood the territory.

## Installation
With pip ```pip install mfdfloods```

From source ```python setup.py install```

## Dependencies
The script requires GDAL installed on your system and python-gdal as a python dependency.

To install GDAL execute `apt install gdal-bin libgdal-dev`.

## Test
Execute test.py from inside the folder to test the algorithm.

The test.py is a script that call the class MFD and execute its modelization with the datasource from the `data/` folder. There you have to place your GeoJSON files with the modelized line geometry.

`python test.py <area::string> <X::float> <Y::float> <break_flow::int> <base_flow::int> <break_time::int> [cellsize::int]{5} [radius::int]{2000}`

Arguments:

1. **Area** is the file on the `data/` folder with your DTM in GeoTiff format.
2. **X** is the longitude in your reference dtm distance units.
3. **Y** is the latitude in your reference dtm distance units.
4. **break_flow** is the start income flow.
5. **base_flow** is the base income flow after the initial pic.
6. **break_time** is the time that has to pass in seconds to go from the break_flow to the base_flow.
7. **cellsize** is the size of your DTM cellsize.
8. **radius** is the maximum extension of the output flood in your dtm distance units.

The output will be placed in your `data/` floder as three files with the name *drainpaths_{x}-{y}_(draft|flood|speed).tif*.

## Use

Include mfdfloods as a module on your scripts with `from mfdfloods import MFD` then instantiate the class MFD to execute its drainpaths method.
