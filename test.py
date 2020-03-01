import numpy as np
from osgeo import gdal
from mfd import MFD
import sys


def gen_model (data, pxsize):
    return MFD(data, pxsize)
    

# def open_file (filename):
#     ds = gdal.Open("data/" + filename)
#     band = ds.GetRasterBand(1)
    # return band.ReadAsArray()

    
def geotransform (filename):
    ds = gdal.Open("data/" + filename)
    return ds.GetGeoTransform()


def write_file (filename, data, ref_file):
    ref_ds = gdal.Open("data/"+ref_file)
    geotransform = ref_ds.GetGeoTransform()
    projection = ref_ds.GetProjection()
    
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(
        "data/"+filename,
        data.shape[1],
        data.shape[0],
        1,
        gdal.GDT_Float32
    )
    ds.SetGeoTransform(geotransform)
    band = ds.GetRasterBand(1)
    band.WriteArray(data)
    band.SetNoDataValue(0.0)
    ds.SetProjection(projection)
    band.FlushCache()
    
    print("data saved as %s" % "data/"+filename)


def test (lng, lat, maxflow, flow):
    ds = gdal.Open("data/dtm.tif")
    geotransform = ds.GetGeoTransform()
    row = int((geotransform[3] - lat)/5)
    col = int((lng - geotransform[0])/5)
    band = ds.GetRasterBand(1)
    dtm = band.ReadAsArray()
    model = gen_model(dtm, 5)
    drainpath = model.drainpaths((row, col), maxflow, flow)
    write_file("drainpath_%s-%s.tif" % (lat, lng), drainpath, "dtm.tif")


if __name__ == "__main__":
    lng = float(sys.argv[1])
    lat = float(sys.argv[2])
    maxflow = int(sys.argv[3])
    flow = int(sys.argv[4])
    test(lng, lat, maxflow, flow)
