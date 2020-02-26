from osgeo import gdal
from mfd import MFD
import sys


def gen_model (data, pxsize):
    return MFD(data, pxsize)
    

def open_file (filename):
    ds = gdal.Open("data/" + filename)
    band = ds.GetRasterBand(1)
    return band.ReadAsArray()
    
    
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
    ds.SetProjection(projection)
    band.FlushCache()
    
    print("data saved as %s" % "data/"+filename)
    

def test (width, height, flow):
    dtm = open_file("dtm.tif")
    model = gen_model(dtm, 5)
    drainpath = model.drainpaths((width, height), flow)
    write_file("drainpath_%s-%s.tif" % (width, height), drainpath, "dtm.tif")


if __name__ == "__main__":
    width = int(sys.argv[1])
    height = int(sys.argv[2])
    flow = int(sys.argv[3])
    test(width, height, flow)
