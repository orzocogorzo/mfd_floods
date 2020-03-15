import numpy as np
from osgeo import gdal
from mfd import MFD
import sys
    

def openf (filename):
    return gdal.Open("data/" + filename)


def get_geotrans (filename=None, ds=None):
    ds = ds or gdal.Open("data/" + filename)
    return ds.GetGeoTransform()


def get_rowcol (lng, lat, ds=None, filename=None):
    ds = ds or openf(filename)
    geotransform = get_geotrans(filename=filename, ds=ds)
    return (int((geotransform[3] - lat)/5), int((lng - geotransform[0])/5))


def writef (filename, data, ref_file):
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


def test (area, lng, lat, break_flow, base_flow, break_time):
    try:
        dtm = openf(area + "_dtm.tif")
        manning = openf(area + "_mannings.tif")
        rowcol = get_rowcol(lng, lat, ds=dtm)
        dtm_array = dtm.GetRasterBand(1).ReadAsArray()
        model = MFD(
            dtm=dtm_array,
            mannings=manning.GetRasterBand(1).ReadAsArray(),
            cellsize=5
        )
        floods, drafts, speeds = model.drainpaths(rowcol, break_flow, base_flow, break_time)
    except KeyboardInterrupt as e:
        print("Keyboard Interruption")
    finally:
        writef(area + "_floods_%s-%s.tif" % rowcol, floods, area + "_dtm.tif")
        writef(area + "_drafts_%s-%s.tif" % rowcol, drafts, area + "_dtm.tif")
        writef(area + "_speeds_%s-%s.tif" % rowcol, speeds, area + "_dtm.tif")


if __name__ == "__main__":
    area = str(sys.argv[1])
    lng = float(sys.argv[2])
    lat = float(sys.argv[3])
    break_flow = int(sys.argv[4])
    base_flow = int(sys.argv[5])
    break_time = int(sys.argv[6])
    test(area, lng, lat, break_flow, base_flow, break_time)
