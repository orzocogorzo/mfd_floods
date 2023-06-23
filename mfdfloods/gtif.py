from osgeo import gdal


def openf (filename):
    return gdal.Open(filename)

def as_array (df, band=1):
    return df.GetRasterBand(band).ReadAsArray()


def get_geotrans (filename=None, ds=None):
    ds = ds or gdal.Open(filename)
    return ds.GetGeoTransform()


def get_rowcol(lng: float, lat: float, ds=None, filename=None):
    ds = ds or openf(filename)
    geotransform = get_geotrans(filename=filename, ds=ds)
    return (int((geotransform[3] - lat) / -geotransform[5]), int((lng - geotransform[0]) / geotransform[1]))


def writef (filename, data, ref_file):
    ref_ds = gdal.Open(ref_file)
    geotransform = ref_ds.GetGeoTransform()
    projection = ref_ds.GetProjection()
    
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(
        filename,
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
    
    print("data saved as %s" % filename)
