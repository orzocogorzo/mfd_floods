class GeoTransformFit:

    def __init__(self, array, geotransform, ref_geotransform, nodata=None):
        self._array = array
        self.x_delta = int((geotransform[0] - ref_geotransform[0]) / abs(geotransform[1]))
        self.y_delta = int((geotransform[3] - ref_geotransform[3]) / abs(geotransform[5]))
        self._nodata = nodata

    def __getitem__(self, rc):
        try:
            return self._array[(rc[0] + self.x_delta, rc[1] + self.y_delta)]
        except IndexError:
            return self._nodata

    def __setitem__(self, rc, value):
        self._array[(rc[0] + self.x_delta, rc[1] + self.y_delta)] = value
