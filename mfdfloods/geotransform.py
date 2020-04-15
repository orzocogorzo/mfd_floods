class GeoTransformFit (object):

    def __init__ (self, array, geotransform, ref_geotransform):
        self._array = array
        self.x_delta = (geotransform[0] - ref_geotransform[0])/geotransform[1]
        self.y_delta = (geotransform[3] - ref_geotransform[3])/(geotransform[5] * -1)
        self.bookmark = 0

    def __getitem__(self, key):
        return self._array[(int(round(key[0] + self.x_delta)), int(round(key[1] + self.y_delta)))]

    def __setitem__(self, key, value):
        self.array[(int(round(key[0] + self.x_delta)), int(round(key[1] + self.y_delta)))] = value