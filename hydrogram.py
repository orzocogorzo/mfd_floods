def hydrogram (f):
    f = float(f)
    
    def _range (start, stop, step):
        x = start
        while x <= stop:
             yield x
             x += step
    
    time = _range(0.0, 1.0, 1e-2)
    
    def _hydrogram ():
        for t in time:
            yield f**-t * f + 60

    return _hydrogram()