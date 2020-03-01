# def hydrogram (f):
#     f = float(f)
    
#     def _range (start, stop, step):
#         x = start
#         while x <= stop:
#              yield x
#              x += step
    
#     time = _range(0.0, 1.0, 1e-2)
    
#     def _hydrogram ():
#         count = 0
#         while True:
#             yield f/count * f + 60
#             count += 1

#     return _hydrogram()


def hydrogram (init, base):
    count = 0
    while True:
        if count == 0: 
            yield init
        else:
            yield count**-1/2 * init + base
        count += 1