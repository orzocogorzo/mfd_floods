# SYS
import sys
import linecache

# VENDOR
import numpy as np

# sys.setrecursionlimit(3000)
# np.seterr(all="raise")
np.set_printoptions(precision=2)


def print_exception ():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))