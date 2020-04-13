# SYS
import sys, os
import linecache
from time import time

# VENDOR
import numpy as np

# sys.setrecursionlimit(3000)
np.seterr(all="warn")
np.set_printoptions(precision=2)


def print_exception ():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print('EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj))


def crono (fn):

    def wrapper (*args, **kwargs):
        start = time()
        res = fn(*args, **kwargs)
        wasted = time() - start
        print("[CRONO]: ", wasted)
        return res

    return wrapper

def truncate (val, precission=2):
    return int(10**precission*val)/10**precission


def progress_bar (length):
    length = float(length)
    start = time()

    sys.stdout.write("\n    0% [")
    sys.stdout.write(" " * (50))
    sys.stdout.write("]\r")
    sys.stdout.flush()

    def _progress_bar (index):
        # index = float(index)
        progress = int(index / length * 100)
        sys.stdout.write("  {!s}% [".format(str(progress)))
        sys.stdout.write("#" * (progress))
        sys.stdout.write(" " * (100 - progress))
        sys.stdout.write("] {!s}s  \r".format(str(truncate(time() - start, 0))))
        sys.stdout.flush()

    return _progress_bar


def progress_counter (name):

    sys.stdout.write("\n[ITER]: 0 [SECONDS]: 0s  \r")
    start = time()

    def _progress_iter (iter, val):
        delta = int(time() - start)  # truncate(time() - start, 0)
        sys.stdout.write(f"[ITER]: {iter} [SECONDS]: {delta}s [{name}]: {truncate(val, 4)}       \r")
        sys.stdout.flush()

    return _progress_iter
    
