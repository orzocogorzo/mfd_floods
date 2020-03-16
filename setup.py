from distutils.core import setup
from Cython.Build import cythonize

setup(
  ext_modules=cythonize([
    "mfd.py",
    "matrix.py",
    "hydrogram.py"
 ], annotate=True),
)
