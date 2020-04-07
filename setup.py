import os
from setuptools import setup

def read (fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name="mfdfloods",
    version="0.1",
    scripts=["mfd/__init__.py"],
    author="Orzo Cogorzo",
    author_email="orzocogorzo@hotmail.com",
    description="A python script to modelate hidrologic behavior of downstream drainpaths",
    long_description=read("README.md"),
    license="MIT",
    keywords="mfd drainpaths hidrology downstreams multiple flow direction",
    url="https://github.com/orzocogorzo/mfdfloods",
    requires=[
        "numpy",
        "richdem",
        "GDAL"
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Hydrology"
    ],
)
