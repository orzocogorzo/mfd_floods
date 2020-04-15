import os
from setuptools import setup, find_packages

def read (fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name="MFD Floods",
    version="0.1.12",
    author="Orzo Cogorzo",
    author_email="orzocogorzo@hotmail.com",
    description="A python script to modelate hidrologic behavior of downstream drainpaths",
    long_description=read("README.md"),
    long_description_content_type="text/markdown",
    license="MIT",
    keywords="mfdfloods mfd flood drainpaths hidrology downstreams multiple flow direction",
    url="https://github.com/orzocogorzo/mfdfloods",
    packages=find_packages(),
    install_requires=[
        "numpy==1.18.2",
        "richdem==0.3.4",
        "GDAL==3.0.4"
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
