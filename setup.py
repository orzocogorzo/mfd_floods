import setuptools

setuptools.setup(
    name="mfdfloods",
    version="0.1",
    scripts=["mfd"],
    author="Orzo Cogorzo",
    author_email="orzocogorzo@hotmail.com",
    description="A python script to modelate hidrologic behavior of downstream drainpaths",
    url="https://github.com/orzocogorzo/mfdfloods",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)