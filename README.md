# geotools

`geotools` is a set of tools for working with LiDAR point clouds and hyperspectral data, developed during my time with the (now-defunct) Hyperspectral-LiDAR Research Group and the Water and Climate Impacts Research at the University of Victoria. The links below go to the Wiki pages for each tool.

* [contrem](https://github.com/rskelly/geotools/wiki/contrem) -- The convex hull continuum removal tool.
* [convolve](https://github.com/rskelly/geotools/wiki/convolve) -- The spectral convolution tool.
* [reflectance](https://github.com/rskelly/geotools/wiki/reflectance) -- Tool for calculating reflectance from radiance images.
* [refl_regress](https://github.com/rskelly/geotools/wiki/refl_regress) -- Tool for calculating reflectance regression coefficients. 
* [pc2grid](https://github.com/rskelly/geotools/wiki/pc2grid) -- Calculates statistics etc. on point clouds of any size.
* [pcnorm](https://github.com/rskelly/geotools/wiki/pcnorm) -- Normalizes point clouds (subtracts the ground elevation from each point).
* [voidfill](https://github.com/rskelly/contrem/wiki/voidfill) -- Fill voids in LiDAR-derived rasters, using "alpha shapes."

Production instructions for regression coefficients are [here](https://github.com/rskelly/geotools/wiki/Reflectance-Coefficients).

# Installation

Install the program by doing the usual:

## Linux
1) Checkout from git
2) mkdir build
3) cd build
4) cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
5) make
6) sudo make install

The program requires GDAL and GEOS so these must be installed first.

## OSX

* Note: The OSX install script isn't really maintained. *

1) ./INSTALL_OSX.sh

The program requires GDAL and GEOS so these must be installed first. The install script attempts this
using homebrew.

## Windows
1) Nope. (There is a Dockerfile in /docker, so that might work. It remains untested except on Linux.)
