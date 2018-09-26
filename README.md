This is the Wiki page for the HLRG tools, including:

* [contrem](https://github.com/rskelly/contrem/wiki/contrem) -- The convex hull continuum removal tool.
* [convolve](https://github.com/rskelly/contrem/wiki/convolve) -- The spectral convolution tool.
* [reflectance](https://github.com/rskelly/contrem/wiki/reflectance) -- Tool for calculating reflectance from radiance images.

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
1) ./INSTALL_OSX.sh

The program requires GDAL and GEOS so these must be installed first. The install script attempts this
using homebrew.

## Windows
1) Nope.

