# geotools

`geotools` is a set of tools for working with LiDAR point clouds and hyperspectral data, developed during my time with the (now-defunct) Hyperspectral-LiDAR Research Group and the Water and Climate Impacts Research at the University of Victoria. The links below go to the Wiki pages for each tool.

* [contrem](https://github.com/rskelly/geotools/wiki/contrem) -- The convex hull continuum removal tool.
* [convolve](https://github.com/rskelly/geotools/wiki/convolve) -- The spectral convolution tool.
* [reflectance](https://github.com/rskelly/geotools/wiki/reflectance) -- Tool for calculating reflectance from radiance images.
* [refl_regress](https://github.com/rskelly/geotools/wiki/refl_regress) -- Tool for calculating reflectance regression coefficients. 
* [pc2grid](https://github.com/rskelly/geotools/wiki/pc2grid) -- Calculates statistics etc. on point clouds of any size.
* [pcnorm](https://github.com/rskelly/geotools/wiki/pcnorm) -- Normalizes point clouds (subtracts the ground elevation from each point).

Production instructions for regression coefficients are [here](https://github.com/rskelly/geotools/wiki/Reflectance-Coefficients).

## Installation

First, install Git. If you do not have it, it can be installed on Linux using the package manager, or on OSX by downloading from [here](https://git-scm.com/download/mac).

Next, open a terminal window, navigate to a directory (type 'cd' and then drag the target folder into the window, otherwise just type the full path) and check out the repository. For example:
    
    cd ~/Documents
    git clone https://github.com/rskelly/geotools
    
This will create a folder called geotools in the Documents folder with the source code, etc., inside. Navigate into the folder by typing,

    cd geotools
    
Then follow the instructions below.

### Linux
1) `mkdir build`
3) `cd build`
4) `cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..`
5) `make`
6) `sudo make install`

Note: The program requires GDAL and GEOS so these must be installed first, preferably using the package manager.

### OSX

*Note, these instructions are not current, and not likely to be updated soon.*

1) `./INSTALL_OSX.sh`

Note: If the script fails with a message about permissions, it may need to have the execute permissions set to run. Change them using the following command in Terminal:

    chmod +x INSTALL_OSX.sh

Note: The program requires GDAL, GEOS and some other libraries to run so these must be installed first. The install script attempts this using Homebrew. It is a common point of failure.

### Windows
1) Nope.

