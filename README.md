This is the Wiki page for the HLRG tools, including:

* [contrem](https://github.com/rskelly/contrem/wiki/contrem) -- The convex hull continuum removal tool.
* [convolve](https://github.com/rskelly/contrem/wiki/convolve) -- The spectral convolution tool.
* [reflectance](https://github.com/rskelly/contrem/wiki/reflectance) -- Tool for calculating reflectance from radiance images.
* [refl_regress](https://github.com/rskelly/contrem/wiki/refl_regress) -- Tool for calculating reflectance regression coefficients. 

Production instructions for regression coefficients are [here](https://github.com/rskelly/contrem/wiki/process).

# Installation

First, install Git. If you do not have it, it can be installed on Linux using the package manager, or on OSX by downloading from [here](https://git-scm.com/download/mac).

Next, open a terminal window, navigate to a directory (type 'cd' and then drag the target folder into the window, otherwise just type the full path) and check out the repository. For example:
    
    cd ~/Documents
    git clone https://github.com/rskelly/geotools
    
This will create a folder called geotools in the Documents folder with the source code, etc., inside. Navigate into the folder by typing,

    cd geotools
    
Then follow the instructions below.

## Linux
1) mkdir build
3) cd build
4) cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
5) make
6) sudo make install

Note: The program requires GDAL and GEOS so these must be installed first, preferably using the package manager.

## OSX
1) ./INSTALL_OSX.sh

Note: the INSTALL_OSX.sh file must have execute permissions set to run. Change them using the following command in Terminal:

    chmod +x INSTALL_OSX.sh

Note: The program requires GDAL and GEOS so these must be installed first. The install script attempts this
using homebrew.

## Windows
1) Nope.

