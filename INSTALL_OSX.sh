#!/bin/sh

# Must have homebrew installed. 
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

# Install some necessaries.
sudo brew install cmake
sudo brew install geos
sudo brew install gdal
sudo brew install qt5
sudo brew link geos
sudo brew link gdal
sudo brew link qt5

# Build the software
mkdir -p build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
make 

# Install
sudo make install
