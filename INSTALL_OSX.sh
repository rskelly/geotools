#!/bin/sh

# Must have homebrew installed.

# Install some necessaries.
brew install cmake
brew install geos
brew install gdal
brew link geos
brew link gdal

# Build the software
mkdir -p build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
make 

# Install
sudo make install