#!/bin/sh

# Must have homebrew installed. 
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

# Install some necessaries.
brew install cmake
brew install geos
brew install gdal
brew install qt5
brew link geos
brew link gdal
brew link qt5

# Build the software
mkdir -p build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
make 

# Install
sudo make install
