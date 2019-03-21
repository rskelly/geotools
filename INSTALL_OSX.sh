#!/bin/sh

export PATH="/usr/local/opt/qt/bin:$PATH"
sudo chmod +t /private/tmp
sudo chown -R "$(whoami)":admin /usr/local

# Must have homebrew installed. 
if ! [ -x "$(command -v brew)" ]; then
	/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)";
fi;

# Install some necessaries.
brew update && brew upgrade

brew install cmake
brew install geos
brew install gdal
brew install qt5
brew install eigen
brew install ccache
brew upgrade cmake
brew upgrade geos
brew upgrade gdal
brew upgrade qt5
brew upgrade eigen
brew upgrade ccache
brew unlink cmake && brew link cmake --force
brew unlink geos && brew link geos --force
brew unlink gdal && brew link gdal --force
brew unlink qt && brew link qt5 --force
brew unlink eigen && brew link eigen --force
brew unlink ccache && brew link ccache --force

# Build the software
mkdir -p build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
make 

# Install
sudo make install
