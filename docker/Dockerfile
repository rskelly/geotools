FROM ubuntu
RUN apt-get update && apt-get install -y wget build-essential libgdal-dev liblas-dev libpdal-dev libcgal-dev git qt5-default libeigen3-dev python3-pip libpython3-dev gfortran ccache liblas-dev liblas-c-dev 
RUN pip3 install --upgrade cmake
RUN git clone --recursive https://github.com/rskelly/geotools && mkdir -p geotools/build
RUN cd geotools/build && cmake .. && make && make install

