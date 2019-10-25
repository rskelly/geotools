#!/usr/bin/env python3

'''
Compute the histogram for the given raster with the given number of buckets.
'''

import os
import sys
import gdal
import numpy as np

try:
    infile = sys.argv[1]
    bins = int(sys.argv[2])
    band = int(sys.argv[3])
except:
    print('Usage: histo.py <raster file> <buckets> <band>')
    sys.exit(1)

ds = gdal.Open(infile)
cols = ds.RasterXSize
rows = ds.RasterYSize
band = ds.GetRasterBand(band)
nodata = band.GetNoDataValue() 

values = band.ReadAsArray()
values = values[values > 0.]
values = np.log(values[values != nodata]


hist = np.histogram(values, bins = bins)

print(hist)
        