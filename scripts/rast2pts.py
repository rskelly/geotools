#!/usr/bin/env python3

import os
import sys
import gdal

filename = sys.argv[1]
band = int(sys.argv[2])

ds = gdal.Open(filename)
data = ds.GetRasterBand(band).ReadAsArray()
x, rx, _, y, _, ry = ds.GetGeoTransform()
nd = ds.GetRasterBand(band).GetNoDataValue()

print('x,y,z')
for r in range(ds.RasterYSize):
    for c in range(ds.RasterXSize):
        if data[r,c] != nd:
            row = list(map(str, [c * rx + x, r * ry + y, data[r,c]]))
            print(','.join(row))

