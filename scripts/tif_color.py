#!/usr/bin/env python3

'''
Colour a map using its range and a selected scheme.
'''

import gdal
import sys
import os
import numpy as np
import math

def colour(infile, scheme, band, outfile):
    
    ds = gdal.Open(infile)
    band = ds.GetRasterBand(band)
    nd = band.GetNoDataValue()
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    vmin = sys.float_info.max
    vmax = -min 
    for r in range(rows):
        buf = band.ReadAsArray(0, r, cols, 1)
        vmin = min(vmin, np.min(buf))
        vmax = max(vmax, np.max(buf))
        
    
    
if __name__ == '__main__':
    
    infile = sys.argv[1]
    scheme = sys.argv[2]
    band = int(sys.argv[3])
    outfile = sys.argv[4]
    
    colour(infile, scheme, band, outfile)
