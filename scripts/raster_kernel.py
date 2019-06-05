#!/usr/bin/env python3

'''
Call this program from the command line using the syntax:

./raster_kernel.py <input file> <output file>

All bands in the input will be processed using a 9-element moving average
kernel and written to the output file.

You must have scipy and numpy installed. You can do this on Mac and Linux
using:

pip3 install scipy

That should install scipy and numpy which is a dependency.
'''

import sys
import gdal
import numpy as np
from scipy.signal import convolve2d

def apply_kernel(infile, outfile, kernel):
	'''
	Apply the given kernel (a numpy 2-d array) to the raster, and save
	the result to the given output file.
	'''
	ds1 = gdal.Open(infile)
	bands = ds1.RasterCount
	cols = ds1.RasterXSize
	rows = ds1.RasterYSize

	ds2 = ds1.GetDriver().CreateCopy(outfile, ds1, strict = 0)
	for b in range(1, bands + 1):
		band = ds1.GetRasterBand(b)
		data = band.ReadAsArray()
		data = convolve2d(data, kernel, 'same')
		band = ds2.GetRasterBand(b)
		band.WriteArray(data)

if __name__ == '__main__':

	kernel = np.array([[1,1,1],[1,1,1],[1,1,1]]) / 9.

	infile = sys.argv[1]
	outfile = sys.argv[2]

	apply_kernel(infile, outfile, kernel)

