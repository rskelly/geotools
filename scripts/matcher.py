#!/usr/bin/env python3

import os
import sys

inputs_list = [
[
	('/home/rob/Desktop/ec/leth_2016_n/L008-1-M01-S1-C1_r.tif', 2),
	('/home/rob/Desktop/ec/leth_2016_n/L007-1-M01-S1-C1_r.tif', 2),
	('/home/rob/Desktop/ec/leth_2016_n/L006-1-M01-S1-C1_r.tif', 2),
	('/home/rob/Desktop/ec/leth_2016_n/L005-1-M01-S1-C1_r.tif', 2),
	('/home/rob/Desktop/ec/leth_2016_n/L004-1-M01-S1-C1_r.tif', 2),
	('/home/rob/Desktop/ec/leth_2016_n/L003-1-M01-S1-C1_r.tif', 2),
	('/home/rob/Desktop/ec/leth_2016_n/L002-1-M01-S1-C1_r.tif', 2),
	('/home/rob/Desktop/ec/leth_2016_n/L001-1-M01-S1-C1_r.tif', 2),
],
[
	('/home/rob/Desktop/ec/leth_2016_n/L008-1-M01-S1-C1_r.tif', 2),
	('/home/rob/Desktop/ec/leth_2016_n/L009-1-M01-S1-C1_r.tif', 2),
	('/home/rob/Desktop/ec/leth_2016_n/L010-1-M01-S1-C1_r.tif', 2),
	('/home/rob/Desktop/ec/leth_2016_n/L011-1-M01-S1-C1_r.tif', 2),
	('/home/rob/Desktop/ec/leth_2016_n/L012-1-M01-S1-C1_r.tif', 2),
	('/home/rob/Desktop/ec/leth_2016_n/L013-1-M01-S1-C1_r.tif', 2),
	('/home/rob/Desktop/ec/leth_2016_n/L014-1-M01-S1-C1_r.tif', 2),
	('/home/rob/Desktop/ec/leth_2016_n/L015-1-M01-S1-C1_r.tif', 2),
	('/home/rob/Desktop/ec/leth_2016_n/L016-1-M01-S1-C1_r.tif', 2),
	('/home/rob/Desktop/ec/leth_2016_n/L017-1-M01-S1-C1_r.tif', 2),
]
]

size = 100

for inputs in inputs_list:
	for i in range(len(inputs) - 1):

		f1, b1 = inputs[i]
		f2, b2 = inputs[i + 1]
		
		tmp = 'tmp_{i}.tif'.format(i = i)

		print('Computing difference...')
		cmd = 'rastermerge -t 6 -s {s} {f1} {b1} {f2} {b2} {o}'.format(s = size, f1 = f1, b1 = b1, f2 = f2, b2 = b2, o = tmp)
		if os.system(cmd) != 0:
			print('Canceled')
			sys.exit(1)

		outfile = os.path.splitext(f2)[0] + '_shift.tif'
		print('Calculating shift...', outfile)
		cmd = 'gdal_calc.py -A {f2} --A_band={b2} -B {t} --calc="(A + B)" --overwrite --outfile="{f}"'.format(f2 = f2, b2 = b2, f = outfile, t = tmp)
		if os.system(cmd) != 0:
			print('Canceled')
			sys.exit(1)

		inputs[i + 1] = (outfile, 1)
