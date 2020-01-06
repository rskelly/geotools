#!/usr/bin/env python3

import os
import sys
import shutil

mask = '/media/rob/robdata/work/ec/DEM/mask_sd0_1_ct20_f35.tif'

try:
	method = sys.argv[1]
	resamp = sys.argv[2]
	inputs_file = sys.argv[3]
	outputs_list = sys.argv[4]
	window_size = int(sys.argv[5])
	thread_count = int(sys.argv[6])
except: 
	print('Usage: matcher.py <type> <resample (1 for no effect)> <input list> <output list> <window size> <thread count>')
	print(' types: idw, gauss, cosine. Gauss, cosine decay to zero, idw doesn\'t.')
	sys.exit(1)

# Load the inputs list. A text file with entries in the form:
#
# filename:band
# 
# One on each line. The first item is the anchor, each item is matched to it, and the new
# file replaces the old one in the list (except for the first).
# 
# A blank line signifies the start of a new list of inputs.

inputs_list = []

with open(inputs_file, 'r') as f:
	inputs = []
	line = f.readline()
	while line:
		line = line.strip()
		if line == '':
			if len(inputs):
				inputs_list.append(inputs)
				inputs = []
		elif line.startswith('#'):
			pass # Comment
		else:
			inputs.append(line.split(':'))
		line = f.readline()

	if len(inputs):
		inputs_list.append(inputs)


# Outputs list -- same format as inputs but with the matche filesnames replacing the input
# filenames.

ol = open(outputs_list, 'w')

# Add processed files to the anchor set for the next iteration.
accum = True

# Iterate over the blocks.
for inputs in inputs_list:

	# Iterate over the pairs of files.
	for i in range(len(inputs) - 1):

		anchors = []

		if accum:
			for j in range(i + 1):
				f1, b1 = inputs[j]
				anchors.extend([f1, b1])
		else:
			f1, b1 = inputs[i]
			anchors.extend([f1, b1])

		f2, b2 = inputs[i + 1]

		outfile = os.path.splitext(f2)[0] + '_shift.tif'

		# Write out the first one unchanged, then every subsequent one with the shift.
		if i == 0:
			ol.write('{f}:{b}\n'.format(f = f1, b = b1))
		ol.write('{f}:{b}\n'.format(f = outfile, b = 1))

		# Configure params from the anchor list.
		print(anchors)
		fparam = ' '.join(list(map(str, anchors)))

		# Compute the difference between the rasters.
		print('Computing difference...')
		tmp = outfile #'tmp.tif'
		cmd = 'rastermerge -d -r {r} -m {m} -t {t} -s {s} -k {k} {kb} {f} {f2} {b2} {o}'.format(k = mask, kb = 1, r = resamp, m = method, t = thread_count, s = window_size, f = fparam, b1 = b1, f2 = f2, b2 = b2, o = tmp)
		print(cmd)

		if os.system(cmd) != 0:
			print('Canceled')
			sys.exit(1)
		
		# Add the difference back into the source file to create the new matched file.
		#print('Calculating shift...', outfile)
		#cmd = 'gdal_calc.py -A {f2} --A_band={b2} -B {t} --calc="(A + B)" --overwrite --outfile="{f}"'.format(f2 = f2, b2 = b2, f = outfile, t = tmp)
		#if os.system(cmd) != 0:
		#	print('Canceled')
		#	sys.exit(1)

		#shutil.move(tmp, outfile)

		inputs[i + 1] = (outfile, 1)

	ol.write('\n')