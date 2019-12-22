#!/usr/bin/env python3

import os
import sys

try:
	inputs_file = sys.argv[1]
	outputs_list = sys.argv[2]
	window_size = int(sys.argv[3])
	thread_count = int(sys.argv[4])
except: 
	print('Usage: matcher.py <input list> <output list> <window size> <thread count>')
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
		else:
			inputs.append(line.split(':'))
		line = f.readline()

	if len(inputs):
		inputs.append(inputs)


# Outputs list -- same format as inputs but with the matche filesnames replacing the input
# filenames.

ol = open(outputs_list, 'w')

# Iterate over the blocks.
for inputs in inputs_list:

	# Iterate over the pairs of files.
	for i in range(len(inputs) - 1):

		f1, b1 = inputs[i]
		f2, b2 = inputs[i + 1]

		outfile = os.path.splitext(f2)[0] + '_shift.tif'

		# Write out the first one unchanged, then every subsequent one with the shift.
		if i == 0:
			ol.write('{f}:{b}\n'.format(f = f1, b = b1))
		ol.write('{f}:{b}\n'.format(f = outfile, b = 1))

		# Compute the difference between the rasters.
		print('Computing difference...')
		tmp = 'tmp.tif'
		cmd = 'rastermerge -t {t} -s {s} {f1} {b1} {f2} {b2} {o}'.format(t = thread_count, s = window_size, f1 = f1, b1 = b1, f2 = f2, b2 = b2, o = tmp)
		if os.system(cmd) != 0:
			print('Canceled')
			sys.exit(1)

		# Add the difference back into the source file to create the new matched file.
		print('Calculating shift...', outfile)
		cmd = 'gdal_calc.py -A {f2} --A_band={b2} -B {t} --calc="(A + B)" --overwrite --outfile="{f}"'.format(f2 = f2, b2 = b2, f = outfile, t = tmp)
		if os.system(cmd) != 0:
			print('Canceled')
			sys.exit(1)

		inputs[i + 1] = (outfile, 1)

	ol.write('\n')