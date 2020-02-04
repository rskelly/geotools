#!/usr/bin/env python3

import os
import sys
import shutil
import json

mask = '/media/rob/robdata/work/ec/DEM/mask_sd0_1_ct20_f35.tif'

try:
	method = sys.argv[1]
	inputs_file = sys.argv[2]
	outputs_list = sys.argv[3]
	window_size = int(sys.argv[4])
	neighbours = int(sys.argv[5])
	thread_count = int(sys.argv[6])
except: 
	print('Usage: matcher.py <method> <input list> <output list> <window size> <neighbours> <thread count>')
	print(' types: dw, idw, gauss, cosine. Gauss, cosine decay to zero, idw doesn\'t.')
	sys.exit(1)

# Load the inputs list. A text file with entries in the form:
#
# filename:band
# 
# One on each line. The first item is the anchor, each item is matched to it, and the new
# file replaces the old one in the list (except for the first).
# 
# A blank line signifies the start of a new list of inputs.


with open(inputs_file, 'r') as f:
	inputs_list = json.loads(f.read())


# Outputs list -- same format as inputs but with the matche filesnames replacing the input
# filenames.

ol = open(outputs_list, 'w')

# Add processed files to the anchor set for the next iteration.
accum = True
skip = False

# Iterate over the blocks.
for obj in inputs_list:
	
	if not obj.get('run', False):
		continue

	mergefile = obj['outfile']
	inputs = obj['infiles']

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

		# Rename the first file as shifted (but with no shift); has only one band ,
		# which makes mosaicing easier.
		if i == 0:
			outfile = os.path.splitext(f1)[0] + '_shift.tif'
			cmd = 'gdal_translate -b {b} {f1} {f2}'.format(b = b1, f1 = f1, f2 = outfile)
			print(cmd)
			os.system(cmd)
			ol.write('{f}:{b}\n'.format(f = f1, b = b1))
			inputs[i] = (outfile, 1)
	
		outfile = os.path.splitext(f2)[0] + '_shift.tif'
		ol.write('{f}:{b}\n'.format(f = outfile, b = 1))
		inputs[i + 1] = (outfile, 1)

		# Configure params from the anchor list.
		# print(anchors)
		fparam = ' '.join(list(map(str, anchors)))

		# Compute the difference between the rasters.
		print('Computing difference...')
		msk = ' -k {k} {kb} '.format(k = mask, kb = 1)
		cmd = 'rastermerge -m {m} -t {t} -s {s} -c {c}{mask}{f} {f2} {b2} {o}'.format(c = neighbours, m = method, t = thread_count, s = window_size, f = fparam, b1 = b1, f2 = f2, b2 = b2, o = outfile, mask = msk)
		print(cmd)

		if not skip and os.system(cmd) != 0:
			print('Canceled')
			sys.exit(1)
		
	cmd = 'gdalbuildvrt /tmp/tmp.vrt {i}'.format(i = ' '.join([x[0] for x in inputs]))
	print(cmd)

	if not skip and os.system(cmd) != 0:
		print('Canceled')
		sys.exit(1)

	cmd = 'gdal_translate /tmp/tmp.vrt {o}'.format(o = mergefile)
	print(cmd)
	
	if not skip and os.system(cmd) != 0:
		print('Canceled')
		sys.exit(1)

	ol.write('\n')

