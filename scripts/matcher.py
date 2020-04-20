#!/usr/bin/env python3

import os
import sys
import shutil
import json

try:
	inputs_file = sys.argv[1]
	outputs_list = sys.argv[2]
	thread_count = int(sys.argv[3])
except: 
	print('Usage: matcher.py <input list> <output list> <thread count>')
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
	input_obj = json.loads(f.read())

inputs_list = input_obj['jobs']
g_usemask = input_obj.get('usemask', None)
mask = input_obj.get('mask', None)

# Outputs list -- same format as inputs but with the matche filesnames replacing the input
# filenames.

ol = open(outputs_list, 'w')

# Add processed files to the anchor set for the next iteration.
accum = True

# Skip the command invocation
skip = False

# Check that input files exist or will be created.
if True:
    print('Checking inputs...')
    available = set()
    for obj in inputs_list:
        if not obj.get('run', False):
            continue
        mergefile = obj['outfile']
        inputs = obj['infiles']
        for f,b in inputs:
            if not os.path.exists(f) and not f in available:
                print('file', f, 'does not exist and will not be created')
                sys.exit(1)
            available.add(f)
        available.add(mergefile)

# Iterate over the blocks.
for obj in inputs_list:
	
	if not obj.get('run', False):
		continue

	mergefile = obj['outfile']
	inputs = obj['infiles']
	usemask = obj.get('usemask', False) if g_usemask is None else g_usemask # Use if set on the object but override if global is set.
	
	trans_from, trans_to = obj.get('transform', [False, False])
	do_trans = trans_from and trans_to

	param_iters = len(obj['params'])
	param_iter = 0
	inter_files = []

	print('iterations', param_iters)

	for param_obj in obj['params']:

		param_iter += 1

		step_size = param_obj['step']
		radius = param_obj['radius']
		window_radius = param_obj['window']
		edges = param_obj.get('edges', False)
		method = param_obj['type']
		ncount = param_obj['count']
		exponent = param_obj['exponent']

		# Rename the first file as shifted (but with no shift); has only one band ,
		# which makes mosaicing easier.
		f1, b1 = inputs[0]
		if param_iter < param_iters:
			outfile = os.path.splitext(f1)[0] + '_shift_iter_{}.tif'.format(param_iter)
			inter_files.append(outfile)
		else:
			outfile = os.path.splitext(f1)[0] + '_shift.tif'

		if param_iter == 1:
			cmd = 'gdal_translate -b {b} {f1} {f2}'.format(b = b1, f1 = f1, f2 = outfile)
			print(cmd)
			os.system(cmd)
			ol.write('{f}:{b}\n'.format(f = f1, b = b1))
			inputs[0] = (outfile, 1)

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
		
			if param_iter < param_iters:
				outfile = os.path.splitext(f2)[0] + '_shift_iter_{}.tif'.format(param_iter)
				inter_files.append(outfile)
			else:
				outfile = os.path.splitext(f2)[0] + '_shift.tif'

			ol.write('{f}:{b}\n'.format(f = outfile, b = 1))
			inputs[i + 1] = (outfile, 1)

			# Configure params from the anchor list.
			# print(anchors)
			fparam = ' '.join(list(map(str, anchors)))

			# Compute the difference between the rasters.
			print('Computing difference...')

			cmd = 'rastermerge {x} {n} {r} {m} {t} {s} {d} {mask} {fp} {f2} {b2} {o}'.format(
				fp = fparam, 
				m = '-m {}'.format(method), 
				t = '-t {}'.format(thread_count), 
				s = '-s {}'.format(step_size), 
				r = '-r {}'.format(radius), 
				d = '-d {}'.format(window_radius),
				n = '-c {}'.format(ncount),
				x = '-x' if edges else '',
				b1 = b1, 
				f2 = f2, 
				b2 = b2, 
				o = outfile, 
				mask = '-k {k} {kb}'.format(k = mask, kb = 1) if usemask else ' '
			)
			print(cmd)

			if not skip and os.system(cmd) != 0:
				print('Canceled')
				sys.exit(1)
		
	cmd = 'gdalbuildvrt /tmp/tmp.vrt -tr 10 10 {i}'.format(i = ' '.join([x[0] for x in inputs[::-1]])) # TODO: Use averaging tool.
	print(cmd)

	if not skip and os.system(cmd) != 0:
		print('Canceled')
		sys.exit(1)

	cmd = 'gdal_translate /tmp/tmp.vrt {o}'.format(
		o = '/tmp/tmp.tif' if do_trans else mergefile
	)
	print(cmd)
	
	if not skip and os.system(cmd) != 0:
		print('Canceled')
		sys.exit(1)

	if do_trans:
		print('Transform', trans_from, trans_to)
		cmd = 'gdalwarp -s_srs {s} -t_srs {t} {f} {o}'.format(
			s = trans_from,
			t = trans_to,
			f = '/tmp/tmp.tif', 
			o = mergefile
		)
		print(cmd)
		if not skip and os.system(cmd) != 0:
			print('Transform Cancelled')
			sys.exit(1)

	#for filename in inter_files:
	# os.unlink(filename)

	ol.write('\n')

