#!/usr/bin/env python3

import os
import sys
import gdal
import matplotlib.pyplot as plt
import numpy as np

def load_csv_col(filename, col, where = []):
	'''
	Loads a column of data from the spreadsheet, where 
	col is the zero-based column index. Where is an optional
	list of tuples containing a column index and a value that
	must be equal to that column's value for a given row
	to be accepted. For example, a column might have a category
	name or ID that can be filtered on.

	Returns a list of floats.
	'''
	data = []
	with open(filename, 'rU') as f:
		for line in f:
			line = line.strip().split(',')
			keep = True
			for c, v in where:
				if line[c] != v:
					keep = False
					break
			if keep:
				data.append(line[col])
	return list(map(float, data))


def load_csv_head(filename):
	'''
	Load the header row from a csv file.
	'''
	with open(filename, 'rU') as f:
		return f.readline().strip().split(',')

def sample_raster(rfilename, pfilename, ofilename, labels = ('light', 'dark'), label_idx = 3):
	'''
	Extract the pixel values from the given raster, from the
	given spreadsheet. The spreadsheet will have the format, 
		x, y, id, target
	where target is the name of the target used (e.g. light or dark).
	
	Coordinates are pixel values and are truncated before use.

	Writes the sampled values to the ofilename. Format is the same as the point
	file with bands appended.
	'''
	if os.path.exists(ofilename):
		raise Exception('Error! The output filename, ' + ofilename + ', already exists.')
	
	ds = gdal.Open(rfilename)
	cols = ds.RasterXSize
	rows = ds.RasterYSize
	bands = ds.RasterCount

	samples = []
	first = True
	with open(pfilename, 'rU') as f:
		for row in f:
			if first:
				first = False
				continue
			samples.append(row.strip().split(','))

	# Read the sample band values into an array.
	for i in range(len(samples)):
		x, y, id, target = samples[i]
		for b in range(bands):
			band = ds.GetRasterBand(b + 1)
			v = band.ReadAsArray(int(float(x)), int(rows - float(y)), 1, 1)
			samples[i].append(v[0,0])
	
	# For light and dark, compute the means of all samples for each band.
	for t in labels:
		sums = [0.] * len(samples[0])
		counts = [0] * len(samples[0])
		for i in range(len(samples)):
			if samples[i][label_idx] == t:
				for j in range(4, len(samples[i])):
					sums[j] += samples[i][j]
					counts[j] += 1
		sums[label_idx] = 'mean_' + t
		for i in range(4, len(samples[0])):
			sums[i] /= counts[i]
		samples.append(sums)

	with open(ofilename, 'w') as f:
		f.write('x,y,id,label,' + ','.join(['b_' + str(x) for x in range(bands)]) + '\n')
		f.write('\n'.join([','.join(list(map(str, r))) for r in samples]))


def mean(lst):
	'''
	Compute the mean of a list.
	'''
	return float(sum(lst)) / len(lst)


def load_asd_nano(asd_filename, nano_filename, 
	asd_offset, nano_offset, asd_label_idx, nano_label_idx,
	asd_labels, nano_labels, num_bands):
	'''
	Load data from the ASD and NANO and return the data.
	'''
	# Load the headers.
	nano_head = load_csv_head(nano_filename)
	asd_head = load_csv_head(asd_filename)

	# Load selected rows and put the data in the arrays.
	nano_y_light = []
	nano_y_dark = []
	asd_y_light = []
	asd_y_dark = []
	for b in range(num_bands):
		nano_y_light.append(mean(load_csv_col(nano_filename, b + nano_offset, [(nano_label_idx, nano_labels[0],)])))
		nano_y_dark.append(mean(load_csv_col(nano_filename, b + nano_offset, [(nano_label_idx, nano_labels[1])])))

		asd_y_light.append(mean(load_csv_col(asd_filename, b + asd_offset, [(asd_label_idx, asd_labels[0],)])))
		asd_y_dark.append(mean(load_csv_col(asd_filename, b + asd_offset, [(asd_label_idx, asd_labels[1],)])))

	count = len(nano_head[nano_offset:])
	return (asd_head, nano_head, asd_y_light, asd_y_dark, nano_y_light, nano_y_dark, count)


def plot(nano_filename, asd_filename, 
	asd_offset, nano_offset, asd_label_idx, nano_label_idx,
	asd_labels, nano_labels, num_bands):
	'''
	Plot the spectra from the nano and the asd.
	'''

	asd_head, nano_head, asd_y_light, asd_y_dark, nano_y_light, nano_y_dark, count = load_asd_nano(asd_filename, nano_filename,
		asd_offset, nano_offset, asd_label_idx, nano_label_idx, 
		asd_labels, nano_labels, num_bands)

	ax0 = plt.subplot()
	#ax0.set_ylim([-10., 40.])

	head = list(map(float, asd_head[asd_offset:asd_offset + count]))
	ax0.plot(head, nano_y_light[:count], label = 'Light')
	ax0.plot(head, nano_y_dark[:count], label = 'Dark')

	ax1 = ax0.twinx()

	ax1.plot(head, asd_y_light[:count], label = 'Light')
	ax1.plot(head, asd_y_dark[:count], label = 'Dark')

	plt.legend()
	plt.show()


def regress(nano_filename, asd_filename, out_filename, asd_offset, nano_offset, asd_label_idx, nano_label_idx,
	asd_labels, nano_labels, num_bands):
	'''
	Plot the light-dark connections between pairs for each band, where
	each pair consists of a member from ASD and the corresponding member
	from Nano (by band).

	The Nano is on the y-axis because it's the dependent outcome: its
	value is a function of the ASD reflectance (which is assumed to have no
	atmospheric component), path radiance and some coefficient (the slope.)

	Computes and writes the slope and y-intercept for each pair to the output database.
	
	nano_filename - The radiance image produced by the spectrometer.
	asd_filename - The convolved ASD spectra.
	out_filename - The output filename, which will contain the coefficients.
	'''
	
	if os.path.exists(out_filename):
		raise Exception('Error! The output file, ' + out_filename + ', already exists.')

	asd_head, nano_head, asd_y_light, asd_y_dark, nano_y_light, nano_y_dark, count = load_asd_nano(asd_filename, nano_filename,
		asd_offset, nano_offset, asd_label_idx, nano_label_idx, 
		asd_labels, nano_labels, num_bands)

	xl = []
	yl = []
	xd = []
	yd = []

	ax0 = plt.subplot()
	#ax0.set_xlim([-10., 60.])

	head = list(map(float, asd_head[asd_offset:]))

	coeffs = []
	for b in range(161, 169): #range(count):
		xl = asd_y_light[b]
		yl= nano_y_light[b]
		xd = asd_y_dark[b]
		yd = nano_y_dark[b]
		
		if head[b] <= 900.:
			ax0.plot([xd, xl], [yd, yl], label = '{w}nm'.format(w = head[b]))
			ax0.plot([xd], [yd], 'ob', label = 'Dark')
			ax0.plot([xl], [yl], 'or', label = 'Light')

		try:
			m = (yl - yd) / (xl - xd)
			yi = yd - m * xd

			coeffs.append((b, head[b], xd, yd, xl, yl, m, yi))
		except: pass

	with open(out_filename, 'w') as f:
		f.write('band_idx,band_wl,asd_dark,nano_dark,asd_light,nano_light,slope,y_int\n')
		f.write('\n'.join([','.join(list(map(str, x))) for x in coeffs]))

	plt.legend()
	plt.show()

def transform_raster(irasterfile, coeffile, orasterfile):
	'''
	Apply the coefficients determined in the regress step to the raster.
	irasterfile - The apparent reflectance raster.
	coeffile - The file containing coefficients, produced by the regress step.
	orasterfile - The transformed raster output file.
	'''
	
	if os.path.exists(orasterfile):
		raise Exception('Error! The output file, ' + orasterfile + ', already exists.')

	with open(coeffile, 'rU') as f:
		f.readline()
		coeffs = [x.split(',') for x in f.read().split('\n')]

	dsin = gdal.Open(irasterfile, gdal.GA_ReadOnly)
	cols = dsin.RasterXSize
	rows = dsin.RasterYSize
	bands = dsin.RasterCount

	cbands = len(coeffs)

	drv = gdal.GetDriverByName('ENVI')
	dsout = drv.Create(orasterfile, cols, rows, bands, gdal.GDT_Float32)
	dsout.SetMetadata(dsin.GetMetadata())
	
	for b in range(cbands):
		# Input and output bands.
		iband = dsin.GetRasterBand(b + 1)
		oband = dsout.GetRasterBand(b + 1)
		idata = iband.ReadAsArray()
		# Slope and intercept from coefficients for this band.
		m = float(coeffs[b][6])
		y = float(coeffs[b][7])
		odata = (idata - y) / m
		oband.WriteArray(odata)


def usage():
	print('''Run the program multiple times with the following arguments, in order.
			* The raster file is apparent reflectance, generated from the spectometer,
			  and divided by the irradiance.
		1) refl_regress.py sample <raster file> <sample point file> <sample value file>
			Samples the raster for every band at each point. The format for the
			sample file is x, y, id, target. Writes to the sample value file.
		2) refl_regress.py plot <sample value file> <asd file>
			Plots the correspondence between the sample value file and the ASD reflectance.
		3) refl_regress.py regress <sample value file> <asd file> <coefficients file>
			Performs the linear regression analysis on the samples and ASD reflectance
			values.
		4) refl_regress.py transform <raster file> <transformed raster file> <coefficients file>
			Transforms the apparent reflectance file, removing path radiance and scaling by the slope.
			Outputs a new raster.
			
		Help for this tool and the processing steps can be found at 
			https://github.com/rskelly/contrem/wiki/refl_regress
			https://github.com/rskelly/contrem/wiki/Reflectance-Coefficients
	''')

if __name__ == '__main__':

	try:
		cmd = sys.argv[1]
	except:
		usage()
		sys.exit(1)

	asd_offset = 3
	nano_offset = 4
	asd_label_idx = 0
	nano_label_idx = 3
	asd_labels = ('mean_light', 'mean_dark')
	nano_labels = ('mean_light', 'mean_dark')
	num_bands = 272

	try:
		if cmd == 'sample':
			sample_raster(*sys.argv[2:5])
		elif cmd == 'plot':
			plot(sys.argv[2], sys.argv[3], asd_offset, nano_offset, asd_label_idx, nano_label_idx, 
			asd_labels, nano_labels, num_bands)
		elif cmd == 'regress':
			regress(sys.argv[2], sys.argv[3], sys.argv[4], asd_offset, nano_offset, asd_label_idx, nano_label_idx, asd_labels, nano_labels, num_bands)
		elif cmd == 'transform':
			transform_raster(*sys.argv[2:5])
		else:
			usage()
	except Exception as e:
		print(e)
		usage()
