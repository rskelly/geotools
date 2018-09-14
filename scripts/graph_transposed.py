#!/usr/bin/env python3

import matplotlib.pyplot as plt
import os
import sys
import csv

def transpose(params):
	'''
	Transpose a head and data row from a csv file so that they become two columns.
	csvfile -- The csv file.
	headrow -- The row index of the header (0-based)
	datarow -- The row index of the data (0-based)
	col -- The column start index (0-based)
	'''
	headrow = int(params.head_offset)
	datarow = int(params.data_offset)
	col = int(params.col_offset)
	scale = float(params.scale)
	shift = float(params.shift)

	delim = params.delim
	if delim == 't':
		delim = '\t'
	elif delim == 's':
		delim = ' '

	a = []
	b = []
	with open(params.file, 'rU') as f:
		db = csv.reader(f, delimiter = delim)
		r = 0
		for row in db:
			if r == headrow:
				r += 1
				a.extend(list(map(lambda x: float(x) + shift, row[col:])))
				break
			r += 1
		for row in db:
			if r == datarow:
				b.extend(list(map(lambda x: scale * float(x), row[col:])))
				break
			r += 1
	return (params.label, params.file, a, b)


def load(params):
	'''
	csvfile -- The csv file.
	headrow -- The row index of the header (0-based)
	datarow -- The row index of the data (0-based)
	col -- The column start index (0-based)
	'''
	headrow = int(params.head_offset)
	datarow = int(params.data_offset)
	col = int(params.col_offset)
	scale = float(params.scale)
	shift = float(params.shift)

	delim = params.delim
	if delim == 't':
		delim = '\t'
	elif delim == 's':
		delim = ' '

	a = []
	b = []
	with open(params.file, 'rU') as f:
		db = csv.reader(f, delimiter = delim)
		r = 0
		for row in db:
			if r < col:
				r += 1
				continue
			try:
				av = float(row[headrow]) + shift
				bv = float(row[datarow]) * scale
				a.append(av)
				b.append(bv)
			except Exception as e: 
				print(e)
				pass

	return (params.label, params.file, a, b)

def smooth(x, y):
	'''
	Use local linear regression to produce a smoothed plot which can be used for normalization.
	'''
	from scipy.interpolate import interp1d
	import statsmodels.api as sm

	lowess = sm.nonparametric.lowess(y, x, frac = .3)
	lx = list(zip(*lowess))[0]
	ly = list(zip(*lowess))[1]

	return lx, ly

def normalize(y, ybase):
	'''
	Just subtracts ybase from y, element-wise.
	'''
	return [y0 - ybase0 for y0, ybase0 in zip(y, ybase)]

def graph(params):

	columns = []

	for param in params:
		if param.type == 't':
			columns.append(transpose(param))
		elif param.type == 'tl':
			label, csvfile, a, b = transpose(param)
			a, b = smooth(a, b)
			columns.append((label, csvfile, a, b))
		elif param.type == 'tln':
			label, csvfile, a, b = transpose(param)
			a0, b0 = smooth(a, b)
			b1 = normalize(b0, b)
			columns.append((label, csvfile, a, b1))
		elif param.type == 'll':
			label, csvfile, a, b = load(param)
			a, b = smooth(a, b)
			columns.append((label, csvfile, a, b))
		elif param.type == 'tlne':
			label, csvfile, a, b = transpose(param)
			a0, b0 = smooth(a, b)
			b1 = normalize(b0, b)
			a2, b2, a3, b3 = find_extrema(25., a0, b1)
			columns.append((label + '_min', csvfile, a2, b2))
			columns.append((label + '_max', csvfile, a3, b3))
		else:
			columns.append(load(param))

	minx = 99999999999.
	maxx = -99999999999.
	for label, file, a, b in columns:
		minx = min(min(a), minx)
		maxx = max(max(a), maxx)

	# Plot the x axis.
	plt.plot([minx, maxx], [0, 0])

	for label, file, a, b in columns:
		print(label, len(a), len(b))
		plt.plot(a, b, label = label)

	plt.legend()
	plt.show()


def find_extrema(window_size, x, y):
	'''
	'''
	minx0 = []
	miny0 = []
	maxx0 = []
	maxy0 = []

	i = 0
	while i < len(x):

		j = i + 1
		if j >= len(x):
			break
		while (j < len(x) - 1) and (x[j] <= x[i] + window_size):
			j += 1

		maxy = -99999999.
		maxx = 0
		miny = 99999999.
		minx = 0
		for k in range(i, j):
			if y[k] < minx:
				miny = y[k]
				minx = x[k]
		for k in range(i, j):
			if y[k] > maxx:
				maxy = y[k]
				maxx = x[k]

		if miny < 99999999.:
			minx0.append(minx)
			miny0.append(miny)
		if maxy > -99999999.:
			maxx0.append(maxx)
			maxy0.append(maxy)

		i = j

	return (minx0, miny0, maxx0, maxy0)

class Params:

	def __init__(self, type, label, file, head_offset, data_offset, col_offset, delim = ',', scale = 1., shift = 0.):
		self.type = type
		self.label = label
		self.file = file
		self.head_offset = head_offset
		self.data_offset = data_offset
		self.col_offset = col_offset
		self.delim = delim
		self.scale = scale
		self.shift = shift


if __name__ == '__main__':
	'''
	Call in the form:
	graph_transposed.py <t/n> <data file> <head row> <data row> <col> <delimiter> <mult> <shift> [<t/n> <data file> <head row> <data row> <col> <delimiter> <mult>]
	'''
	if len(sys.argv) == 1:
		args = []
		#args.append(['t', 'irr_norm', 'roof6_AbsoluteIrradiance_09-47-30-406_60s.txt', 13, 37, 2, 't', 1, -9.])
		#args.append(['t', 'irr_raw', 'roof6_AbsoluteIrradiance_09-47-30-406.txt', 14, 15, 2, 't', 1, 0.])
		#args.append(['t', 'shift_test', 'shift_AbsoluteIrradiance_09-44-34-037.txt', 14, 15, 2, 't', 100., 0.])
		#args.append(['t', 'shift_test', 'shift2_FLMS128791_10-14-01-727.txt', 14, 15, 2, 't', .01, 0.])

		#args.append(['t', 'conv', 'convolved/roo6_convolved_nano_bands2.csv', 0, 1, 2, ',', 1, 0.])

		# Sample NANO spectrum.
		#args.append(['t', 'nano', 'roof_irradiance_convolved/nano1.csv', 0, 1, 0, ',', 3, 0.])

		# ASTM Reference spectrum scaled by 15.
		args.append(Params('n', 'ASTMG173', 'roof_irradiance_convolved/ASTMG173.csv', 0, 2, 2, ',', 15, 0.))

		# Bartier 1, with shift applied.
		#args.append(Params('t', 'bartier_1_1', '/home/rob/Documents/field/2018_sept/flame/bartier1_AbsoluteIrradiance_10-09-26-788.txt', 14, 15, 2, 't', 1., -7.97))

		# Bartier 1, shifted and convolved.
		#args.append(Params('t', 'bartier_1_c_1', '/home/rob/Documents/field/2018_sept/flame/bartier1_AbsoluteIrradiance_10-09-26-788_conv_nano.csv', 0, 1, 2, ',', 1., 0.))

		# Bartier 2, with shift applied.
		#args.append(Params('t', 'bartier_2_1', '/home/rob/Documents/field/2018_sept/flame/bartier2_AbsoluteIrradiance_11-20-47-135.txt', 14, 15, 2, 't', 1., -7.97))

		# Bartier 2, shifted and convolved.
		#args.append(Params('t', 'bartier_2_c_1', '/home/rob/Documents/field/2018_sept/flame/bartier2_AbsoluteIrradiance_11-20-47-135_conv_nano.txt', 0, 1, 2, ',', 1., 0.))

		# Bartier 3, with shift applied.
		#args.append(Params('t', 'bartier_3_1', '/home/rob/Documents/field/2018_sept/flame/bartier3_AbsoluteIrradiance_13-05-18-419.txt', 14, 15, 2, 't', 1., -7.97))

		# Bartier 3, shifted and convolved.
		#args.append(Params('t', 'bartier_3_c_1', '/home/rob/Documents/field/2018_sept/flame/bartier3_AbsoluteIrradiance_13-05-18-419_conv_nano.txt', 0, 1, 2, ',', 1., 0.))

		# SF 1, with shift applied.
		args.append(Params('t', 'sf_1', '/home/rob/Documents/field/2018_sept/flame/sf10_1_AbsoluteIrradiance_10-40-58-450.txt', 14, 1015, 2, 't', 1., -7.97))

		# SF 1, shifted and convolved.
		args.append(Params('t', 'sf_c_1', '/home/rob/Documents/field/2018_sept/flame/sf10_1_AbsoluteIrradiance_10-40-58-450_conv_nano.txt', 0, 1001, 2, ',', 1., 0.))

		# SF, with shift applied.
		args.append(Params('t', 'sf', '/home/rob/Documents/field/2018_sept/flame/sf10_AbsoluteIrradiance_13-02-17-800.txt', 14, 1015, 2, 't', 1., -7.97))

		# SF, shifted and convolved.
		args.append(Params('t', 'sf_c', '/home/rob/Documents/field/2018_sept/flame/sf10_AbsoluteIrradiance_13-02-17-800_conv_nano.txt', 0, 1000, 2, ',', 1., 0.))

		graph(args)
	else:
		graph(sys.argv[1:])