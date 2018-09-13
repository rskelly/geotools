#!/usr/bin/env python3

import matplotlib.pyplot as plt
import os
import sys
import csv

# def load(csvfile):
# 	result = {}
# 	with open(csvfile, 'rU') as f:
# 		db = csv.reader(f)
# 		first = True
# 		for row in db:
# 			if first:
# 				for i in range(len(row)):
# 					result[i] = []
# 				first = False
# 			for i in range(0, len(row), 2):
# 				try:
# 					a = float(row[i])
# 					b = float(row[i + 1])
# 					result[i].append(a)
# 					result[i + 1].append(b)
# 				except Exception as e:
# 					#print(e) 
# 					#print(row)
# 					pass
# 	return result

def transpose(label, csvfile, headrow, datarow, col, delim, scale, xoffset):
	'''
	Transpose a head and data row from a csv file so that they become two columns.
	csvfile -- The csv file.
	headrow -- The row index of the header (0-based)
	datarow -- The row index of the data (0-based)
	col -- The column start index (0-based)
	'''
	headrow = int(headrow)
	datarow = int(datarow)
	col = int(col)
	scale = float(scale)

	if delim == 't':
		delim = '\t'
	elif delim == 's':
		delim = ' '

	a = []
	b = []
	with open(csvfile, 'rU') as f:
		db = csv.reader(f, delimiter = delim)
		r = 0
		for row in db:
			if r == headrow:
				r += 1
				a.extend(list(map(lambda x: float(x) + xoffset, row[col:])))
				break
			r += 1
		for row in db:
			if r == datarow:
				b.extend(list(map(lambda x: scale * float(x), row[col:])))
				break
			r += 1
	return (label, csvfile, a, b)


def load(label, csvfile, headrow, datarow, col, delim, scale, xoffset):
	'''
	csvfile -- The csv file.
	headrow -- The row index of the header (0-based)
	datarow -- The row index of the data (0-based)
	col -- The column start index (0-based)
	'''
	headrow = int(headrow)
	datarow = int(datarow)
	col = int(col)
	scale = float(scale)

	if delim == 't':
		delim = '\t'
	elif delim == 's':
		delim = ' '

	a = []
	b = []
	with open(csvfile, 'rU') as f:
		db = csv.reader(f, delimiter = delim)
		r = 0
		for row in db:
			if r < col:
				r += 1
				continue
			try:
				av = float(row[headrow]) + xoffset
				bv = float(row[datarow]) * scale
				a.append(av)
				b.append(bv)
			except Exception as e: 
				print(e)
				pass

	return (label, csvfile, a, b)

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

	for i in range(0, len(params), 9):
		t, label, csvfile, headrow, datarow, col, delim, scale, xoffset = params[i:i + 9]
		if t == 't':
			columns.append(transpose(label, csvfile, headrow, datarow, col, delim, scale, xoffset))
		elif t == 'tl':
			label, csvfile, a, b = transpose(label, csvfile, headrow, datarow, col, delim, scale, xoffset)
			a, b = smooth(a, b)
			columns.append((label, csvfile, a, b))
		elif t == 'tln':
			label, csvfile, a, b = transpose(label, csvfile, headrow, datarow, col, delim, scale, xoffset)
			a0, b0 = smooth(a, b)
			b1 = normalize(b0, b)
			columns.append((label, csvfile, a, b1))
		elif t == 'll':
			label, csvfile, a, b = load(label, csvfile, headrow, datarow, col, delim, scale, xoffset)
			a, b = smooth(a, b)
			columns.append((label, csvfile, a, b))
		elif t == 'tlne':
			label, csvfile, a, b = transpose(label, csvfile, headrow, datarow, col, delim, scale, xoffset)
			a0, b0 = smooth(a, b)
			b1 = normalize(b0, b)
			a2, b2, a3, b3 = find_extrema(25., a0, b1)
			columns.append((label + '_min', csvfile, a2, b2))
			columns.append((label + '_max', csvfile, a3, b3))
		else:
			columns.append(load(label, csvfile, headrow, datarow, col, delim, scale, xoffset))

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

		print(i, j, len(x))
		print(x[i], x[j])
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


if __name__ == '__main__':
	'''
	Call in the form:
	graph_transposed.py <t/n> <data file> <head row> <data row> <col> <delimiter> <mult> [<t/n> <data file> <head row> <data row> <col> <delimiter> <mult>]
	'''
	if len(sys.argv) == 1:
		args = []
		#args.extend(['t', 'irr', 'roof6_AbsoluteIrradiance_09-47-30-406_60s.txt', 13, 37, 2, 't', 1, 0.])
		#args.extend(['tl', 'irr_loess', 'roof6_AbsoluteIrradiance_09-47-30-406_60s.txt', 13, 37, 2, 't', 1, 0.])
		#args.extend(['tln', 'irr_norm_-8', 'roof6_AbsoluteIrradiance_09-47-30-406_60s.txt', 13, 37, 2, 't', 1, -8.])
		#args.extend(['tln', 'irr_norm_-8.5.', 'roof6_AbsoluteIrradiance_09-47-30-406_60s.txt', 13, 37, 2, 't', 1, -8.5])
		args.extend(['tln', 'irr_norm_-9', 'roof6_AbsoluteIrradiance_09-47-30-406_60s.txt', 13, 37, 2, 't', 1, -9.])
		#args.extend(['tlne', 'irr_extrema', 'roof6_AbsoluteIrradiance_09-47-30-406_60s.txt', 13, 37, 2, 't', 1, 0.])
		#args.extend(['t', 'convolved/roo6_convolved_nano_bands.csv', 0, 1, 2, ',', 1, 0.])
		args.extend(['tln', 'conv', 'convolved/roo6_convolved_nano_bands2.csv', 0, 1, 2, ',', 1, 0.])
		#args.extend(['t', 'nano', 'nano1.csv', 0, 1, 0, ',', 1, 0.])
		#args.extend(['tl', 'nano_loess', 'nano1.csv', 0, 1, 0, ',', 1, 0.])
		args.extend(['tln', 'nano_norm', 'nano1.csv', 0, 1, 0, ',', 3, 0.])
		#args.extend(['tlne', 'nano_extrema', 'nano1.csv', 0, 1, 0, ',', 10, 0.])
		#args.extend(['n', 'ASTMG173.csv', 0, 2, 2, ',', 15, 0.])
		graph(args)
	else:
		graph(sys.argv[1:])