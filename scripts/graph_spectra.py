#!/usr/bin/env python3

'''
This script is primarly for graphing spectra from delimited text files. It can
load column-based files, and load and transpose row-based files.

Each method takes a Params object with properties that may or may not be used
by that method. In some cases the properties may be interpreted differently by
the method, so pay attention to the doc comments.
'''

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

	Return the label, the filename, the header column and data column.
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
	Load a head and data row from a csv file so that they become two columns.
	
	csvfile -- The csv file.
	headrow -- The row index of the header (0-based)
	datarow -- The row index of the data (0-based)
	col -- The column start index (0-based)

	Return the label, the filename, the header column and data column.
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

	x - The column containing x values.
	y - The column containing y values.
	'''
	from scipy.interpolate import interp1d
	import statsmodels.api as sm

	lowess = sm.nonparametric.lowess(y, x, frac = .3)
	lx = list(zip(*lowess))[0]
	ly = list(zip(*lowess))[1]

	return lx, ly

def normalize(y, ybase):
	'''
	Just subtracts ybase from y, element-wise. The ybase
	argument must correspond element-wise to the y argument.

	y - The y values to be normalized.
	ybase - The y values to normalize against.
	'''
	return [y0 - ybase0 for y0, ybase0 in zip(y, ybase)]

def graph(params):
	'''
	Produce a graph using the spectra loaded using the params given in the list.

	The type property of each param object determines how the data are loaded
	and processed. The following types are available.

	[none; default] - Load the dataset.
	t - Load and transpose the dataset.
	tl - Load, transpose and smooth the dataset using loess.
	tln - Load, transpose, smooth and normalize.
	ll - Load and smooth the dataset.
	tlne - Load, transpose, smooth, normalize and graph the extrema.
	'''
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
	Find the extrema within a given window along the
	x and y data. Not really useful now.
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
	''' 
	A parameter container for the methods in this script.
	'''

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
	graph_transposed.py <config file>

	The configuration file is a text file. Each line contains a space-delimited list of arguments. This
	program is terribly naive, so none of the arguments can have spaces. The arguments are, 

	<type> <data file> <head row> <data row> <col> <delimiter> <mult> <shift>

	Any row that begins with a space or # is ignored.
	'''
	params = []
	with open(sys.argv[1], 'r') as f:
		for line in f:
			line = line.strip()
			if line == '\n' or line == '' or line.startswith(' ') or line.startswith('#'):
				continue
			args = line.split(' ')
			print(args)
			params.append(Params(*args))
	graph(params)