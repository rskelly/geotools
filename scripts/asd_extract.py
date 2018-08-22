#!/usr/bin/env python

'''
Phase 1: convert ASD ViewSpec output file(s) to a csv file, columns to rows.
Phase 2: crunch the output from phase 1 to take the means of related spectral bands.

To execute phase 1, run the program as follows:

asd_extract.py phase1 dirname > output.csv

Before executing phase 2, two columns must be updated:
 - The set column is updated with a numerical ID that identifies groups of related rows. These rows
   will be averaged together band-by-band. 
 - The note column gives some information about the contents of the set, for example 'rad' for radiance.
   If the value of the note column is set to discard, that row will be ignored. A set ID need not be
   assigned for discarded columns.

To execute phase 2, run

asd_extract.py phase2 output.csv > final.csv
'''

import os
import sys
import re
import numpy as np

delim = '\t'
fext = '.txt'

def usage():
	print '''
Phase 1: convert ASD ViewSpec output file(s) to a csv file, columns to rows.
Phase 2: crunch the output from phase 1 to take the means of related spectral bands.

To execute phase 1, run the program as follows:

asd_extract.py phase1 dirname > output.csv

Before executing phase 2, two columns must be updated:
 - The set column is updated with a numerical ID that identifies groups of related rows. These rows
   will be averaged together band-by-band. 
 - The note column gives some information about the contents of the set, for example 'rad' for radiance.
   If the value of the note column is set to discard, that row will be ignored. A set ID need not be
   assigned for discarded columns.

To execute phase 2, run

asd_extract.py phase2 output.csv > final.csv
'''

def read_asd_file(filename):
	try:
		base, num, ext = os.path.basename(filename).split('.')
	except Exception as e:
		print(e)
		print filename
		raise e
	header = []
	spectrum = []
	with open(filename, 'rU') as f:
		for line in f:
			line = line.strip()
			if line.startswith('Wavelength'):
				break
			elif line:
				header.append(line)
		for line in f:
			spectrum.append(map(float, map(str.strip, line.strip().split(delim))))
	header = map(lambda x: x.replace('\x00', ''), header)
	return num, header, spectrum

def process_asd_dir(dirname):
	headrow = ['id', 'label', 'set', 'note', 'date', 'time', 'mean', 'stddev']

	first = True
	for f in [x for x in os.listdir(dirname) if x.endswith(fext)]:
		try:
			num, header, spectrum = read_asd_file(os.path.join(dirname, f))
		except:
			continue
		if first:
			headrow.extend([x[0] for x in spectrum])
			print(','.join(map(str, headrow)))
			first = False
		spec = np.array([x[1] for x in spectrum])
		mean = np.mean(spec)
		stddev = np.std(spec)
		dateg = re.search('Spectrum saved: (.+?) at (.+)', header[5])
		row = [num, header[2], '', '', dateg.group(1), dateg.group(2), mean, stddev]
		row.extend([x[1] for x in spectrum])
		print(','.join(map(str, row)))

def process_phase2(filename):
	'''
	Process the output from phase 1 which is the output of process_asd_dir.
	'''
	sets = {} # Contains the rows of data that will be computed as sets.
	firstwl = -1
	header = None

	# Load the data into the sets object. TODO: Memory intensive.
	with open(filename, 'rU') as f:
		for line in f:
			if not header:
				# Get the header. Will be used for locating elements in the data rows.
				header = map(str.strip, line.split(','))
			else:
				row = map(str.strip, line.split(','))
				note = row[header.index('note')]
				if note == 'discard':
					continue
				label = row[header.index('label')]
				setid = int(row[header.index('set')])
				curset = sets.get(setid, False)
				if not curset:
					curset = sets[setid] = {'label': '', 'note': '', 'setid': setid, 'rows': []}
				if label:
					curset['label'] = label
				if note:
					curset['note'] = note
				curset['rows'].append(row)

	# Find the first wl index
	for i in range(len(header)):
		try:
			int(header[i])
			firstwl = i
			break
		except: pass

	outhead = ['id', 'label', 'note', 'date', 'time', 'mean', 'stddev'] + header[firstwl:]
	print(','.join(outhead))

	# Get the column index of the first wavelength (the first integer field)
	for setid, obj in sets.iteritems():
		outrow = [setid, obj['label'], obj['note'], obj['rows'][0][header.index('date')], obj['rows'][0][header.index('time')]]

		# Create a multi-dimensional array of the values.
		vals = []
		for row in obj['rows']:
			vals.append(map(float, row[firstwl:]))
		vals = np.array(vals)

		# Calculate the mean of the set, columnwise, then take the mean and stddev of the 
		# resulting row.
		means = np.mean(vals, axis = 0)
		mean = np.mean(means)
		stddev = np.std(means)

		outrow.append(mean)
		outrow.append(stddev)
		outrow.extend(list(means))
		print(','.join(map(str, outrow)))


if __name__ == '__main__':

	phase = sys.argv[1]

	if phase == 'phase1':
		dirname = sys.argv[2]
		process_asd_dir(dirname)
	elif phase == 'phase2':
		filename = sys.argv[2]
		process_phase2(filename)
	else:
		usage()