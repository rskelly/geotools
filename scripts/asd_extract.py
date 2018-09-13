#!/usr/bin/env python3

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
import shutil
import json
import numpy as np

delim = ',' #\t'
fext = '.txt'

def usage():
	print('''
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
''')


def read_asd_file(filename):
	try:
		base, num, ext = os.path.basename(filename).split('.')
	except Exception as e:
		print(e)
		print(filename)
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
			spectrum.append(list(map(float, list(map(str.strip, line.strip().split(delim))))))
	header = list(map(lambda x: x.replace('\x00', ''), header))
	return num, header, spectrum

def process_phase1(dirname):
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
	Process the output from phase 1 which is the output of process_phase1.
	'''
	sets = {} # Contains the rows of data that will be computed as sets.
	firstwl = -1
	header = None

	# Load the data into the sets object. TODO: Memory intensive.
	with open(filename, 'rU') as f:
		for line in f:
			if not header:
				# Get the header. Will be used for locating elements in the data rows.
				header = list(map(str.strip, line.split(',')))
			else:
				row = list(map(str.strip, line.split(',')))
				note = row[header.index('note')]
				if note == 'discard':
					continue
				label = row[header.index('label')]
				try:
					setid = int(row[header.index('set')])
				except:
					setid = 0
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
	for setid, obj in sets.items():
		outrow = [setid, obj['label'], obj['note'], obj['rows'][0][header.index('date')], obj['rows'][0][header.index('time')]]

		# Create a multi-dimensional array of the values.
		vals = []
		for row in obj['rows']:
			vals.append(list(map(float, row[firstwl:])))
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

def process_web(filename, outdir):
	'''
	Create a navigable page from phase2 data.
	'''
	sets = {} # Contains the rows of data that will be computed as sets.
	firstwl = -1
	header = None

	# Load the data into the sets object. TODO: Memory intensive.
	with open(filename, 'rU') as f:
		for line in f:
			if not header:
				# Get the header. Will be used for locating elements in the data rows.
				header = list(map(str.strip, line.split(',')))
			else:
				row = list(map(str.strip, line.split(',')))
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

	try:
		for f in os.listdir(os.path.join(outdir, 'plots')):
			shutil.rmtree(os.path.join(outdir, 'plots', f))
	except: pass

	output = []

	# Get the column index of the first wavelength (the first integer field)
	for setid, obj in sets.iteritems():
		outrow = {
			'id' : setid, 
			'label' : obj['label'], 
			'note' : obj['note'], 
			'date' : obj['rows'][0][header.index('date')], 
			'time' : obj['rows'][0][header.index('time')],
			'count' : len(obj['rows'])
		}

		vals = [] # Multi-dimensional array of the values.
		for row in obj['rows']:
			vals.append(list(map(float, row[firstwl:])))
		vals = np.array(vals)

		# Calculate the mean of the set, columnwise, then take the mean and stddev of the resulting row.
		means = np.mean(vals, axis = 0)
		mean = np.mean(means)
		stddev = np.std(means)

		# Compute the trend between rows using the mean along axis 1
		trend = np.mean(vals, axis = 1)

		# Append stats.
		outrow['mean'] = mean
		outrow['stddev'] = stddev
		#outrow['means'] = list(means)

		# Generate an ID that will be used for the plots.
		rowid = gen_id()
		outrow['rowid'] = rowid

		plot_dir = os.path.join(outdir, 'plots', rowid)
		try:
			os.makedirs(plot_dir)				
		except: pass
		gen_plot(os.path.join(plot_dir, 'trend.svg'), trend, 50, 10)
		gen_plot(os.path.join(plot_dir, 'means.svg'), means, 50, 10)
		for i in range(len(vals)):
			gen_plot(os.path.join(plot_dir, 'spec_{i}.svg'.format(i = i)), vals[i], 50, 10)

		output.append(outrow)

	print(json.dumps(output))

def process_db(name, filename, dbname, user, password = None):
	'''
	Load data into a database table
	'''

	rows = []
	firstwl = -1
	header = None

	# Load the data into the sets object. TODO: Memory intensive.
	with open(filename, 'rU') as f:
		for line in f:
			if not header:
				# Get the header. Will be used for locating elements in the data rows.
				header = list(map(str.strip, line.split(',')))
				# Find the first wl index
				for i in range(len(header)):
					try:
						int(header[i])
						firstwl = i
						break
					except: pass
			else:
				row = list(map(str.strip, line.split(',')))
				note = row[header.index('note')]
				if note == 'discard':
					continue
				row0 = [name, int(row[header.index('set')]), row[header.index('label')], note, row[header.index('date')], row[header.index('time')]]
				row0.extend(list(map(float, row[firstwl:])))
				rows.append(row0)

	import psycopg2

	db = psycopg2.connect('dbname={dbname} user={user}'.format(dbname = dbname, user = user))
	cur = db.cursor()

	# Create db if necessary
	try:
		wls = ', '.join(['wl_{wl} double precision'.format(wl = wl) for wl in header[firstwl:]])
		cur.execute('create table asd_spectra (id serial primary key, name varchar, setid integer, label varchar, note varchar, date varchar, time varchar, ' + wls + ')')
		db.commit()
	except Exception as e:
		print(e)
		db.rollback()

	# Insert
	stmt = 'insert into asd_spectra (name, setid, label, note, date, time, ' + ', '.join(['wl_{wl}'.format(wl = wl) for wl in header[firstwl:]]) + ') values (' + ', '.join(['%s'] * len(rows[0])) + ')'
	for row in rows:
		cur.execute(stmt, row)

	db.commit()

def gen_id():
	import uuid
	return uuid.uuid4().hex

def gen_plot(outfile, data, width, height):
	import matplotlib.pyplot as plt
	plt.figure(num=None, figsize=(width, height), dpi=1, facecolor='w', edgecolor='k')
	plt.axis('off')
	plt.plot(data)
	plt.savefig(outfile, bbox_inches = 'tight')
	plt.close()

if __name__ == '__main__':

	try:
		phase = sys.argv[1]

		if phase == 'phase1':
			dirname = sys.argv[2]
			process_phase1(dirname)
		elif phase == 'phase2':
			filename = sys.argv[2]
			process_phase2(filename)
		elif phase == 'web':
			filename = sys.argv[2]
			outdir = sys.argv[3]
			process_web(filename, outdir)
		elif phase == 'db':
			name, filename, dbname, user = sys.argv[2:6]
			password = None
			try:
				password = sys.argv[6]
			except: pass
			process_db(name, filename, dbname, user, password)
		else:
			usage()
	except Exception as e:
		import traceback
		print(traceback.format_exc())
		usage()