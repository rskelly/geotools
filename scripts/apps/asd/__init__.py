import json
import psycopg2

def get_names(db):
	names = []
	cur = db.cursor()
	for row in cur.execute('select distinct name from asd_spectra'):
		names.append(row[0])
	return names

def run(web):
	input = web.input(data = '')
	db = None
	if data:
		db = psycopg2.connect("dbname=web user=rob")
	if data == 'names':
		return json.dumps(get_names(db))

	return 'running ' + input.var