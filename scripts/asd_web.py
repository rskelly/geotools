#!/usr/bin/env python3

import web
import os
import importlib

urls = (
	'/(.*)', 'index'
)

render = web.template.render('/home/rob/Documents/git/contrem/scripts/templates/')
app = web.application(urls, globals())

appsdir = '/home/rob/Documents/git/contrem/scripts/apps/'

def appsmenu():
	apps = [x for x in os.listdir(appsdir) if os.path.isdir(os.path.join(appsdir, x))]
	return str(render.appsmenu(apps))

class index:
	def GET(self, appname):
		if not appname:
			return appsmenu()
		else:
			try:
				_app = importlib.import_module('apps.' + appname)
				return _app.run(web)
			except Exception as e:
				output = 'Failed to load app, ' + appname + ': ' + e.__str__()
				output += appsmenu()
				return output

if __name__ == "__main__":
	app.run()