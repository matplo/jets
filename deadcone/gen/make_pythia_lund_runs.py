#!/usr/bin/env python

import os

def grab_module_environment():
	rets = []
	module_home = os.environ.get('MODULESHOME')
	init_scipt = os.path.join(module_home, 'init/bash')
	if os.path.isfile(init_scipt):
		rets.append('. {}'.format(init_scipt))
	module_path = os.environ.get('MODULEPATH')
	splits = module_path.split(':')
	for sp in splits:
		rets.append('module use {}'.format(sp))
	module_loaded = os.environ.get('LOADEDMODULES')
	splits = module_loaded.split(':')
	for sp in splits:
		rets.append('module load {}'.format(sp))
	return '\n'.join(rets)


class PythiaSJLund(object):
	jetty_command = 'jetty_subjets_exe --ca-task --nev={nev} --pTHatMin={minpt} --jptcut={minpt} Beams:eCM={ecm} {process} {level} {seed}'
	def __init__(self, outputdir, process_short_name):
		self.outputdir_base = outputdir
		self.seed = 'Random:setSeed=on Random:seed=0'
		self.nev = 10000
		self.pthatmin = 20.
		self.jptcut = self.pthatmin
		self.ecm = 13000
		self.process_short_name = process_short_name

	def make_job(self):
		self.outputdir = os.path.join(self.outputdir_base, self.process_short_name, 'ecm_' + str(self.ecm), 'pthatmin_' + str(self.pthatmin))
		print self.outputdir


def main():
	sjl = PythiaSJLund('/tmp', 'test')
	sjl.make_job()
	print grab_module_environment()

if __name__ == '__main__':
	main()
