#!/usr/bin/env python

import os
import sys


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
	def _setting(self, key, default_value=None):
		try:
			return self.kwargs[key]
		except KeyError:
			return default_value
		return None

	def __init__(self, outputdir, **kwargs):
		self.kwargs = kwargs
		self.outputdir_base = outputdir
		# self.seed = 'Random:setSeed=on Random:seed=0'
		self.seed = self._setting('seed', '--time-seed')
		self.pyseed = ''
		if self.seed == '--time-seed':
			self.pyseed = self.seed
		else:
			if self.seed == '--no-seed':
				self.pyseed = ''
			else:
				self.pyseed = 'Random:setSeed=on Random:seed={}'.format(self.seed)
		self.nev = self._setting('nev', 10000)
		self.pthatmin = self._setting('pthatmin', 20.)
		self.jptcut = self.pthatmin
		self.ecm = self._setting('ecm', 13000)
		self.R = self._setting('R', 0.4)
		self.level = self._setting('level', 'hadron')
		self.pylevel = ''
		if self.level == 'parton':
			self.pylevel = 'HadronLevel:all=off'
		else:
			self.pylevel = ''
		self.process_short_name = None
		process_short_name = self._setting('process', None)
		if process_short_name in ['hardQCD', 'hardQCDlf', 'hardQCDuds', 'hardQCDgluons', 'hardQCDquarks', 'hardQCDbeauty', 'hardQCDcharm']:
			self.process_short_name = process_short_name
		self.dir_seed = 0
		self.jetty_command = 'jetty_subjets_exe --ca-task --R={} --nev={} --pTHatMin={} --jptcut={} --eCM={} --{} {} {} {}'.format(self.R, self.nev, self.pthatmin, self.jptcut, self.ecm, self.process_short_name, self.pyseed, self.pylevel, self._setting('extra', ''))
		self.outputfname_base = self._setting('fname', 'job.sh')

	def get_output_dir(self):
		outputdir_base = os.path.join(self.outputdir_base, 'ecm_' + str(self.ecm), self.process_short_name, self.level, 'pthatmin_' + str(self.pthatmin))
		if self.seed == '--no-seed':
			return outputdir_base
		if self.seed == '--time-seed':
			self.dir_seed = 0
			outputdir = os.path.join(outputdir_base, 's' + str(self.dir_seed))
			self.outputfname = os.path.join(outputdir, self.outputfname_base)
			# while os.path.exists(outputdir):
			while os.path.exists(self.outputfname):
				self.dir_seed = self.dir_seed + 1
				outputdir = os.path.join(outputdir_base, 's' + str(self.dir_seed))
				self.outputfname = os.path.join(outputdir, self.outputfname_base)
			return outputdir
		self.dir_seed = self.seed
		outputdir = os.path.join(outputdir_base, 's' + str(self.dir_seed))
		return outputdir

	def make_job(self):
		self.outputdir = self.get_output_dir()
		if not self.process_short_name:
			print '[e] unable to make the job with process:', self.process_short_name
		try:
			os.makedirs(self.outputdir)
		except OSError as e:
			if e.errno != 17:
				print >> sys.stderr, e
				return
		modenv = grab_module_environment()
		try:
			with open(self.outputfname, 'w') as f:
				f.write('#!/bin/bash\n')
				f.writelines([modenv, '\n', 'cd {}\n'.format(self.outputdir), self.jetty_command])
			print '[i] written:', self.outputfname
		except OSError as e:
			print >> sys.stderr, e
			return


def main():
	for proc in ['hardQCDbeauty', 'hardQCDuds', 'hardQCDcharm', 'hardQCDgluons']:
		for pthm in [200, 20, 40, 80, 100]:
			extra_s = ''
			extra_s = 'PartonLevel:MPI=off PartonLevel:ISR=off'
			sjl = PythiaSJLund('./', nev=100000, pthatmin=pthm, process=proc, fname='jobR07.sh', level='parton', extra=extra_s, R=0.7)
			sjl.make_job()
			sjl = PythiaSJLund('./', nev=100000, pthatmin=pthm, process=proc, fname='jobR04.sh', level='parton', extra=extra_s)
			sjl.make_job()

			extra_s = ''
			if 'charm' in proc:
				extra_s = '411:maydecay=no 421:maydecay=no'
			if 'beauty' in proc:
				extra_s = '511:maydecay=no'
			sjl = PythiaSJLund('./', nev=100000, pthatmin=pthm, process=proc, fname='jobR07.sh', level='hadron', extra=extra_s, R=0.7)
			sjl.make_job()
			sjl = PythiaSJLund('./', nev=100000, pthatmin=pthm, process=proc, fname='jobR04.sh', level='hadron', extra=extra_s)
			sjl.make_job()


if __name__ == '__main__':
	main()
