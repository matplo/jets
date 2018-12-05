#!/usr/bin/env python

import pyutils
import sys


class JobInfo(object):
	def __init__(self, s):
		self.s = 'id = ' + s
		self.d = {}
		lines = self.s.split('\n   ')
		# print lines
		for i, l in enumerate(lines):
			splits = l.split(' = ')
			key = splits[0].lstrip(' ').strip('\n').lstrip('  ')
			val = splits[1].replace('\n', '').replace('\t', '')
			self.d[key] = val
		# for k in self.d:
		# 	print '-[' + k + '] = [' + self.d[k] + ']'
		#print self.d['submit_args']
		if len(self.d['id']) > 0:
			self.valid = True
		else:
			self.valid = False
			# print '[w] Unknown job? : ', self.s   # first line?


class JobStats(object):
	def __init__(self, cmnd):
		self.out, self.err = pyutils.call_cmnd(cmnd)
		lines = self.out.split('Job Id: ')
		self.jobs = []
		for l in lines:
			j = JobInfo(l)
			if j.valid:
				self.jobs.append(j)

def main():
	jstats = JobStats('ssh mpubu qstat -f')

	if '-h' in sys.argv:
		for j in jstats.jobs:
			print '[i] got keys...'
			for k in j.d:
				print k
			break
		return

	jr = []
	jr_o = []
	for j in jstats.jobs:
		if j.d['job_state'] == 'R':
			jr.append(j.d['submit_args'].split(' ')[-1])
			jr_o.append(j.d['Output_Path'])
	print '[i] runnin...', len(jr)
	for i, j in enumerate(jr):
		print j
		print 'output path:', jr_o[i]

	jq = []
	for j in jstats.jobs:
		if j.d['job_state'] == 'Q':
			jq.append(j.d['submit_args'].split(' ')[-1])
	print '[i] queued...', len(jq)
	for j in jq:
		print j

if __name__ == '__main__':
	main()
