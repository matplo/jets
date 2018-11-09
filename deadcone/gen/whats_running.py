#!/usr/bin/env python

import pyutils

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
	jr = []
	for j in jstats.jobs:
		if j.d['job_state'] == 'R':
			jr.append(j.d['submit_args'].split(' ')[-1])
	print '[i] runnin...', len(jr)
	for j in jr:
		print j

	jq = []
	for j in jstats.jobs:
		if j.d['job_state'] == 'Q':
			jq.append(j.d['submit_args'].split(' ')[-1])
	print '[i] queued...', len(jq)
	for j in jq:
		print j

if __name__ == '__main__':
	main()
