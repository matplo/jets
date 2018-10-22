#!/usr/bin/env python

# import tutils
import ROOT as r
import os
import dlist
import sys

from lund_draw import file_draw_lund
from lund_draw import make_diffs

def lund_draw():
	fnames = [  "vac_80/subjets_ca_sjR10_R0.4_A2_r0.1_sjA1_sdzcut0.1_sdbeta0_sdr0.4_maxEta3_minpt80_maxpt10000.root",
				"med_80/subjets_ca_sjR10_R0.4_A2_r0.1_sjA1_sdzcut0.1_sdbeta0_sdr0.4_maxEta3_minpt80_maxpt10000.root",
				"vac_200/subjets_ca_sjR10_R0.4_A2_r0.1_sjA1_sdzcut0.1_sdbeta0_sdr0.4_maxEta3_minpt200_maxpt10000.root",
				"med_200/subjets_ca_sjR10_R0.4_A2_r0.1_sjA1_sdzcut0.1_sdbeta0_sdr0.4_maxEta3_minpt200_maxpt10000.root"]
	pt0 = 80
	pt1 = 120
	for fname in fnames:
		if '_200' in fname:
			pt0 = 200
			pt1 = 250
		file_draw_lund(fname, pt0, pt1)

def lund_draw_extra():
	fnames = []
	iterations = True
	for lundpt in [[0, 1e5]]:
		for fname in fnames:
			file_draw_lund(fname, 80, 120, lundpt[0], lundpt[1], iterations=iterations)

def diffs():
	fref = "vac_80/hout_subjets_ca_sjR10_R0.4_A2_r0.1_sjA1_sdzcut0.1_sdbeta0_sdr0.4_maxEta3_minpt80_maxpt10000.root"
	fname = "med_80/hout_subjets_ca_sjR10_R0.4_A2_r0.1_sjA1_sdzcut0.1_sdbeta0_sdr0.4_maxEta3_minpt80_maxpt10000.root"
	make_diffs(fname, fref)
	fref = "vac_200/hout_subjets_ca_sjR10_R0.4_A2_r0.1_sjA1_sdzcut0.1_sdbeta0_sdr0.4_maxEta3_minpt200_maxpt10000.root"
	fname = "med_200/hout_subjets_ca_sjR10_R0.4_A2_r0.1_sjA1_sdzcut0.1_sdbeta0_sdr0.4_maxEta3_minpt200_maxpt10000.root"
	make_diffs(fname, fref)

if __name__ == '__main__':
	if '--lund' in sys.argv:
		lund_draw()
	if '--diffs' in sys.argv:
		diffs()
