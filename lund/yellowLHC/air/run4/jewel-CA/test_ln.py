#!/usr/bin/env python

import tutils
import ROOT as r
import argparse
import os


def main():

	tutils.setup_basic_root()

	_e = r.TMath.E()

	tc = r.TCanvas('ln1oD', 'ln1oD')
	tc.Divide(2, 2)
	tc.cd(1)
	f1oD = r.TF1('f1oD', 'TMath::Log(1./x[0])', 0, 0.4)
	f1oD.Draw()
	r.gPad.SetGridx()
	r.gPad.SetGridy()

	tc.cd(2)
	fk = r.TF1('fk', 'TMath::Log(x[0])', 0, 0.4)
	fk.Draw()
	r.gPad.SetGridx()
	r.gPad.SetGridy()

	tc.cd(3)
	fk1 = r.TF1('fk1', 'TMath::Log(x[0])', 0.006, 0.05)
	fk1.Draw()
	fk1.SetRange(r.TMath.Power(_e, -5), r.TMath.Power(_e, -3))
	r.gPad.SetGridx()
	r.gPad.SetGridy()

	tc.cd(4)
	fk2 = r.TF1('fk2', 'TMath::Log(x[0])', 0.0007, 0.012)
	fk2.SetRange(r.TMath.Power(_e, -7), r.TMath.Power(_e, -5))
	fk2.Draw()
	r.gPad.SetGridx()
	r.gPad.SetGridy()

	r.gPad.Update()
	tutils.gList.append(f1oD)
	tutils.gList.append(fk)
	tutils.gList.append(fk1)
	tutils.gList.append(fk2)
	tutils.gList.append(tc)

	tc2 = r.TCanvas('lund', 'lund')
	for dz in range(1, 10):
		z = 1. / (dz)
		# z = 1
		# flund = r.TF1('flund_{}'.format(z), 'TMath::Log([0]/x[0]) + [1]', 1, 6)
		# flund = r.TF1('flund_{}'.format(z), 'TMath::Power(TMath::E(), - [0] * x[0]) + [1]', 1, 6)
		# make a graph instead of a function...
		flund = r.TF1('flund_{}'.format(z), '- x[0] - [0]', 1, 6)
		# flund = r.TF1('flund', 'TMath::Log([0] / x[0]) + [1]', 1/0.4, 0)
		flund.SetParameter(0, z)
		flund.SetParameter(1, 0.0)
		if dz == 1:
			flund.Draw()
		else:
			flund.Draw('same')
		tutils.gList.append(flund)
	# r.gPad.SetLogy()
	# r.gPad.SetLogx()
	r.gPad.Update()
	tutils.gList.append(tc2)


if __name__ == '__main__':
	# https://docs.python.org/2/howto/argparse.html
	parser = argparse.ArgumentParser(description='another starter not much more...', prog=os.path.basename(__file__))
	parser.add_argument('-b', '--batch', help='batchmode - do not end with IPython prompt', action='store_true')
	parser.add_argument('-i', '--prompt', help='end with IPython prompt', action='store_true')
	args = parser.parse_args()
	main()
	if not args.batch:
		tutils.run_app()
	if args.prompt:
		tutils.run_app()
