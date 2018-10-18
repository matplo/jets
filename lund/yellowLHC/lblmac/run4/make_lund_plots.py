#!/usr/bin/env python

import tutils
import ROOT as t
import IPython
import argparse
import os
import dlist

def main(args):

	r = int(args.radius)
	if args.radius < 1:
		r = int(args.radius * 10.)
	sr = '0' + '{}'.format(r)[0]
	hnspname = 'fHLundIterative_fuJet_R{}_Sum'.format(sr)
	stype = 'full'
	if args.charged:
		hnspname = 'fHLundIterative_chJet_R{}_Sum'.format(sr)
		stype = 'charged'
	fin = t.TFile(args.filename)
	if not fin.IsOpen():
		print '[e] unable to open', args.filename
		return None

	hn = fin.Get(hnspname)
	print hnspname,'at',hn
	sfoutname = 'proj_r{}_{}.root'.format(sr, stype)
	fout = t.TFile(sfoutname, 'recreate')
	if fout.IsOpen():
		ptaxis = hn.GetAxis(2)
		# hn.GetAxis(1).SetRange(hn.GetAxis(1).FindBin(1.), hn.GetAxis(1).FindBin(5.))
		for b in range(1, ptaxis.GetNbins()):
			ptmin = ptaxis.GetBinLowEdge(b)
			ptmax = ptaxis.GetBinLowEdge(b) + ptaxis.GetBinWidth(b)
			print b, ptmin, ptmax
			ptaxis.SetRange(b, b)
			prj = hn.Projection(0, 1, 'EAO')
			prj.SetName('lund_r{}_{}_pt_{}_{}'.format(sr, stype, int(ptmin), int(ptmax)))
			# prj.SetTitle('lund_r{}_{}_pt_{}_{};{};{}'.format(sr, stype, int(ptmin), int(ptmax), hn.GetAxis(0).GetTitle(), hn.GetAxis(1).GetTitle()))
			# prj.Smooth(1)
			prj.Write()
			px = prj.ProfileX()
			py = prj.ProfileY()
			px.Write()
			py.Write()
			pyT = dlist.h_to_graph(py, drop_zero_entries=True, xerror=True, transpose=True)
			pyT.Write(py.GetName() + '_T')
		fout.Close()

	# tutils.setup_basic_root()
	# tutils.gList.append(h)

if __name__=="__main__":
	# https://docs.python.org/2/howto/argparse.html
	parser = argparse.ArgumentParser(description='make lund plots in pT', prog=os.path.basename(__file__))
	# parser.add_argument('-i', '--prompt', help='end with IPython prompt', action='store_true')
	parser.add_argument('-r', '--radius', help='jet radius', type=float, default=0.4)
	parser.add_argument('-b', '--rebin', help='rebin projections', type=int, default=0)
	parser.add_argument('--charged', help='use charged jets', action='store_true')
	parser.add_argument('filename', help='input root file with THnSparse', type=str, default='merged.root')
	args = parser.parse_args()
	print args
	main(args)
	# if args.prompt:
	# 	tutils.run_app()
