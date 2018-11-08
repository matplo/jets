#!/usr/bin/env python

import tutils
import ROOT as r
import IPython
import argparse
import os
import draw_utils as du

g = None
i = 0
nfile = 0

def myexec():
	# get event information
	event = r.gPad.GetEvent()
	px    = r.gPad.GetEventX()
	py    = r.gPad.GetEventY()
	# some magic to get the coordinates...
	xd = r.gPad.AbsPixeltoX(px);
	yd = r.gPad.AbsPixeltoY(py);
	x = r.gPad.PadtoX(xd);
	y = r.gPad.PadtoY(yd);
	# left mouse button click
	# print x, y

	global g, i, nfile


	if event == 1:
		D = 1./r.TMath.Exp(x)
		k = r.TMath.Exp(y)
		z = k / D
		print 'ln(1/D)=', x,  'D=', 1./r.TMath.Exp(x), 'ln(k)=', y, 'k=', r.TMath.Exp(y), 'z=', z
		print 'ln(1/D)=', r.TMath.Log(1./D),  'D=', 1./r.TMath.Exp(x), 'ln(k)=', r.TMath.Log(z*D), 'k=', r.TMath.Exp(y), 'z=', z

		if g is None:
			g = r.TGraph()
			tutils.gList.append(g)
		g.SetPoint(i,x,y)
		if i == 0:
			g.Draw("L")
		i = i + 1
		r.gPad.Update();

	if event == 24:
		fout = r.TFile("click_tgraph_{}.root".format(nfile), "recreate");
		g.Write();
		print "dumpding graph... to " , fout.GetName()
		fout.Close();
		g = None
		i = 0
		nfile = nfile + 1
		r.gPad.Update();

	return

class TestCall( r.TObject ):
	def __init__(self):
		r.TObject()
	def callback(self):
		print 'call...'

def line_intersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1]) #Typo was here

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       raise Exception('lines do not intersect')

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y

def log_k_for_d_z(d, z):
	return r.TMath.Log(z * d)

def z_for_log1od_logk(log1od, logk):
	d = 1./r.TMath.Exp(log1od)
	z = r.TMath.Exp(logk) / d
	return z

def z_for_d_logk(d, logk):
	z = r.TMath.Exp(logk) / d
	return z

def k_for_d_z(d, z):
	return z * d

def z_for_d_k(d, logk):
	z = k / d
	return z

def max_logk_for_d(d):
	return log_k_for_d_z(d, 0.5)

def max_d_for_logk(logk):
	d = r.TMath.Exp(logk) / 0.5
	return d

def adjust_pad_margins(_left=0.17, _right=0.05, _top=0.1, _bottom=0.17+0.03):
	du.adjust_pad_margins(_left, _right, _top, _bottom)

def main(draw_z = True, draw_tau = True, draw_theta = True, pdf = False):

	tutils.setup_basic_root(-1)
	# tutils.setup_basic_root()

	tc = r.TCanvas('lund_z_{}_tau_{}_theta_{}'.format(int(draw_z), int(draw_tau), int(draw_theta)),
	               'lund_z_{}_tau_{}_theta_{}'.format(int(draw_z), int(draw_tau), int(draw_theta)),
	               1000, 1000)
	tutils.gList.append(tc)
	# r.gPad.SetGridx()
	# r.gPad.SetGridy()

	r.gPad.AddExec("myexec",'TPython::Exec("myexec()");')

	log1od_low = 0.9
	log1od_hi = 6.1
	# log1od_low = 0.9
	log1od_low = 0.0
	log1od_hi = 8.9
	h = r.TH2D('h_lund_z_{}_tau_{}_theta_{}'.format(int(draw_z), int(draw_tau), int(draw_theta)),
	           'Lund Jet Plane;ln(1/#Delta);ln(#kappa) := ln(z #Delta)', 8, log1od_low, log1od_hi, 8, -8, 0)
	h.SetLineColor(0)
	h.SetMarkerColor(0)
	h.Draw("colz")
	tutils.gList.append(h)

	# fin = r.TFile("jewel/run3/med_200/hout_subjets_ca_sjR10_R0.4_A2_r0.1_sjA1_sdzcut0.1_sdbeta0_sdr0.4_maxEta3_minpt200_maxpt10000.root")
	# hjewel = fin.Get("med_vac")
	fname = "/Volumes/mp128/data/run4/pythia-CA/lightf/hout_subjets_ca_sjR10_R0.4_A2_r0.1_sjA1_sdzcut0.1_sdbeta0_sdr0.4_maxEta3_minpt80_maxpt10000.root"
	fin = r.TFile(fname)
	hjewel = fin.Get("hlund2D_tw")
	#hjewel.Smooth()
	hjewel.SetMaximum(0.25)
	hjewel.SetMinimum(-0.8)
	hjewel.DrawCopy("colz same")
	tutils.gList.append(hjewel)

	adjust_pad_margins(0.1, 0.15, 0.05, 0.1)

	# ic = [0, 1, 2, 3, 4, 6, 7, 8]
	ic = [0, 32, 5, 2, 3, 4, 6, 7, 8, 9]
	ic.reverse()

	# tleg = r.gPad.BuildLegend(0.5,0.3922187,0.8797595,0.8801262)
	tleg = r.gPad.BuildLegend(0.4624248,0.4574132,0.8421844,0.9453207)
	tleg.SetBorderSize(0)

	if draw_z:
		iz = 0
		for z in [0.01, 0.05, 0.1, 0.2, 0.25, 0.3, 0.5, 1.0]:
			iz = iz + 1
			d_low = 1./r.TMath.Exp(log1od_low)
			d_hi = 1./r.TMath.Exp(log1od_hi)
			fname = 'k_z_{}'.format(z)
			# f = r.TF1(fname, 'TMath::Log({} * x[0])'.format(z), d_low, d_hi)
			f = r.TF1(fname, 'TMath::Log({} * 1./TMath::Exp(x[0]))'.format(z), log1od_low, log1od_hi)
			# f.SetTitle('z = {}'.format(z))
			f.SetLineColor(ic[iz])
			tleg.AddEntry(f, 'z = {}'.format(z), 'l')
			f.Draw("same")
			tutils.gList.append(f)


	if draw_tau:
		iz = 0
		# for z in [-10, -9, -8, -7, -6]:
		for z in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]:
			iz = iz + 1
			d_low = 1./r.TMath.Exp(log1od_low)
			d_hi = 1./r.TMath.Exp(log1od_hi)
			fname = 'k_c_{}'.format(z)
			_log1od_hi = log1od_hi
			# print r.TMath.Log(z)
			# print log1od_hi * r.TMath.Log(z)
			maxd = r.TMath.Log(0.5 / z) / 2.
			_log1od_hi = maxd
			# print log1od_hi, r.TMath.Log(z), z, maxd, _log1od_hi
			f = r.TF1(fname, 'x[0] + TMath::Log({})'.format(z), log1od_low, _log1od_hi)
			#f.SetTitle('ln(1/#Delta #times {}) #approx ln(1/#Delta) {:.1f}'.format(z, r.TMath.Log(z)))
			f.SetLineColor(ic[iz])
			f.SetLineStyle(9)
			tleg.AddEntry(f, 'ln(1/#Delta #times {}) #approx ln(1/#Delta) {:.1f}'.format(z, r.TMath.Log(z)), 'l')
			f.Draw("same")
			tutils.gList.append(f)

	if draw_theta:
		iz = 0
		for z in [0.4, 0.3, 0.2, 0.1, 0.05, 0.01]:
			iz = iz + 1
			maxY = r.TMath.Log(0.5 * z)
			l = r.TLine(r.TMath.Log(1/z), -8, r.TMath.Log(1/z), maxY)
			l.SetLineColor(ic[iz])
			l.SetLineStyle(2)
			l.Draw()
			tleg.AddEntry(l, '#Delta = {}'.format(z), 'l')
			tutils.gList.append(l)
			print '#Delta = {}'.format(z)
			print '#line  {}, {}, {}, {},  38, 1, 1'.format(r.TMath.Log(1/z), -8, r.TMath.Log(1/z), maxY)
	r.gPad.Update()
	if pdf:
		tc.Print(tc.GetName()+'.pdf', 'pdf')

	if draw_theta and draw_tau and draw_z:
		fout = r.TFile('curves.root', 'recreate')
		for o in tutils.gList:
			try:
				if o.InheritsFrom('TLine') or o.InheritsFrom('TF1'):
					o.Write()
			except:
				pass
		fout.Purge()
		fout.Close()


if __name__=="__main__":
	# https://docs.python.org/2/howto/argparse.html
	parser = argparse.ArgumentParser(description='another starter not much more...', prog=os.path.basename(__file__))
	parser.add_argument('-b', '--batch', help='batchmode - do not end with IPython prompt', action='store_true')
	parser.add_argument('-i', '--prompt', help='end with IPython prompt', action='store_true')
	args = parser.parse_args()
	main(False, False, False, True)
	main(True, False, False, True)
	main(True, True, False, True)

	main(False, True, False, True)
	main(False, False, True, True)
	main(True, False, True, True)

	main(True, True, True, True)

	if not args.batch:
		tutils.run_app()
	if args.prompt:
		tutils.run_app()
