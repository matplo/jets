#!/usr/bin/env python

# import tutils
import ROOT as r
import os
import dlist
import sys


def flatten_2D(h, zmax=None):
	if zmax is None:
		zmax = h.GetMaximum()
	for ibx in range(1, h.GetNbinsX()):
		for iby in range(1, h.GetNbinsY()):
			if h.GetBinContent(ibx, iby) > 0:
				h.SetBinContent(ibx, iby, zmax)
	h.SetContour(2)
	h.SetLineColor(2)


def med_vac(hlund2D, fname, foutname):
	fout = r.TFile(foutname, 'update')
	fout.cd()
	hmed_diff = hlund2D.Clone("{}_med_vac".format(hlund2D.GetName()))
	hmed_ratio = hlund2D.Clone("{}_med_div_vac".format(hlund2D.GetName()))
	fvac_name = fname.replace("med_", "vac_")
	fvac = r.TFile(fvac_name)
	hvac = fvac.Get(hlund2D.GetName())
	zmin_diff = 0
	for ix in range(1, hmed_diff.GetXaxis().GetNbins() + 1):
		for iy in range(1, hmed_diff.GetYaxis().GetNbins() + 1):
			v = hlund2D.GetBinContent(ix, iy) - hvac.GetBinContent(ix, iy)
			if zmin_diff > v * 1.0:
				zmin_diff = v * 1.0
			hmed_diff.SetBinContent(ix, iy, v * 1.0)
			try:
				v = hlund2D.GetBinContent(ix, iy) / hvac.GetBinContent(ix, iy)
				hmed_ratio.SetBinContent(ix, iy, v * 1.0)
			except ZeroDivisionError:
				hmed_ratio.SetBinContent(ix, iy, 1.0)
	for ix in range(1, hmed_diff.GetXaxis().GetNbins() + 1):
		for iy in range(1, hmed_diff.GetYaxis().GetNbins() + 1):
			if hlund2D.GetBinContent(ix, iy) == 0 and hvac.GetBinContent(ix, iy) == 0:
				# v = hmed_diff.GetMinimum()
				v = zmin_diff
				if v <= 0:
					hmed_diff.SetBinContent(ix, iy, 1.2 * v) # just to make white zeros
				else:
					hmed_diff.SetBinContent(ix, iy, 0.8 * v) # just to make white zeros
				# hmed_diff.SetBinContent(ix, iy, zmin_diff * 0.9)
				hmed_ratio.SetBinContent(ix, iy, 1.0)
	# hproj1_vac = fvac.Get("hlund_proj1")
	# hproj2_vac = fvac.Get("hlund_proj2")
	fvac.Close()

	fout.cd()
	if hmed_diff:
		hmed_diff.Write()
		hprofx_diff = hmed_diff.ProfileX()
		hprofx_diff.Write()
		_hprofy_diff = hmed_diff.ProfileY()
		_hprofy_diff.Write()
		hprofy_diff = dlist.h_to_graph(_hprofy_diff, drop_zero_entries=False, xerror=True, transpose=True)
		hprofy_diff.Write()
	if hmed_ratio:
		hmed_ratio.Write()
		hprofx_ratio = hmed_ratio.ProfileX()
		hprofx_ratio.Write()
		_hprofy_ratio = hmed_ratio.ProfileY()
		_hprofy_ratio.Write()
		hprofy_ratio = dlist.h_to_graph(_hprofy_ratio, drop_zero_entries=False, xerror=True, transpose=True)
		hprofy_ratio.Write()
	hmed_diff.Write()
	hmed_ratio.Write()
	fout.Purge()
	fout.Close()

def med_vac_1d(hlund2D, fname, foutname):
	fout = r.TFile(foutname, 'update')
	fout.cd()
	hmed_diff = hlund2D.Clone("{}_med_vac".format(hlund2D.GetName()))
	hmed_diff_relat = hlund2D.Clone("{}_med_vac_relat".format(hlund2D.GetName()))
	hmed_ratio = hlund2D.Clone("{}_med_div_vac".format(hlund2D.GetName()))
	fvac_name = fname.replace("med_", "vac_")
	fvac = r.TFile(fvac_name)
	hvac = fvac.Get(hlund2D.GetName())

	hmed_diff.Add(hvac, -1)
	hmed_diff_relat.Add(hvac, -1)
	hmed_diff_relat.Divide(hvac)
	hmed_ratio.Divide(hvac)

	fvac.Close()

	fout.cd()
	hmed_diff.Write()
	hmed_diff_relat.Write()
	hmed_ratio.Write()
	fout.Purge()
	fout.Close()

def make_diffs(fname):
	foutname = fname.replace('.root', '_diff.root')
	try:
		os.unlink(foutname)
	except OSError:
		pass
	print '[i] diffs for', fname
	print '[i] output:', foutname
	if 'med_' in os.path.dirname(fname):
		f = r.TFile(fname)
		klist = f.GetListOfKeys()
		for k in klist:
			h = f.Get(k.GetName())
			if 'hlund2D' not in h.GetName():
				continue
			print 'make diff for', h.GetName()
			if h.InheritsFrom('TH2'):
				med_vac(h, fname, foutname)
			if h.InheritsFrom('TH1'):
				med_vac_1d(h, fname, foutname)
		f.Close()


def file_draw_lund(fname, ptmin, ptmax):
	print fname, ptmin, ptmax

	f = r.TFile(fname)
	jt = f.Get("jt")

	scond = "abs(j_jeta) < 2. && j_jpt > {} && j_jpt < {}".format(ptmin, ptmax)
	jt.Draw("j_lund_logzdr:j_lund_log1odr>>htmp2D", scond, "e")
	htmp = r.gDirectory.Get("htmp2D")
	xmin = htmp.GetXaxis().GetXmin()
	xmax = htmp.GetXaxis().GetXmax()
	ymin = htmp.GetYaxis().GetXmin()
	ymax = htmp.GetYaxis().GetXmax()

	xmin = 0.9
	xmax = 6.
	# xmax = 10.

	ymin = -8
	# ymin = -10
	ymax = 0.

	nbins = 100

	jt.Draw("j_jpt>>hnjets({}, {}, {})".format(nbins, ptmin, ptmax), scond, "e")
	hnjets = r.gDirectory.Get("hnjets")
	njets = hnjets.Integral()
	print '[i] number of jets', njets

	st = 'j_lund_logzdr:j_lund_log1odr>>hlund2D({}, {}, {}, {}, {}, {})'.format(nbins, xmin, xmax, nbins, ymin, ymax)
	bwx = (xmax - xmin) / (nbins * 1.0)
	bwy = (ymax - ymin) / (nbins * 1.0)
	jt.Draw(st, scond, "e")
	hlund2D = r.gDirectory.Get("hlund2D")
	hlund2D.GetXaxis().SetTitle("j_lund_log1odr")
	hlund2D.GetYaxis().SetTitle("j_lund_logzdr")
	hlund2D.Sumw2()
	hlund2D.Scale(1. / njets / (bwx * bwy))
	hlund2D_flat = hlund2D.Clone("hlund2D_flat")
	flatten_2D(hlund2D_flat, 0.01)

	st = 'j_lund_logzdr:j_lund_log1odr>>hlund2D_dsjL({}, {}, {}, {}, {}, {})'.format(nbins, xmin, xmax, nbins, ymin, ymax)
	scond_dsj = " && j_sj_delta < 0.2 && j_sj_dR < 20 && j_sd_dR > 0.1"
	jt.Draw(st, scond + scond_dsj, "e")
	hlund2D_dsjL = r.gDirectory.Get("hlund2D_dsjL")
	hlund2D_dsjL.Sumw2()
	hlund2D_dsjL.Scale(1. / njets / (bwx * bwy))

	st = 'j_lund_logzdr:j_lund_log1odr>>hlund2D_dsjH({}, {}, {}, {}, {}, {})'.format(nbins, xmin, xmax, nbins, ymin, ymax)
	scond_dsj = " && j_sj_delta > 0.8 && j_sj_dR < 20 && j_sd_dR > 0.1"
	jt.Draw(st, scond + scond_dsj, "e")
	hlund2D_dsjH = r.gDirectory.Get("hlund2D_dsjH")
	hlund2D_dsjH.Sumw2()
	hlund2D_dsjH.Scale(1. / njets / (bwx * bwy))

	st = 'TMath::Log(j_sd_zg * j_sd_dR):TMath::Log(1./j_sd_dR)>>hlund2D_SDzg({}, {}, {}, {}, {}, {})'.format(nbins, xmin, xmax, nbins, ymin, ymax)
	# st = 'TMath::Log(j_sd_zg * j_sd_dR):TMath::Log(1./j_sd_dR)>>hlund2D_SDzg'
	scond_dsj = " && j_sd_zg > 0.1 && j_sd_dR < 0.4"
	jt.Draw(st, scond + scond_dsj, "e")
	hlund2D_SDzg = r.gDirectory.Get("hlund2D_SDzg")
	hlund2D_SDzg.Sumw2()
	hlund2D_SDzg.Scale(1. / njets / (bwx * bwy))
	hlund2D_SDzg_flat = hlund2D_SDzg.Clone("hlund2D_SDzg_flat")
	flatten_2D(hlund2D_SDzg_flat)

	st = 'j_lund_logzdr:j_lund_log1odr>>hlund2D_zg({}, {}, {}, {}, {}, {})'.format(nbins, xmin, xmax, nbins, ymin, ymax)
	scond_dsj = " && j_lund_z > 0.1 && j_lund_dR < 0.4"
	jt.Draw(st, scond + scond_dsj, "e")
	hlund2D_zg = r.gDirectory.Get("hlund2D_zg")
	hlund2D_zg.Sumw2()
	hlund2D_zg.Scale(1. / njets / (bwx * bwy))
	hlund2D_zg_flat = hlund2D_zg.Clone("hlund2D_zg_flat")
	flatten_2D(hlund2D_zg_flat)

	st = 'j_lund_logzdr:j_lund_log1odr>>hlund2D_zgH({}, {}, {}, {}, {}, {})'.format(nbins, xmin, xmax, nbins, ymin, ymax)
	scond_dsj = " && j_lund_z > 0.25 && j_lund_dR < 0.1"
	jt.Draw(st, scond + scond_dsj, "e")
	hlund2D_zgH = r.gDirectory.Get("hlund2D_zgH")
	hlund2D_zgH.Sumw2()
	hlund2D_zgH.Scale(1. / njets / (bwx * bwy))
	hlund2D_zgH_flat = hlund2D_zgH.Clone("hlund2D_zgH_flat")
	flatten_2D(hlund2D_zgH_flat)

	st = 'j_lund_logzdr:j_lund_log1odr>>hlund2D_zgHdRH({}, {}, {}, {}, {}, {})'.format(nbins, xmin, xmax, nbins, ymin, ymax)
	scond_dsj = " && j_lund_z > 0.25 && j_lund_dR > 0.2"
	jt.Draw(st, scond + scond_dsj, "e")
	hlund2D_zgHdRH = r.gDirectory.Get("hlund2D_zgHdRH")
	hlund2D_zgHdRH.Sumw2()
	hlund2D_zgHdRH.Scale(1. / njets / (bwx * bwy))
	hlund2D_zgHdRH_flat = hlund2D_zgHdRH.Clone("hlund2D_zgHdRH_flat")
	flatten_2D(hlund2D_zgHdRH_flat)

	# formation time cuts:
	# TMath::Log(1/(pT x t))
	st = 'j_lund_logzdr:j_lund_log1odr>>hlund2Dtf0({}, {}, {}, {}, {}, {})'.format(nbins, xmin, xmax, nbins, ymin, ymax)
	bwx = (xmax - xmin) / (nbins * 1.0)
	bwy = (ymax - ymin) / (nbins * 1.0)
	scond_tf0 = " && (j_lund_logzdr > j_lund_log1odr + TMath::Log(1e-2)) && (j_lund_logzdr < j_lund_log1odr + TMath::Log(1e-1))"
	jt.Draw(st, scond + scond_tf0, "e")
	hlund2Dtf0 = r.gDirectory.Get("hlund2Dtf0")
	hlund2Dtf0.GetXaxis().SetTitle("j_lund_log1odr")
	hlund2Dtf0.GetYaxis().SetTitle("j_lund_logzdr")
	hlund2Dtf0.Sumw2()
	hlund2Dtf0.Scale(1. / njets / (bwx * bwy))
	hlund2Dtf0_flat = hlund2Dtf0.Clone("hlund2Dtf0_flat")
	flatten_2D(hlund2Dtf0_flat)

	st = 'j_lund_logzdr:j_lund_log1odr>>hlund2Dtf1({}, {}, {}, {}, {}, {})'.format(nbins, xmin, xmax, nbins, ymin, ymax)
	bwx = (xmax - xmin) / (nbins * 1.0)
	bwy = (ymax - ymin) / (nbins * 1.0)
	scond_tf1 = " && (j_lund_logzdr > j_lund_log1odr + TMath::Log(1e-3)) && (j_lund_logzdr < j_lund_log1odr + TMath::Log(1e-2))"
	jt.Draw(st, scond + scond_tf1, "e")
	hlund2Dtf1 = r.gDirectory.Get("hlund2Dtf1")
	hlund2Dtf1.GetXaxis().SetTitle("j_lund_log1odr")
	hlund2Dtf1.GetYaxis().SetTitle("j_lund_logzdr")
	hlund2Dtf1.Sumw2()
	hlund2Dtf1.Scale(1. / njets / (bwx * bwy))
	hlund2Dtf1_flat = hlund2Dtf1.Clone("hlund2Dtf1_flat")
	flatten_2D(hlund2Dtf1_flat)

	st = 'j_lund_logzdr:j_lund_log1odr>>hlund2Dtf2({}, {}, {}, {}, {}, {})'.format(nbins, xmin, xmax, nbins, ymin, ymax)
	bwx = (xmax - xmin) / (nbins * 1.0)
	bwy = (ymax - ymin) / (nbins * 1.0)
	scond_tf2 = " && (j_lund_logzdr > j_lund_log1odr + TMath::Log(1e-4)) && (j_lund_logzdr < j_lund_log1odr + TMath::Log(1e-3))"
	jt.Draw(st, scond + scond_tf2, "e")
	hlund2Dtf2 = r.gDirectory.Get("hlund2Dtf2")
	hlund2Dtf2.GetXaxis().SetTitle("j_lund_log1odr")
	hlund2Dtf2.GetYaxis().SetTitle("j_lund_logzdr")
	hlund2Dtf2.Sumw2()
	hlund2Dtf2.Scale(1. / njets / (bwx * bwy))
	hlund2Dtf2_flat = hlund2Dtf2.Clone("hlund2Dtf2_flat")
	flatten_2D(hlund2Dtf2_flat)

	st = 'j_lund_logzdr:j_lund_log1odr>>hlund2Dtf3({}, {}, {}, {}, {}, {})'.format(nbins, xmin, xmax, nbins, ymin, ymax)
	bwx = (xmax - xmin) / (nbins * 1.0)
	bwy = (ymax - ymin) / (nbins * 1.0)
	scond_tf3 = " && (j_lund_logzdr > j_lund_log1odr + TMath::Log(1e-5)) && (j_lund_logzdr < j_lund_log1odr + TMath::Log(1e-4))"
	jt.Draw(st, scond + scond_tf3, "e")
	hlund2Dtf3 = r.gDirectory.Get("hlund2Dtf3")
	hlund2Dtf3.GetXaxis().SetTitle("j_lund_log1odr")
	hlund2Dtf3.GetYaxis().SetTitle("j_lund_logzdr")
	hlund2Dtf3.Sumw2()
	hlund2Dtf3.Scale(1. / njets / (bwx * bwy))
	hlund2Dtf3_flat = hlund2Dtf3.Clone("hlund2Dtf3_flat")
	flatten_2D(hlund2Dtf3_flat)

	st = 'j_lund_logzdr:j_lund_log1odr>>hlund2Dtf4({}, {}, {}, {}, {}, {})'.format(nbins, xmin, xmax, nbins, ymin, ymax)
	bwx = (xmax - xmin) / (nbins * 1.0)
	bwy = (ymax - ymin) / (nbins * 1.0)
	scond_tf4 = " && (j_lund_logzdr > j_lund_log1odr + TMath::Log(1e-6)) && (j_lund_logzdr < j_lund_log1odr + TMath::Log(1e-5))"
	jt.Draw(st, scond + scond_tf4, "e")
	hlund2Dtf4 = r.gDirectory.Get("hlund2Dtf4")
	hlund2Dtf4.GetXaxis().SetTitle("j_lund_log1odr")
	hlund2Dtf4.GetYaxis().SetTitle("j_lund_logzdr")
	hlund2Dtf4.Sumw2()
	hlund2Dtf4.Scale(1. / njets / (bwx * bwy))
	hlund2Dtf4_flat = hlund2Dtf4.Clone("hlund2Dtf4_flat")
	flatten_2D(hlund2Dtf4_flat)

	# formation time cuts - plot z
	nbins_z = 10
	xmax_z = 0.5
	xmin_z = 0
	bwx_z = (xmax_z - xmin_z) / (nbins_z * 1.0)

	st = 'j_lund_z>>hlund2Dtf0_z({}, {}, {})'.format(nbins_z, 0, 0.5)
	scond_tf0_z = scond_tf0
	jt.Draw(st, scond + scond_tf0_z, "e")
	hlund2Dtf0_z = r.gDirectory.Get("hlund2Dtf0_z")
	hlund2Dtf0_z.GetXaxis().SetTitle("z")
	#hlund2Dtf0_z.Sumw2()
	hlund2Dtf0_z.Scale(1. / njets / (bwx_z))

	st = 'j_lund_z>>hlund2Dtf1_z({}, {}, {})'.format(nbins_z, 0, 0.5)
	scond_tf1_z = scond_tf1
	jt.Draw(st, scond + scond_tf1_z, "e")
	hlund2Dtf1_z = r.gDirectory.Get("hlund2Dtf1_z")
	hlund2Dtf1_z.GetXaxis().SetTitle("z")
	#hlund2Dtf1_z.Sumw2()
	hlund2Dtf1_z.Scale(1. / njets / (bwx_z))

	st = 'j_lund_z>>hlund2Dtf2_z({}, {}, {})'.format(nbins_z, 0, 0.5)
	scond_tf2_z = scond_tf2
	jt.Draw(st, scond + scond_tf2_z, "e")
	hlund2Dtf2_z = r.gDirectory.Get("hlund2Dtf2_z")
	hlund2Dtf2_z.GetXaxis().SetTitle("z")
	#hlund2Dtf2_z.Sumw2()
	hlund2Dtf2_z.Scale(1. / njets / (bwx_z))

	st = 'j_lund_z>>hlund2Dtf3_z({}, {}, {})'.format(nbins_z, 0, 0.5)
	scond_tf3_z = scond_tf3
	jt.Draw(st, scond + scond_tf3_z, "e")
	hlund2Dtf3_z = r.gDirectory.Get("hlund2Dtf3_z")
	hlund2Dtf3_z.GetXaxis().SetTitle("z")
	#hlund2Dtf3_z.Sumw2()
	hlund2Dtf3_z.Scale(1. / njets / (bwx_z))

	st = 'j_lund_z>>hlund2Dtf4_z({}, {}, {})'.format(nbins_z, 0, 0.5)
	scond_tf4_z = scond_tf4
	jt.Draw(st, scond + scond_tf4_z, "e")
	hlund2Dtf4_z = r.gDirectory.Get("hlund2Dtf4_z")
	hlund2Dtf4_z.GetXaxis().SetTitle("z")
	#hlund2Dtf4_z.Sumw2()
	hlund2Dtf4_z.Scale(1. / njets / (bwx_z))

	# profiles
	hprofx = hlund2D.ProfileX()
	_hprofy = hlund2D.ProfileY()
	hprofy = dlist.h_to_graph(_hprofy, drop_zero_entries=False, xerror=True, transpose=True)

	# some projections on x
	# y: j_lund_logzdr
	nbins_proj = nbins / 4.
	st = 'j_lund_log1odr>>hlund_proj1({}, {}, {})'.format(nbins_proj, xmin, xmax)
	scond_proj1 = "abs(j_jeta) < 2. && j_jpt > {} && j_jpt < {} && j_lund_logzdr < {} && j_lund_logzdr > {}".format(ptmin, ptmax, -3., -5.)
	jt.Draw(st, scond_proj1, "e")
	hlund_proj1 = r.gDirectory.Get("hlund_proj1")
	hlund_proj1.GetXaxis().SetTitle("j_lund_log1odr")
	hlund_proj1.GetYaxis().SetTitle("counts")
	hlund_proj1.Sumw2()
	hlund_proj1.Scale(1. / njets / hlund_proj1.GetBinWidth(1))

	st = 'j_lund_log1odr>>hlund_proj2({}, {}, {})'.format(nbins_proj, xmin, xmax)
	scond_proj2 = "abs(j_jeta) < 2. && j_jpt > {} && j_jpt < {} && j_lund_logzdr < {} && j_lund_logzdr > {}".format(ptmin, ptmax, -5., -7.)
	jt.Draw(st, scond_proj2, "e")
	hlund_proj2 = r.gDirectory.Get("hlund_proj2")
	hlund_proj2.GetXaxis().SetTitle("j_lund_log1odr")
	hlund_proj2.GetYaxis().SetTitle("counts")
	hlund_proj2.Sumw2()
	hlund_proj2.Scale(1. / njets / hlund_proj2.GetBinWidth(1))
	# end projections

	foutname = '{}/hout_{}'.format(os.path.dirname(fname), os.path.basename(fname))

	# make a diff and a ratio
	hmed_diff = None
	hmed_ratio = None
	if 'med_' in os.path.dirname(foutname):
		hmed_diff = hlund2D.Clone("med_vac")
		hmed_ratio = hlund2D.Clone("med_div_vac")
		fvac_name = foutname.replace("med_", "vac_")
		fvac = r.TFile(fvac_name)
		hvac = fvac.Get("hlund2D")
		zmin_diff = 0
		for ix in range(1, hmed_diff.GetXaxis().GetNbins() + 1):
			for iy in range(1, hmed_diff.GetYaxis().GetNbins() + 1):
				v = hlund2D.GetBinContent(ix, iy) - hvac.GetBinContent(ix, iy)
				if zmin_diff > v * 1.0:
					zmin_diff = v * 1.0
				hmed_diff.SetBinContent(ix, iy, v * 1.0)
				try:
					v = hlund2D.GetBinContent(ix, iy) / hvac.GetBinContent(ix, iy)
					hmed_ratio.SetBinContent(ix, iy, v * 1.0)
				except ZeroDivisionError:
					hmed_ratio.SetBinContent(ix, iy, 1.0)
		for ix in range(1, hmed_diff.GetXaxis().GetNbins() + 1):
			for iy in range(1, hmed_diff.GetYaxis().GetNbins() + 1):
				if hlund2D.GetBinContent(ix, iy) == 0 and hvac.GetBinContent(ix, iy) == 0:
					v = zmin_diff
					if v <= 0:
						hmed_diff.SetBinContent(ix, iy, 1.2 * v) # just to make white zeros
					else:
						hmed_diff.SetBinContent(ix, iy, 0.8 * v) # just to make white zeros
					hmed_ratio.SetBinContent(ix, iy, 0.0)
		# hproj1_vac = fvac.Get("hlund_proj1")
		# hproj2_vac = fvac.Get("hlund_proj2")
		fvac.Close()

	fout = r.TFile(foutname, "recreate")
	fout.cd()
	htmp.Write()
	hlund2D.Write()
	hlund2D_flat.Write()

	hlund2D_dsjL.Write()
	hlund2D_dsjH.Write()
	hlund2D_SDzg.Write()
	hlund2D_SDzg_flat.Write()
	hlund2D_zg.Write()
	hlund2D_zg_flat.Write()
	hlund2D_zgH.Write()
	hlund2D_zgH_flat.Write()
	hlund2D_zgHdRH.Write()
	hlund2D_zgHdRH_flat.Write()

	hlund2Dtf0.Write()
	hlund2Dtf1.Write()
	hlund2Dtf2.Write()
	hlund2Dtf3.Write()
	hlund2Dtf4.Write()

	hlund2Dtf0_flat.Write()
	hlund2Dtf1_flat.Write()
	hlund2Dtf2_flat.Write()
	hlund2Dtf3_flat.Write()
	hlund2Dtf4_flat.Write()


	hlund2Dtf0_z.Write()
	hlund2Dtf1_z.Write()
	hlund2Dtf2_z.Write()
	hlund2Dtf3_z.Write()
	hlund2Dtf4_z.Write()

	hnjets.Write()
	hprofy.Write()
	_hprofy.Write()
	hprofx.Write()
	if hmed_diff:
		hmed_diff.Write()
		hprofx_diff = hmed_diff.ProfileX()
		hprofx_diff.Write()
		_hprofy_diff = hmed_diff.ProfileY()
		_hprofy_diff.Write()
		hprofy_diff = dlist.h_to_graph(_hprofy_diff, drop_zero_entries=False, xerror=True, transpose=True)
		hprofy_diff.Write()
	if hmed_ratio:
		hmed_ratio.Write()
		hprofx_ratio = hmed_ratio.ProfileX()
		hprofx_ratio.Write()
		_hprofy_ratio = hmed_ratio.ProfileY()
		_hprofy_ratio.Write()
		hprofy_ratio = dlist.h_to_graph(_hprofy_ratio, drop_zero_entries=False, xerror=True, transpose=True)
		hprofy_ratio.Write()
	hlund_proj1.Write()
	hlund_proj2.Write()
	fout.Close()


def lund_draw():
	fname = "vac_80/subjets_ca_sjR10_R0.4_A2_r0.1_sjA1_sdzcut0.1_sdbeta0_sdr0.4_maxEta3_minpt80_maxpt10000.root"
	file_draw_lund(fname, 80, 120)
	fname = "med_80/subjets_ca_sjR10_R0.4_A2_r0.1_sjA1_sdzcut0.1_sdbeta0_sdr0.4_maxEta3_minpt80_maxpt10000.root"
	file_draw_lund(fname, 80, 120)
	fname = "vac_200/subjets_ca_sjR10_R0.4_A2_r0.1_sjA1_sdzcut0.1_sdbeta0_sdr0.4_maxEta3_minpt200_maxpt10000.root"
	file_draw_lund(fname, 200, 250)
	fname = "med_200/subjets_ca_sjR10_R0.4_A2_r0.1_sjA1_sdzcut0.1_sdbeta0_sdr0.4_maxEta3_minpt200_maxpt10000.root"
	file_draw_lund(fname, 200, 250)


def diffs():
	fname = "med_80/hout_subjets_ca_sjR10_R0.4_A2_r0.1_sjA1_sdzcut0.1_sdbeta0_sdr0.4_maxEta3_minpt80_maxpt10000.root"
	make_diffs(fname)
	fname = "med_200/hout_subjets_ca_sjR10_R0.4_A2_r0.1_sjA1_sdzcut0.1_sdbeta0_sdr0.4_maxEta3_minpt200_maxpt10000.root"
	make_diffs(fname)


if __name__ == '__main__':
	if '--lund' in sys.argv:
		lund_draw()
	if '--diffs' in sys.argv:
		diffs()
