void file_draw_lund(const char *fname)
{
	TFile *f = new TFile(fname);
	TTree *jt = (TTree*)f->Get("jt");
	jt->Draw("j_lund_logzdr:j_lund_log1odr>>htmp2D", "", "colz");
	TH2 *htmp = gDirectory->Get("htmp2D");
	Double_t xmin = htmp->GetXaxis()->GetXmin();
	Double_t xmax = htmp->GetXaxis()->GetXmax();

	Double_t ymin = htmp->GetYaxis()->GetXmin();
	Double_t ymax = htmp->GetYaxis()->GetXmax();

	cout << xmin << " " << xmax << endl;
	cout << ymin << " " << ymax << endl;

	Int_t nbins = 100;
	TString st = TString::Format("j_lund_logzdr:j_lund_log1odr>>hlund2D(%d, %f, %f, %d, %f %f)", nbins, xmin, xmax, nbins, ymin, ymax);
	jt->Draw(st.Data(), "(j_jpt>80 && j_jpt<120)", "colz");
	TH2 *htmp = (TH2*)gDirectory->Get("hlund2D");
	htmp->GetXaxis()->SetTitle("j_lund_log1odr");
	htmp->GetYaxis()->SetTitle("j_lund_logzdr");

	TString foutname = TString::Format("hout_%s", fname);

	TFile *fout = new TFile(foutname.Data(), "recreate");
	fout->cd();
	htmp->Write();
	fout->Close();
}

void lund_draw()
{
	const char *fname = "subjets_ca_sjR10_R0.4_A2_r0.1_sjA1_sdzcut0.1_sdbeta0_sdr0.4_maxEta3_minpt80_maxpt10000.root";
	file_draw_lund(fname);
}
