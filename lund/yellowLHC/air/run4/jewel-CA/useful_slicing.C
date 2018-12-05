void useful_slicing()
{
TFile *_file0 = TFile::Open("vac_80/subjets_ca_sjR10_R0.4_A2_r0.1_sjA1_sdzcut0.1_sdbeta0_sdr0.4_maxEta3_minpt80_maxpt10000.root");
// TBrowser *s = new TBrowser();
TTree *jt = (TTree*)_file0->Get("jt");
TCanvas *tc = new TCanvas();
TProfile *profdR = new TProfile("profdR", "profdR;n_{iter};#LT #Delta #GT", 11, 0, 11);
jt->Draw("(j_lund_dR):Iteration$>>profdR", "abs(j_jeta) < 2. && j_jpt > 80 && j_jpt < 120", "colz");

TCanvas *tc = new TCanvas();
TProfile *proft = new TProfile("proft", "proft;p_{T,a}+p_{T,b};#LT t #GT = #LT 1/(k_{T} #Delta) #GT = #LT 1/(p_{T,b} #Delta^{2}) #GT", 12*3, 0, 120);
jt->Draw("1/(j_lund_dR*j_lund_dR*j_lund_pt1):(j_lund_pt1 + j_lund_pt2)>>proft", "abs(j_jeta) < 2. && j_jpt > 80 && j_jpt < 120", "colz");
gPad->SetLogy();
TCanvas *tc = new TCanvas();
TProfile *proft2 = new TProfile("proft2", "proft2;ln(1/#Delta);#LT t #GT = #LT 1/(k_{T} #Delta) #GT = #LT 1/(p_{T,b} #Delta^{2}) #GT", 10*3, -0.9, 6);
jt->Draw("1/(j_lund_dR*j_lund_dR*j_lund_pt1):j_lund_log1odr>>proft2", "abs(j_jeta) < 2. && j_jpt > 80 && j_jpt < 120 ", "colz");
gPad->SetLogy();
TCanvas *tc = new TCanvas();
TProfile *proft3 = new TProfile("proft3", "proft3;n_{split};#LT t #GT = #LT 1/(k_{T} #Delta) #GT = #LT 1/(p_{T,b} #Delta^{2}) #GT", 11, 0, 11);
jt->Draw("1/(j_lund_dR*j_lund_dR*j_lund_pt2):Iteration$>>proft3", "abs(j_jeta) < 2. && j_jpt > 80 && j_jpt < 120 ", "colz");
gPad->SetLogy();
TCanvas *tc = new TCanvas();
TProfile *proft4 = new TProfile("proft4", "proft4;z;#LT t #GT = #LT 1/(k_{T} #Delta) #GT = #LT 1/(p_{T,b} #Delta^{2}) #GT", 20, 0, 0.5);
jt->Draw("1/(j_lund_dR*j_lund_dR*j_lund_pt2):j_lund_z>>proft4", "abs(j_jeta) < 2. && j_jpt > 80 && j_jpt < 120 ", "colz");
gPad->SetLogy();
}
