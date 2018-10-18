
#include "HepMC/Units.h"
#include "HepMC/WeightContainer.h"
#include "HepMC/CompareGenEvent.h"
#include "HepMC/enable_if.h"
#include "HepMC/Flow.h"
#include "HepMC/GenCrossSection.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenRanges.h"
#include "HepMC/GenVertex.h"
#include "HepMC/HeavyIon.h"
#include "HepMC/HEPEVT_Wrapper.h"
#include "HepMC/HepMCDefs.h"
#include "HepMC/HerwigWrapper.h"
#include "HepMC/IO_AsciiParticles.h"
#include "HepMC/IO_BaseClass.h"
#include "HepMC/IO_Exception.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/IO_HEPEVT.h"
#include "HepMC/IO_HERWIG.h"
#include "HepMC/is_arithmetic.h"
#include "HepMC/IteratorRange.h"
#include "HepMC/PdfInfo.h"
#include "HepMC/Polarization.h"
#include "HepMC/PythiaWrapper6_4.h"
#include "HepMC/PythiaWrapper6_4_WIN32.h"
#include "HepMC/PythiaWrapper.h"
#include "HepMC/SearchVector.h"
#include "HepMC/SimpleVector.h"
#include "HepMC/StreamHelpers.h"
#include "HepMC/StreamInfo.h"
#include "HepMC/TempParticleMap.h"
#include "HepMC/Version.h"

#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "TPDGCode.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "THnSparse.h"  //LUND

#include <iostream>
#include <fstream>
#include <iomanip>
#include <TString.h>
#include <algorithm>  //LUND
#include <vector>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;

//class TProfile;

//______________________________________________________________________________________
class TT{     //wrapper for TT bins settings
   public:
     TT(): triggLow(0),triggHigh(0){ }

     void Set(Int_t lt, Int_t ht){
        triggLow      = lt;
        triggHigh     = ht;
     }

     Int_t L()  const { return triggLow;}
     Int_t H()  const { return triggHigh;}

   private:
     Int_t   triggLow; //low trigger pt range
     Int_t   triggHigh; //high trigger pt range
};

//______________________________________________________________________________________

/*
int is_stable(const HepMC::GenParticle *part) {
  // copied from AliStack::IsStable()	 
  int pdg = abs(part->pdg_id());
  if(pdg>1000000000)return kTRUE;

  const Int_t kNstable = 18;
  Int_t i;
  

  Int_t pdgStable[kNstable] = {
    kGamma,             // Photon
    kElectron,          // Electron
    kMuonPlus,          // Muon 
    kPiPlus,            // Pion
    kKPlus,             // Kaon
    kK0Short,           // K0s
    kK0Long,            // K0l
    kProton,            // Proton 
    kNeutron,           // Neutron
    kLambda0,           // Lambda_0
    kSigmaMinus,        // Sigma Minus
    kSigmaPlus,         // Sigma Plus
    3312,               // Xsi Minus 
    3322,               // Xsi 
    3334,               // Omega
    kNuE,               // Electron Neutrino 
    kNuMu,              // Muon Neutrino
    kNuTau              // Tau Neutrino
  };
    
  Bool_t isStable = kFALSE;
  for (i = 0; i < kNstable; i++) {
    if (pdg == abs(pdgStable[i])) {
      isStable = kTRUE;
      break;
    }
  }
  
  return isStable;
}*/
/*
int is_stable_charged(const HepMC::GenParticle *part) {
  // copied from AliStack::IsStable()	 
  int pdg = abs(part->pdg_id());
  if(pdg>1000000000)return kTRUE;

  const Int_t kNstableCharged = 9;
  Int_t i;
  

  Int_t pdgStableCharged[kNstableCharged] = {
    kElectron,          // Electron
    kMuonPlus,          // Muon 
    kPiPlus,            // Pion
    kKPlus,             // Kaon
    kProton,            // Proton 
    kSigmaMinus,        // Sigma Minus
    kSigmaPlus,         // Sigma Plus
    3312,               // Xsi Minus 
    3334                // Omega
  };
    
  Bool_t isStable = kFALSE;
  for (i = 0; i < kNstableCharged; i++) {
    if (pdg == abs(pdgStableCharged[i])) {
      isStable = kTRUE;
      break;
    }
  }
  
  return isStable;
}
*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++
int isPi0(const HepMC::GenParticle *part) { 
    //Pi0 selection by pdg code
    int abs_kf = abs(part->pdg_id());
    if(abs_kf == 111) return 1;
    else return 0; 
} 
//++++++++++++++++++++++++++++++++++++++++++++++++++++
int is_charged(const HepMC::GenParticle *part) {
   // identify charged particle
   int abs_kf = abs(part->pdg_id());
   
   if (abs_kf==211 || abs_kf==321 || abs_kf==2212 || abs_kf==11 || abs_kf==13)
     return 1;
   else if (abs_kf != 22 && abs_kf!=111 && abs_kf!=130 && abs_kf!=2112 && abs_kf!=311 && abs_kf!=12 && abs_kf !=14 && abs_kf!=16)
     cout << " Unexpected particle: kf=" << abs_kf << endl;
   return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++
Float_t dphi(Float_t phi1, Float_t phi2){
   //calculate difference of azimuthal angles ph1 and phi2
   Float_t pi = 3.14159;

   if(phi1 < -pi){
      phi1+=2*pi;
   }else if(phi1 > pi){
      phi1-=2*pi;
   }

   if(phi2 < -pi){
      phi2+=2*pi;
   }else if(phi2 > pi){
      phi2-=2*pi;
   }


   Float_t dphii = phi1-phi2;

   if(dphii < -pi){
      dphii+=2*pi;

   }else if(dphii > pi){
      dphii-=2*pi;
   }

   return dphii;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++

Float_t phi02pi(Float_t phi){
   //converts azimuthal angle to 0,2pi  range
   Float_t pi = 3.14159;

   if(phi < 0){
      phi+=2*pi;

   }else if(phi > 2*pi){
      phi-=2*pi;
   }

   return phi;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++
void RecursiveParents(fastjet::PseudoJet jet,  THnSparse* hLund, Double_t dphi,  Double_t wevt);

//++++++++++++++++++++++++++++++++++++++++++++++++++++


int main(int argc, char **argv) {
  //
   // Takes 4 arguments: infile (HEPMC format), text Xsection file, low TT boarder,  high TT boarder 
   //

   TString name,title; 
 
   static const int debug = 0;
 
 
   if (argc < 4) {
     cerr << "Need 4 arguments: infile xsectionfile TTlow TThigh" << endl << "infile is HEPMC ascii format; outfile will be root format" << endl;
     return 1;
   }
   
   char *inname = argv[1];  //infile (HEPMC format)
   char *xnname = argv[2];  //xsectionfk.txt text file with pp,np, nn Xsections from log file 1_Xsection/extractXsection.sh
   // specify an input file
   HepMC::IO_GenEvent ascii_in(inname,std::ifstream::in);
   
   TT tt;
   tt.Set(atoi(argv[3]), atoi(argv[4])); // TT bin boarders

   TRandom3 *rnd = new TRandom3();
   rnd->SetSeed(0);

   name = Form("hjet_recoilOFF_output_TT%d_%d.root",tt.L(), tt.H()); //otuput file
   
   TFile fout(name.Data(),"RECREATE");
   static const int nR = 4; // jet cone radii 0.2, 0.3, 0.4, 0.5
   Float_t jetR[] = {0.2, 0.3, 0.4, 0.5};
   Float_t jetRKT= 0.4; //BG jet cone radius for KT jets
 
   //ACCEPTANCE OF DCAL AND EMCAL
   //  /home/fkrizek/alice/AliPhysics/PWGJE/EMCALJetTasks/UserTasks/AliAnalysisTaskEmcalDijetImbalance.cxx
   const Float_t eta_DCAL  = 0.7;
   const Float_t phi_DCAL1 = 260.*TMath::DegToRad();
   const Float_t phi_DCAL2 = 327.*TMath::DegToRad();

   const Float_t eta_EMCAL  = 0.7;
   const Float_t phi_EMCAL1 =  80.*TMath::DegToRad();
   const Float_t phi_EMCAL2 = 187.*TMath::DegToRad();

   //const Float_t eta_PHOS  = 0.13;
   //const Float_t phi_PHOS1 = 250.*TMath::DegToRad();
   //const Float_t phi_PHOS2 = 320.*TMath::DegToRad();


   //HISTOGRAMS 
   name ="hNEvent";
   TH1F *hNEvent = new TH1F(name.Data(),name.Data(),1,0,1); // event counter with jewel weight
   hNEvent->Sumw2();

   name = "Xsection";
   TProfile *hXsection = new TProfile(name.Data(),"X section [mb]",1,0,1,""); //jewel Xsection
   hXsection->Sumw2();

   //tracks
   name  = "hPtTrackEta"; //inclusive track pT versus eta
   title = "inclusive track p_{T} spectrum versus #eta; p_{T} (GeV/c); #eta";
   TH2F *hPtTrackEta = new TH2F(name.Data(),title.Data(),100,0,100,20,-1,1);
   hPtTrackEta->Sumw2();

   //pi0
   name  = "hEtaPhiEmcal"; //eta vesus phi for pi0 in emcal acceptance 
   title = "inclusive #phi vesus #eta in emcal; #phi; #eta";
   TH2F *hEtaPhiEmcal = new TH2F(name.Data(),title.Data(),50,0,2*TMath::Pi(),20,-1,1);
   hEtaPhiEmcal->Sumw2();

   name  = "hEtaPhiDcal"; //eta vesus phi for pi0 in dcal acceptance 
   title = "inclusive #phi vesus #eta in Dcal; #phi; #eta";
   TH2F *hEtaPhiDcal = new TH2F(name.Data(),title.Data(),50,0,2*TMath::Pi(),20,-1,1);
   hEtaPhiDcal->Sumw2();
 
   name  = "hPhiPi0"; // azimuthal angle distribution of inclusive pi0 
   title = "inclusive pi0 #phi; #phi; counts";
   TH1F *hPhiPi0 = new TH1F(name.Data(),title.Data(),100,0,2*TMath::Pi());
   hPhiPi0->Sumw2();
 
   name  = "hPtPi0Eta"; //inclusive eta versus pT distribution of pi0
   title = "inclusive gamma p_{T} spectrum vesus #eta; p_{T} (GeV/c); #eta";
   TH2F *hPtPi0Eta = new TH2F(name.Data(),title.Data(),100,0,100,20,-1,1);
   hPtPi0Eta->Sumw2();

   name  = "hPtPi0Eta_DCAL"; //inclusive eta versus pT distribution of pi0 in DCAL
   title = "inclusive gamma p_{T} spectrum vesus #eta in DCAL; p_{T} (GeV/c); #eta";
   TH2F *hPtPi0Eta_DCAL = new TH2F(name.Data(),title.Data(),100,0,100,20,-1,1);
   hPtPi0Eta_DCAL->Sumw2();

   name  = "hPtPi0Eta_EMCAL"; //inclusive eta versus pT distribution of pi0 in EMCAL
   title = "inclusive gamma p_{T} spectrum vesus #eta in EMCAL; p_{T} (GeV/c); #eta";
   TH2F *hPtPi0Eta_EMCAL = new TH2F(name.Data(),title.Data(),100,0,100,20,-1,1);
   hPtPi0Eta_EMCAL->Sumw2();


   //TRIGGERS
   name  = "hTTSinglIncl_EMCAL"; //pT spectrum of pi0 TT in EMCAL
   title = "trigger pi0 p_{T} spectrum (EMCAL pi0); p_{T} (GeV/c); cross section";
   TH1F* hTTSI_EMCAL = new TH1F(name.Data(),title.Data(),100,0,100); 
   hTTSI_EMCAL->Sumw2();
 
   name  = "hTTSinglIncl_DCAL"; //pT spectrum of pi0 TT in DCAL
   title = "trigger pi0 p_{T} spectrum (DCAL pi0); p_{T} (GeV/c); cross section";
   TH1F* hTTSI_DCAL = new TH1F(name.Data(),title.Data(),100,0,100); 
   hTTSI_DCAL->Sumw2();
 
   name  = "hTTSinglIncl_h"; //pT spectrum of hadron TT
   title = "trigger hadron p_{T} spectrum; p_{T} (GeV/c); cross section";
   TH1F* hTTSI_h = new TH1F(name.Data(),title.Data(),100,0,100); 
   hTTSI_h->Sumw2();
 
   name  = "hTTMult_EMCAL"; //mutiplicity of pi0 TT candidates in EMCAL 
   title = "trigger pi0 multiplicity EMCAL; multiplicity; cross section";
   TH1F* hTTMult_EMCAL = new TH1F(name.Data(),title.Data(),50,0,50); 
   hTTMult_EMCAL->Sumw2();

   name  = "hTTMult_DCAL";  //mutiplicity of pi0 TT candidates in DCAL 
   title = "trigger pi0 multiplicity DCAL; multiplicity; cross section";
   TH1F* hTTMult_DCAL = new TH1F(name.Data(),title.Data(),50,0,50); 
   hTTMult_DCAL->Sumw2();
 
   name  = "hTTMult_h";  //mutiplicity of hadron TT candidates 
   title = "trigger hadron multiplicity; multiplicity; cross section";
   TH1F* hTTMult_h = new TH1F(name.Data(),title.Data(),50,0,50); 
   hTTMult_h->Sumw2();
 
   //jets 
   TH2F *hPtJetEta_ch[nR];    //inclusive charged jet eta versus pT spectrum 
   TH2F *hPtJetEta_ch_Bgsub[nR];    //BG inclusive charged jet eta versus pT spectrum background corrected
   TH2F *hPtJetEta_full[nR];  //inclusive full jet eta versus pT spectrum in EMCAL
  
   for(Int_t ir = 0; ir < nR; ir++) {
      name  = Form("hPtJetEta_R%02d_ch", TMath::Nint(jetR[ir]*10));
      title = Form("charged jet spectrum R=%.1f;p_{T} (GeV/c);#eta", jetR[ir]);
      hPtJetEta_ch[ir] = new TH2F(name.Data(),title.Data(),500,0,500,20,-1,1);
      hPtJetEta_ch[ir]->Sumw2();
   }

   for(Int_t ir = 0; ir < nR; ir++) { //BG
      name  = Form("hPtJetEta_R%02d_ch_Bgsub", TMath::Nint(jetR[ir]*10));
      title = Form("charged jet spectrum R=%.1f bg corrected;p_{T} (GeV/c);#eta", jetR[ir]);
      hPtJetEta_ch_Bgsub[ir] = new TH2F(name.Data(),title.Data(),600,-100,500,20,-1,1);
      hPtJetEta_ch_Bgsub[ir]->Sumw2();
   }


   for(Int_t ir = 0; ir < nR; ir++) {
      name  = Form("hPtJetEta_R%02d_full", TMath::Nint(jetR[ir]*10));
      title = Form("full jet spectrum R=%.1f in EMCAL;p_{T} (GeV/c);#eta", jetR[ir]);
      hPtJetEta_full[ir] = new TH2F(name.Data(),title.Data(),500,0,500,20,-1,1);
      hPtJetEta_full[ir]->Sumw2();
   }

   //recoil jets
   TH2F *hRecoilJetTT_h_JetPt_ch_DPhi[nR];          //TT h            + charged jet   DeltaPhi versus pT
   TH2F *hRecoilJetTT_h_JetPt_ch_DPhi_Bgsub[nR];    //BG //TT h      + charged jet   DeltaPhi versus pT background subtracted
   TH2F *hRecoilJetTT_EMCAL_JetPt_ch_DPhi[nR];      //TT pi0 in EMCAL + charged jet 
   TH2F *hRecoilJetTT_EMCAL_JetPt_ch_DPhi_Bgsub[nR];//BG //TT pi0 in EMCAL + charged jet   (pT background subtracted)
   TH2F *hRecoilJetTT_DCAL_JetPt_full_DPhi[nR];     //TT pi0 in DCAL  + full jet in EMCAL

   for(Int_t ir = 0; ir < nR; ir++) {
      name  = Form("hRecoilJetTT_h_Pt_DPhi_ch_R%02d",TMath::Nint(jetR[ir]*10));
      title = Form("Charged recoil jet p_{T} vs dphi R=%.1f (TT hadron);p_{T,jet} (GeV/c);#Delta#varphi",jetR[ir]);
      hRecoilJetTT_h_JetPt_ch_DPhi[ir] = new TH2F(name.Data(), title.Data(), 500,0,500,48,0,fastjet::pi);
      hRecoilJetTT_h_JetPt_ch_DPhi[ir]->Sumw2();
   }
   for(Int_t ir = 0; ir < nR; ir++) { //BG
      name  = Form("hRecoilJetTT_h_Pt_DPhi_ch_R%02d_Bgsub",TMath::Nint(jetR[ir]*10));
      title = Form("Charged recoil jet p_{T} vs dphi R=%.1f (TT hadron);p_{T,jet}^{corr} (GeV/c);#Delta#varphi",jetR[ir]);
      hRecoilJetTT_h_JetPt_ch_DPhi_Bgsub[ir] = new TH2F(name.Data(), title.Data(), 550,-50,500,48,0,fastjet::pi);
      hRecoilJetTT_h_JetPt_ch_DPhi_Bgsub[ir]->Sumw2();
   }
   for(Int_t ir = 0; ir < nR; ir++) {
      name  = Form("hRecoilJetTT_EMCAL_Pt_DPhi_ch_R%02d",TMath::Nint(jetR[ir]*10));
      title = Form("Charged recoil jet p_{T} vs dphi R=%.1f (TT pi0 EMCAL);p_{T,jet} (GeV/c);#Delta#varphi",jetR[ir]);
      hRecoilJetTT_EMCAL_JetPt_ch_DPhi[ir] = new TH2F(name.Data(), title.Data(), 500,0,500,48,0,fastjet::pi);
      hRecoilJetTT_EMCAL_JetPt_ch_DPhi[ir]->Sumw2();
   } 
   for(Int_t ir = 0; ir < nR; ir++) { //BG
      name  = Form("hRecoilJetTT_EMCAL_Pt_DPhi_ch_R%02d_Bgsub",TMath::Nint(jetR[ir]*10));
      title = Form("Charged recoil jet p_{T} vs dphi R=%.1f (TT pi0 EMCAL);p_{T,jet}^{corr} (GeV/c);#Delta#varphi",jetR[ir]);
      hRecoilJetTT_EMCAL_JetPt_ch_DPhi_Bgsub[ir] = new TH2F(name.Data(), title.Data(), 550,-50,500,48,0,fastjet::pi);
      hRecoilJetTT_EMCAL_JetPt_ch_DPhi_Bgsub[ir]->Sumw2();
   } 
   for(Int_t ir = 0; ir < nR; ir++) {
      name  = Form("hRecoilJetTT_DCAL_Pt_DPhi_full_R%02d",TMath::Nint(jetR[ir]*10));
      title = Form("Full EMCAL recoil jet p_{T} vs dphi R=%.1f (TT pi0 DCAL);p_{T,jet} (GeV/c);#Delta#varphi",jetR[ir]);
      hRecoilJetTT_DCAL_JetPt_full_DPhi[ir] = new TH2F(name.Data(), title.Data(), 500,0,500,48,0,fastjet::pi);
      hRecoilJetTT_DCAL_JetPt_full_DPhi[ir]->Sumw2();
   }

       
   name  = "hRho";  //BG
   title = "background density without TT jet; #rho (GeV); counts";
   TH1F* hRho = new TH1F(name.Data(),title.Data(),200,0,200); 
   hRho->Sumw2();
 
   name  = "hRho2";  //BG
   title = "background density without pi0 jet; #rho (GeV); counts";
   TH1F* hRho2 = new TH1F(name.Data(),title.Data(),200,0,200); 
   hRho2->Sumw2();
 
   name  = "hRhoTT";  //BG
   title = "background density without TT jet (event with TT); #rho (GeV); counts";
   TH1F* hRhoTT = new TH1F(name.Data(),title.Data(),200,0,200); 
   hRhoTT->Sumw2();
 
   name  = "hRho2TT";  //BG
   title = "background density without pi0 jet (event with TT); #rho (GeV); counts";
   TH1F* hRho2TT = new TH1F(name.Data(),title.Data(),200,0,200); 
   hRho2TT->Sumw2();
 
   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //                           LUND PLOT
   const Int_t    dimSpec   = 3;  //inclusive jet
   const Int_t    nBinsSpec[3]  = {50, 50, 10};
   const Double_t lowBinSpec[3] = {0.0,-10, 0};
   const Double_t hiBinSpec[3]  = {5.0,  0,200};

   const Int_t    dimSpecTT   = 4; //recoil jet
   const Int_t    nBinsSpecTT[4]  = {50, 50, 10, 48};
   const Double_t lowBinSpecTT[4] = {0.0,-10, 0, 0};
   const Double_t hiBinSpecTT[4]  = {5.0,  0,200, fastjet::pi};

   THnSparse* fhLundIterative_ch[nR]; //inclusive charged jets
   THnSparse* fhLundIterative_fu[nR]; //inclusive full jets
   THnSparse* fhLundIterative_ch_h[nR]; //recoil charged jets recoiling from h TT
   THnSparse* fhLundIterative_ch_pi0EMCAL[nR]; //recoil charged jets recoiling from EMCAL pi0 TT
   THnSparse* fhLundIterative_fu_pi0DCAL[nR]; //recoil full jets from pi0 DCAl TT 

   for(Int_t ir = 0; ir < nR; ir++) {
      name  = Form("fHLundIterative_chJet_R%02d",TMath::Nint(jetR[ir]*10));
      title = Form("LundIterativePlot [log(1/theta),log(z*theta),pTjet] incl. charged jet R=%.1f",jetR[ir]);
      fhLundIterative_ch[ir] = new THnSparseF(name.Data(), title.Data(),
                                     dimSpec,nBinsSpec,lowBinSpec,hiBinSpec);
      fhLundIterative_ch[ir]->Sumw2(); 
   }
   //---------------
   for(Int_t ir = 0; ir < nR; ir++) {
      name  = Form("fHLundIterative_fuJet_R%02d",TMath::Nint(jetR[ir]*10));
      title = Form("LundIterativePlot [log(1/theta),log(z*theta),pTjet] incl. full jet R=%.1f",jetR[ir]);
      fhLundIterative_fu[ir] = new THnSparseF(name.Data(), title.Data(),
                                     dimSpec,nBinsSpec,lowBinSpec,hiBinSpec);
      fhLundIterative_fu[ir]->Sumw2();
   }
   //---------------
   for(Int_t ir = 0; ir < nR; ir++) {
      name  = Form("fHLundIterative_chJet_R%02d_h",TMath::Nint(jetR[ir]*10));
      title = Form("LundIterativePlot [log(1/theta),log(z*theta),pTjet] recoil charged jet R=%.1f (hadron TT)",jetR[ir]);
      fhLundIterative_ch_h[ir] = new THnSparseF(name.Data(), title.Data(),
                                     dimSpecTT,nBinsSpecTT,lowBinSpecTT,hiBinSpecTT);
      fhLundIterative_ch_h[ir]->Sumw2(); 
   }
   //---------------
   for(Int_t ir = 0; ir < nR; ir++) {
      name  = Form("fHLundIterative_chJet_R%02d_pi0EMCAL",TMath::Nint(jetR[ir]*10));
      title = Form("LundIterativePlot [log(1/theta),log(z*theta),pTjet] recoil charged jet R=%.1f (EMCAL pi0 TT)",jetR[ir]);
      fhLundIterative_ch_pi0EMCAL[ir] = new THnSparseF(name.Data(), title.Data(),
                                            dimSpecTT,nBinsSpecTT,lowBinSpecTT,hiBinSpecTT);
      fhLundIterative_ch_pi0EMCAL[ir]->Sumw2();
   }

   //---------------
   for(Int_t ir = 0; ir < nR; ir++) {
      name  = Form("fHLundIterative_fuJet_R%02d_pi0DCAL",TMath::Nint(jetR[ir]*10));
      title = Form("LundIterativePlot [log(1/theta),log(z*theta),pTjet] recoil full jet R=%.1f (DCAL pi0 TT)",jetR[ir]);
      fhLundIterative_fu_pi0DCAL[ir] = new THnSparseF(name.Data(), title.Data(),
                                           dimSpecTT,nBinsSpecTT,lowBinSpecTT,hiBinSpecTT);
      fhLundIterative_fu_pi0DCAL[ir]->Sumw2();
   }
   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //ACO//   CHARGED DIJET ACCOPLANARITY 
   const Int_t    dimAcpl   = 3;  // pTLead, pT sublead, Deltaphi
   const Int_t    nBinsAcpl[3]  = {50,   50, 30};
   const Double_t lowBinAcpl[3] = {0.0, 0.0,  0};
   const Double_t hiBinAcpl[3]  = {250.0, 250, fastjet::pi};


   THnSparse* fhDijetAccoplanarity_ch_ch[nR]; //recoil full jets from pi0 DCAl TT 
   THnSparse* fhDijetAccoplanarity_ch_ch_Bgsub[nR]; //BG recoil full jets from pi0 DCAl TT 

   for(Int_t ir = 0; ir < nR; ir++) {
      name  = Form("fhDijetAccoplanarity_ch_ch_R%02d",TMath::Nint(jetR[ir]*10));
      title = Form("Charged Dijet Accoplanarity [pTJetL,pTjetSL, dphi] for R=%.1f",jetR[ir]);
      fhDijetAccoplanarity_ch_ch[ir] = new THnSparseF(name.Data(), title.Data(),
                                     dimAcpl, nBinsAcpl, lowBinAcpl, hiBinAcpl);
      fhDijetAccoplanarity_ch_ch[ir]->Sumw2(); 
   }

   for(Int_t ir = 0; ir < nR; ir++) { //BG
      name  = Form("fhDijetAccoplanarity_ch_ch_R%02d_Bgsub",TMath::Nint(jetR[ir]*10));
      title = Form("Charged Dijet Accoplanarity [pTJetL,pTjetSL, dphi] for R=%.1f bg subtracted",jetR[ir]);
      fhDijetAccoplanarity_ch_ch_Bgsub[ir] = new THnSparseF(name.Data(), title.Data(),
                                     dimAcpl, nBinsAcpl, lowBinAcpl, hiBinAcpl);
      fhDijetAccoplanarity_ch_ch_Bgsub[ir]->Sumw2(); 
   }

 
   // https://arxiv.org/pdf/1609.05842.pdf  formula 3.3
   TH2F *hdr_chJet[nR]; //delta_r distribution of charged jet
   TH1F *hdr_sum[nR]; //sum of track pT inside delta_r  annulus

   for(Int_t ir = 0; ir < nR; ir++) {
      name  = Form("hRadialChTrkSubLeadJet_ch_R%02d",TMath::Nint(jetR[ir]*10));
      title = Form("Radial distrib of charged tracks around sub lead jet (p_{T,jet} > 30 GeV,  R=%.1f); #it{r}; #psi^{subleading}_{pT}",jetR[ir]);
      hdr_chJet[ir] = new TH2F(name.Data(), title.Data(), 35,0,0.7,2000,0,20);
      hdr_chJet[ir]->Sumw2();

      name  = Form("hSumPt_R%02d",TMath::Nint(jetR[ir]*10));
      title = Form("Dummy Sum of p_{T}; #it{r}; counts");
      hdr_sum[ir] = new TH1F(name.Data(), title.Data(), 35, 0, 0.7);
      hdr_sum[ir]->Sumw2();
   }


   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // read Xsection from the text file
   ifstream xfile(xnname,std::ifstream::in);

   if(!xfile){
      cout <<"ERROR: "<<xnname<<" doesn't exist !!!"<< endl;
      exit(-1);
   }

   TString tstr; 
   Double_t xsec;
   while(xfile.good()){
      tstr.ReadLine(xfile);
      tstr.Remove(0,tstr.First(":")+1);
      tstr.ReplaceAll("mb","");
      xsec = tstr.Atof();
      if( xsec > 0){
         hXsection->Fill(0.5,xsec); //average pp pn, nn cross sections
      } 
   }


 
   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // get the first event
   HepMC::GenEvent* evt = ascii_in.read_next_event();
   if(!evt)  cerr << "Input file not found " << inname << endl;
  
   Float_t pt_min_cut  = 0.150; // low cut on track pT 
   Float_t max_eta = 0.9;       // ALICE track acceptance cut
   Float_t phi_rnd;             // random angle which rotates the event plane
   TLorentzVector candidate, particle;
   std::vector<TLorentzVector>  candidateTT_EMCAL; // pi0 trigger in EMCAL candidates
   std::vector<TLorentzVector>  candidateTT_DCAL; // pi0 trigger in DCAL candidates
   std::vector<TLorentzVector>  candidateTT_h;    // charged hadron trigger candidates 
   std::vector <fastjet::PseudoJet> fjInputs_ch;  // charged jets 
   std::vector <fastjet::PseudoJet> fjInputs_full;// full jets in EMCAL 
   Double_t mom;
   Float_t phiii; // phi from 0 to 2pi
   Float_t leadJetPt, leadJetPhi, leadJetArea, dphi_jetjet, dphi_jettrk, deta_jettrk, dr_jettrk, dr_pt; // leading jet pT and phi    //ACO
   Float_t subleadJetPt, subleadJetPhi, subleadJetArea; // subleading jet pT and phi    //ACO
   Int_t leadidx, subleadidx, ibin; //ACO  index of the subleading jet
   Float_t dphicut = 5*fastjet::pi/6.0; //cut jet-jet dphi
   Double_t acoarr[3];  //ACO
   Int_t index_ch, index_full; //charged jet constituent index, full jet constituent index

   Double_t rhoKT    = 0.0;    //BG background density excluding jet with TT 
   Double_t rhoKTemc = 0.0;    //BG background density excluding jet that would contain EMCAL TT
   Bool_t   bJetWithTT = 0; //BG flag for jets containing TT
   Double_t frhovec[999];   //BG auxiliary array to store pT/A of kT jets
   Int_t    nJetAcc = 0;    //BG counter for jets  

   // loop until we run out of events
   while(evt){
      // analyze the event
      if(debug)  cout << "Event " << endl;
      candidateTT_EMCAL.clear();
      candidateTT_DCAL.clear();
      candidateTT_h.clear();
      fjInputs_ch.clear();
      fjInputs_full.clear();
      index_ch = 0;
      index_full = 0;

      //JEWEL DOES NOT RANDOMIZE EVENT PLANE  GAMMA IS PREFERENTIALLY EMITTED AT PI
      phi_rnd = TMath::TwoPi() * rnd->Uniform(0,1);

      //cout<<"+++++++++++++++++++++++++++++++++++++++++++"<<endl;
      hNEvent->Fill(0.5, evt->weights()[0]); // count events
    
      for(HepMC::GenEvent::particle_iterator pit = evt->particles_begin();
                pit != evt->particles_end(); ++pit){ //loop over particles

         const HepMC::GenParticle *p = *pit;
         if(!p->end_vertex() && p->status()==1 && fabs(p->momentum().eta()) < max_eta) { // final state charged particle
         //if (p->status()==1 && is_stable_charged(p) && fabs(p->momentum().eta()) < max_eta)  // charged primary (NB heavy flavour decay products probably not included)
            if(p->momentum().perp()<pt_min_cut) continue;

            phiii = phi02pi(p->momentum().phi() + phi_rnd); //randomize  phi  and convert it to 0--2pi range
            particle.SetPtEtaPhiE(p->momentum().perp(), p->momentum().eta(),  phiii, p->momentum().e());

            if(is_charged(p)){
               hPtTrackEta->Fill(particle.Pt(), particle.Eta(), evt->weights()[0]);
               
               
               
               mom = sqrt(particle.Px()*particle.Px() + 
                          particle.Py()*particle.Py() +
                          particle.Pz()*particle.Pz());
               
               //CHARGED JET
               fastjet::PseudoJet jInp(particle.Px(),particle.Py(),particle.Pz(),mom);
               jInp.set_user_index(index_ch);
               fjInputs_ch.push_back(jInp);
               index_ch++;

               // TT candidates charged h 
               if(tt.L()<particle.Pt() && particle.Pt()<tt.H()){    
                  candidate = particle;
                  candidateTT_h.push_back(candidate);
               }
            }//charged partiles


            //FULL JET IN EMCAL 
            if(TMath::Abs(particle.Eta())<eta_EMCAL && phi_EMCAL1 < phiii && phiii < phi_EMCAL2){


               mom = sqrt(particle.Px()*particle.Px() + 
                          particle.Py()*particle.Py() +
                          particle.Pz()*particle.Pz());
               
               
               fastjet::PseudoJet jInpf(particle.Px(),particle.Py(),particle.Pz(),mom);
               jInpf.set_user_index(index_full);
               fjInputs_full.push_back(jInpf);
               index_full++;
            }

            if(isPi0(p)){
               hPtPi0Eta->Fill(particle.Pt(), particle.Eta(), evt->weights()[0]);
               hPhiPi0->Fill(phiii, evt->weights()[0]);

               if(TMath::Abs(particle.Eta())<eta_EMCAL && phi_EMCAL1 < phiii && phiii <phi_EMCAL2){
                  hEtaPhiEmcal->Fill(phiii, particle.Eta()); //inclusive pi0 EMCAL spectrum
                  hPtPi0Eta_EMCAL->Fill(particle.Pt(), particle.Eta(), evt->weights()[0]);
               }

               if(TMath::Abs(particle.Eta())<eta_DCAL && phi_DCAL1 < phiii && phiii <phi_DCAL2){
                  hEtaPhiDcal->Fill(phiii, particle.Eta()); //inclusive pi0 DCAL spectrum
                  hPtPi0Eta_DCAL->Fill(particle.Pt(), particle.Eta(), evt->weights()[0]);
               }

               //pi0 TT  candidate 
               if(tt.L()<particle.Pt() && particle.Pt()<tt.H()){    
                  //PI0 IN EMCAL   
                  if(TMath::Abs(particle.Eta())<eta_EMCAL && phi_EMCAL1 < phiii && phiii <phi_EMCAL2){
                      candidate = particle;
                      candidateTT_EMCAL.push_back(candidate);
                  }
                  //PI0 IN DCAL   
                  if(TMath::Abs(particle.Eta())<eta_DCAL && phi_DCAL1 < phiii && phiii <phi_DCAL2){
                      candidate = particle;
                      candidateTT_DCAL.push_back(candidate);
                  }
               }//tt cut
            }//gamma
         } 
      } //end paticle loop  

 
      //_________________________________________________________________
      //SELECT TRIGGER PION  EMCAL
      Int_t itt_EMCAL=-1;
       
      hTTMult_EMCAL->Fill(candidateTT_EMCAL.size(), evt->weights()[0]);
 
      if(candidateTT_EMCAL.size() >0){ 
    
         //single inclusive
         itt_EMCAL = rnd->Integer(candidateTT_EMCAL.size());
         hTTSI_EMCAL->Fill( candidateTT_EMCAL[itt_EMCAL].Pt(), evt->weights()[0]);
      }

      //_________________________________________________________________
      //SELECT TRIGGER PION  DCAL
      Int_t itt_DCAL=-1;
       
      hTTMult_DCAL->Fill(candidateTT_DCAL.size(), evt->weights()[0]);
 
      if(candidateTT_DCAL.size() >0){ 
    
         //single inclusive
         itt_DCAL = rnd->Integer(candidateTT_DCAL.size());
         hTTSI_DCAL->Fill( candidateTT_DCAL[itt_DCAL].Pt(), evt->weights()[0]);
    
      }

      //_________________________________________________________________
      //SELECT TRIGGER hadron 
      Int_t itt_h=-1;
       
      hTTMult_h->Fill(candidateTT_h.size(), evt->weights()[0]);
 
      if(candidateTT_h.size() >0){ 
    
         //single inclusive
         itt_h = rnd->Integer(candidateTT_h.size());
         hTTSI_h->Fill( candidateTT_h[itt_h].Pt(), evt->weights()[0]);
      } 

      //_________________________________________________________________
 
 
///    hTTSinglIncl2050->Sumw2();
     // Do jet finding
     // Need R =0.2 and R=0.4 later on...
     fastjet::GhostedAreaSpec ghostSpec(max_eta,1,0.01);
     fastjet::Strategy               strategy = fastjet::Best;
     fastjet::RecombinationScheme    recombScheme = fastjet::BIpt_scheme;
     fastjet::AreaType areaType      = fastjet::active_area;
     fastjet::AreaDefinition areaDef = fastjet::AreaDefinition(areaType,ghostSpec);
 
     //fastjet::ClusterSequenceArea clustSeqChBkg(fjInputs_ch, *jetDefChBkg,areaDef);
     
 
 
     /*
       // Background jets -- not needed
     inclusiveJetsChBkg = clustSeqChBkg.inclusive_jets();
     sortedJetsChBkg    = sorted_by_pt(inclusiveJetsChBkg); 
     if(sortedJetsChBkg.size()>2) sortedJetsChBkg.erase(sortedJetsChBkg.begin(),sortedJetsChBkg.begin()+2);
 
     Double_t rho=0;
     Double_t sigma=0.;
     Double_t meanarea=0.;
    


 
     clustSeqChBkg.get_median_rho_and_sigma(sortedJetsChBkg, range, true, rho, sigma, meanarea, true);
     */
     Float_t jetEtaMax_ch; 
     Float_t jetEtaMax_full;
     Float_t jetPhi1_full;   
     Float_t jetPhi2_full;   
     Float_t jet_pt;
     Float_t jet_eta;
     Float_t dphi_jetTT;
     Float_t jet_area; //BG

     //---------------- Lund plot --------------
     fastjet::PseudoJet jj;
     fastjet::PseudoJet j1;
     fastjet::PseudoJet j2;
 
     //+++++++++++++++++++++++/+++++++++++++++++++++++/+++++++++++++++++++++++
     //BG  CALCULATE RHO KT JETS R=0.4  SKIP JET WHERE TT IS A CONSTITUENT 
     //+++++++++++++++++++++++/+++++++++++++++++++++++/+++++++++++++++++++++++
     fastjet::JetDefinition jetDefChKT(fastjet::kt_algorithm, jetRKT, recombScheme, strategy); //BG
     fastjet::ClusterSequenceArea clustSeqChKT(fjInputs_ch, jetDefChKT, areaDef); //BG
     vector <fastjet::PseudoJet> inclusiveJetsChKT = clustSeqChKT.inclusive_jets(); //BG
     Float_t jetEtaMax_ch_KT   = 0.9 - jetRKT; 

     rhoKT    = 0.0; 
     rhoKTemc = 0.0; 
     nJetAcc  = 0;

     for(unsigned int iJet = 0; iJet < inclusiveJetsChKT.size(); iJet++) {
        if(TMath::Abs(inclusiveJetsChKT[iJet].eta()) > jetEtaMax_ch_KT) continue;
        // ACCEPT GHOSTS  
        bJetWithTT = 0;   
        if(itt_h>-1){  //SKIP THE JET WITH TT CONSTITUENT
           std::vector<fastjet::PseudoJet> constituents = inclusiveJetsChKT[iJet].constituents();
           for(unsigned int iCst = 0; iCst < constituents.size(); iCst++) { //loop over constituents to seach for TT
              if(fabs(constituents[iCst].perp() - candidateTT_h[itt_h].Pt()) > 0.001) continue;
              if(fabs(constituents[iCst].eta()  - candidateTT_h[itt_h].Eta()) > 0.001) continue;
              if(fabs(dphi(constituents[iCst].phi(), (Float_t) candidateTT_h[itt_h].Phi())) > 0.001) continue;
                 
              bJetWithTT = 1; //TT found among the jet constituents
              break;   
           }
        }
        if(bJetWithTT){//skip the jet with TT
                 continue; 
        }
        frhovec[nJetAcc] = inclusiveJetsChKT[iJet].perp()/clustSeqChKT.area(inclusiveJetsChKT[iJet]);
        nJetAcc++;
     }
     
     if(nJetAcc){
        rhoKT = TMath::Median(nJetAcc, frhovec);
     } 
     hRho->Fill(rhoKTemc);
     //TT EMCAL
     nJetAcc  = 0;
     for(unsigned int iJet = 0; iJet < inclusiveJetsChKT.size(); iJet++) {
        bJetWithTT = 0;   
        if(TMath::Abs(inclusiveJetsChKT[iJet].eta()) > jetEtaMax_ch_KT) continue;
        // ACCEPT GHOSTS   if(inclusiveJetsChKT[iJet].perp() < pt_min_cut) continue;
        if(itt_EMCAL>-1){  //SKIP THE JET WITH TT
           deta_jettrk = inclusiveJetsChKT[iJet].eta() - candidateTT_EMCAL[itt_EMCAL].Eta();
           dphi_jetTT  = dphi(inclusiveJetsChKT[iJet].phi(), (Float_t) candidateTT_EMCAL[itt_EMCAL].Phi());
           if(sqrt(deta_jettrk*deta_jettrk + dphi_jetTT*dphi_jetTT) < jetRKT){ 
              bJetWithTT = 1; //TT found close to the jet axis 
           }
        }

        if(bJetWithTT){
           continue; 
        }
        frhovec[nJetAcc] = inclusiveJetsChKT[iJet].perp()/clustSeqChKT.area(inclusiveJetsChKT[iJet]);
        nJetAcc++;
     }
     
     if(nJetAcc){
        rhoKTemc = TMath::Median(nJetAcc, frhovec); 
     }
     hRho2->Fill(rhoKTemc);
     //BG end 

     //---------------------------------------------
     for(int ir = 0; ir < nR; ir++) {
        //fastjet::RangeDefinition range(-max_eta+jetR[ir], max_eta-jetR, 0, 2.*fastjet::pi); //REPLACE BY SOMETHING ELSE 

        fastjet::JetDefinition jetDefCh(fastjet::antikt_algorithm, jetR[ir], recombScheme, strategy);
        fastjet::JetDefinition jetDefFull(fastjet::antikt_algorithm, jetR[ir], recombScheme, strategy);
       
        fastjet::ClusterSequenceArea clustSeqCh(  fjInputs_ch, jetDefCh, areaDef); 
        fastjet::ClusterSequenceArea clustSeqFull(fjInputs_full, jetDefFull, areaDef);
       
        vector <fastjet::PseudoJet> inclusiveJetsCh   = clustSeqCh.inclusive_jets();
        vector <fastjet::PseudoJet> inclusiveJetsFull = clustSeqFull.inclusive_jets();
      
        //fiducial cuts on jet acceptance
        jetEtaMax_ch   = 0.9 - jetR[ir]; 
        jetEtaMax_full = 0.7 - jetR[ir]; 
        jetPhi1_full   = phi_EMCAL1 + jetR[ir]; 
        jetPhi2_full   = phi_EMCAL2 - jetR[ir]; 

        leadJetPt     = -1.0;//ACO
        leadJetPhi    = -1.0;//ACO
        leadJetArea   = -1.0;//ACO
        subleadJetPt  = -1.0;//ACO
        subleadJetPhi = -1.0;//ACO
        subleadJetArea= -1.0;//ACO
        leadidx       = -1; //ACO
        subleadidx    = -1; //ACO

        //+++++++++++++++++++++++
        //CHARGED  JETS
        for(unsigned int iJet = 0; iJet < inclusiveJetsCh.size(); iJet++) {
           if(TMath::Abs(inclusiveJetsCh[iJet].eta()) > jetEtaMax_ch) continue;
           if(inclusiveJetsCh[iJet].perp() < pt_min_cut) continue;
            
           jet_pt   = inclusiveJetsCh[iJet].perp();
           jet_eta  = inclusiveJetsCh[iJet].eta();
           jet_area = clustSeqCh.area(inclusiveJetsCh[iJet]); //BG
           hPtJetEta_ch[ir]->Fill(jet_pt, jet_eta, evt->weights()[0]); //inclusive jet spectrum
           hPtJetEta_ch_Bgsub[ir]->Fill(jet_pt -  jet_area*rhoKT, jet_eta, evt->weights()[0]); //inclusive jet spectrum
           RecursiveParents(inclusiveJetsCh[iJet],(THnSparse*)  fhLundIterative_ch[ir], -999., evt->weights()[0]); //LUND

	   //ACO
	   if(leadJetPt<=jet_pt){
              subleadJetPt   = leadJetPt;
              subleadJetPhi  = leadJetPhi;
              subleadJetArea = leadJetArea;//BG
              subleadidx     = leadidx; 

              leadJetPt      = jet_pt;
              leadJetPhi     = inclusiveJetsCh[iJet].phi();
              leadJetArea    = jet_area;//BG 
              leadidx        = iJet;
           }else if(subleadJetPt < jet_pt){
              subleadJetPt   = jet_pt;
              subleadJetPhi  = inclusiveJetsCh[iJet].phi();
              subleadJetArea = jet_area;//BG 
              subleadidx     = iJet;
           }
          
           //-------------------------- 
     
           //charged recoil jet spectrum associated with hadron TT
           if(itt_h>-1){
              dphi_jetTT = fabs(dphi(inclusiveJetsCh[iJet].phi(), (Float_t) candidateTT_h[itt_h].Phi()));
         
              hRecoilJetTT_h_JetPt_ch_DPhi[ir]->Fill(jet_pt, dphi_jetTT, evt->weights()[0]);

              hRecoilJetTT_h_JetPt_ch_DPhi_Bgsub[ir]->Fill(jet_pt - jet_area*rhoKT, dphi_jetTT, evt->weights()[0]);//BG
              hRhoTT->Fill(rhoKT);//BG

              RecursiveParents(inclusiveJetsCh[iJet], (THnSparse*) fhLundIterative_ch_h[ir], dphi_jetTT, evt->weights()[0]); //LUND
           }

           //charged recoil jet spectrum for pi0 EMCAL TT
           if(itt_EMCAL>-1){
              dphi_jetTT = fabs(dphi(inclusiveJetsCh[iJet].phi(), (Float_t) candidateTT_EMCAL[itt_EMCAL].Phi()));
         
              hRecoilJetTT_EMCAL_JetPt_ch_DPhi[ir]->Fill(jet_pt, dphi_jetTT, evt->weights()[0]);

              hRecoilJetTT_EMCAL_JetPt_ch_DPhi_Bgsub[ir]->Fill(jet_pt - jet_area*rhoKTemc, dphi_jetTT, evt->weights()[0]); //BG
              hRho2TT->Fill(rhoKTemc);//BG

              RecursiveParents(inclusiveJetsCh[iJet], (THnSparse*) fhLundIterative_ch_pi0EMCAL[ir], dphi_jetTT, evt->weights()[0]); //LUND
           }
        }
        //++++++++++++++
        //ACO
        if(subleadJetPt > 0.1 &&  leadJetPt > 0.1){
           dphi_jetjet = fabs(dphi(subleadJetPhi, leadJetPhi)); 
 
           acoarr[0] = leadJetPt;  
           acoarr[1] = subleadJetPt;  
           acoarr[2] = dphi_jetjet;  
           fhDijetAccoplanarity_ch_ch[ir]->Fill(acoarr); //recoil full jets from pi0 DCAl TT

           acoarr[0] = leadJetPt - rhoKT*leadJetArea; //BG 
           acoarr[1] = subleadJetPt - rhoKT*subleadJetArea;  //BG
           fhDijetAccoplanarity_ch_ch_Bgsub[ir]->Fill(acoarr); //recoil full jets from pi0 DCAl TT


           //radial distribution of tracks wrt subleading track  
           if(dphi_jetjet > dphicut){
              if(subleadidx > -1){
                 if(leadJetPt > 120 && subleadJetPt > 30){
                    if(TMath::Abs(inclusiveJetsCh[subleadidx].eta()) < 0.2){ // for jet in midrapidity 
                       hdr_sum[ir]->Reset(); 
                       // https://arxiv.org/pdf/1609.05842.pdf  formula 3.3
                       for(UInt_t ik=0; ik < fjInputs_ch.size(); ik++){ //loop over tracks
                             
                          dphi_jettrk = dphi( fjInputs_ch[ik].phi(), subleadJetPhi); 
                          deta_jettrk = fjInputs_ch[ik].eta() - inclusiveJetsCh[subleadidx].eta(); 
                          dr_jettrk = sqrt( dphi_jettrk*dphi_jettrk + deta_jettrk*deta_jettrk);//distance track-jet 
                          ibin = hdr_sum[ir]->FindBin(dr_jettrk);
                  
                          dr_pt = hdr_sum[ir]->GetBinContent(ibin) + fjInputs_ch[ik].perp();//sum pt in dr bin
                          hdr_sum[ir]->SetBinContent(ibin, dr_pt);
                       }
                    
                       for(Int_t ib=1; ib<=hdr_sum[ir]->GetNbinsX(); ib++){ //loop over dr bins 
                          hdr_chJet[ir]->Fill( hdr_sum[ir]->GetXaxis()->GetBinCenter(ib), hdr_sum[ir]->GetBinContent(ib)/subleadJetPt); 
                          //hdr_chJet[ir]->Fill( hdr_sum[ir]->GetXaxis()->GetBinCenter(ib), hdr_sum[ir]->GetBinContent(ib)/subleadJetPt/hdr_sum[ir]->GetBinWidth(ib)); 
                       }  
                    }           
                 }
              }
           }
        }
        //++++++++++++++
        //FULL JETS 
        for(unsigned int iJet = 0; iJet < inclusiveJetsFull.size(); iJet++) {
           if(TMath::Abs(inclusiveJetsFull[iJet].eta()) > jetEtaMax_full) continue;
           if(inclusiveJetsFull[iJet].perp() < pt_min_cut) continue;
           phiii =  phi02pi((Float_t) inclusiveJetsFull[iJet].phi());
           if(jetPhi1_full < phiii && phiii < jetPhi2_full){

              jet_pt  = inclusiveJetsFull[iJet].perp();
              jet_eta = inclusiveJetsFull[iJet].eta();
              hPtJetEta_full[ir]->Fill(jet_pt,jet_eta,evt->weights()[0]); //inclusive EMCAL full jet spectrum 
	      RecursiveParents(inclusiveJetsFull[iJet], (THnSparse*) fhLundIterative_fu[ir],-999., evt->weights()[0]); //LUND 
          
              //--------------------------------------------
              if(itt_DCAL>-1){ // PI0 DCAL TRIGGER
                 dphi_jetTT = fabs(dphi(inclusiveJetsFull[iJet].phi(), (Float_t) candidateTT_DCAL[itt_DCAL].Phi()));
              
                 hRecoilJetTT_DCAL_JetPt_full_DPhi[ir]->Fill(jet_pt, dphi_jetTT, evt->weights()[0]);
              
                 RecursiveParents(inclusiveJetsFull[iJet],(THnSparse*) fhLundIterative_fu_pi0DCAL[ir], dphi_jetTT, evt->weights()[0]); //LUND

              }
           }
        }
     }
    // delete the created event from memory
     delete evt;
     // read the next event
     ascii_in >> evt;    
   }



   for(int ir = 0; ir < nR; ir++) 
      fhLundIterative_ch[ir]->Write();
   for(int ir = 0; ir < nR; ir++) 
      fhLundIterative_fu[ir]->Write();
   for(int ir = 0; ir < nR; ir++)
      fhLundIterative_ch_h[ir]->Write();
   for(int ir = 0; ir < nR; ir++)
      fhLundIterative_ch_pi0EMCAL[ir]->Write();
   for(int ir = 0; ir < nR; ir++) 
      fhLundIterative_fu_pi0DCAL[ir]->Write();
   for(int ir = 0; ir < nR; ir++)  
      fhDijetAccoplanarity_ch_ch[ir]->Write(); //ACO 
   for(int ir = 0; ir < nR; ir++)  
      fhDijetAccoplanarity_ch_ch_Bgsub[ir]->Write(); //BG ACO 


   fout.Write();
   return 0;
} 



void RecursiveParents(fastjet::PseudoJet jet,  THnSparse* hLund, Double_t dphi,  Double_t wevt){
 //---------------- LUND PLOT --------------
        
   std::vector<fastjet::PseudoJet> constituents = jet.constituents();
   fastjet::JetAlgorithm jetalgo(fastjet::cambridge_algorithm);
   fastjet::JetDefinition fJetDef(jetalgo, 1., static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30 ); 
   try {

      fastjet::ClusterSequence fClustSeqSA(constituents, fJetDef);
      std::vector<fastjet::PseudoJet>   fOutputJets;
      fOutputJets.clear();
      fOutputJets = fClustSeqSA.inclusive_jets(0);
    
      fastjet::PseudoJet jj;
      fastjet::PseudoJet j1;
      fastjet::PseudoJet j2;

      jj = fOutputJets[0];
    
      while(jj.has_parents(j1,j2)){
         if(j1.perp() < j2.perp()) std::swap(j1,j2);
         double delta_R = j1.delta_R(j2);
         double z       = j2.perp()/(j1.perp()+j2.perp());
         double y       = log(1.0/delta_R);
         double lnpt_rel=log(z*delta_R);
 
         if(dphi<-99){ 
            Double_t LundEntries[3] = {y, lnpt_rel, fOutputJets[0].perp()};  
            hLund->Fill(LundEntries,  wevt);
         }else{
            Double_t LundEntries2[4] = {y, lnpt_rel, fOutputJets[0].perp(), dphi};  
            hLund->Fill(LundEntries2,  wevt);
         }  
         jj=j1;
      }
   
   } catch (fastjet::Error) {
      cout<<" [w] FJ Exception caught."<<endl;
      //return -1;
   }
   return;
}


  


     
