#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TFile.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMath.h"
#include <stdio.h>
#include <stdlib.h>
#include "TColor.h"
#include <iostream>
#include "TGraphErrors.h"


#include <iostream>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;



Int_t LundPlot_PbPb(){

   Double_t lumi_mb = 10000000.;//10nb-1;   integrated luminosity
   Double_t TPbPb = 23.42; //the mean nuclear thickness function for 0-10% centrality bin of PbPb at sqrt{s}NN=5.02 TeV
   Double_t fcent = 0.1;  //fraction of the total hadronic Xsection simulated by JEWEL  (0-10% centrality bin)
   Double_t XsectionPbPb = 7700; //mb total Xsection of pbpb

   //TRIGGER HADRON RANGE
   Int_t ttl=20; //TT trigger pT low
   Int_t tth=25; //TT trigger pT high


   //HARD BINS
   Int_t hardbins[]={5,10,15,20,30,40,60,80,110,140,180,220,270,320,400};
   const Int_t khb = 14;  //the number of hard bins 

   //+++++++++++++++++++++++++++++++++++++++++
   // JEWEL RECOIL MODE 
   
   TString indir;   //MATEUSZ, HERE CHANGE PATH TO DATA
   indir = Form("/home/fkrizek/ANALYSIS/JEWEL/anal4b/pbpb_TT%d_%d_h_norecoil",ttl,tth);//WITHOUT RECOIL

   //+++++++++++++++++++++++++++++++++++++++++
   TFile *f;
   Int_t b1,b2,B1,B2;
   Double_t norm;
   TString name;

   TH1F* hE;  //event counting histogram in given hard bin
   TProfile *hXsec;  //Xsection in given hard bin
   Double_t nevt, xsec, scale;

   Float_t jetR[]={0.2,0.3,0.4,0.5}; //jet radii
   const Int_t nR = 4;               //number of jet radii


   //LundIterativePlots[log(1/theta),log(z*theta),pTjet] incl. charged jet R=
   THnSparse* fHLundIterative_chJet_R_Sum[nR];// LUND plot inclusive CHARGED jets with given R (sum over hard bins)
   THnSparse* fHLundIterative_chJet_R[nR];    // LUND plot inclusive CHARGED jets with given R (in given hard bin)

   THnSparse* fHLundIterative_fuJet_R_Sum[nR];// LUND plot inclusive FULL jets with given R (sum over hard bins)
   THnSparse* fHLundIterative_fuJet_R[nR];    // LUND plot inclusive FULL jets with given R (in given hard bin)

   THnSparse* fHLundIterative_chJet_R_h_Sum[nR]; // LUND plot CHARGED jets recoiling from TT charged hadron (sum over hard bins)
   THnSparse* fHLundIterative_chJet_R_h[nR]; // LUND plot CHARGED jets recoiling from TT charged hadron (in given hard bin)

   THnSparse* fHLundIterative_chJet_R_pi0EMCAL_Sum[nR]; // LUND plot CHARGED jets recoiling from TT pi0 EMCAL (sum over hard bins)
   THnSparse* fHLundIterative_chJet_R_pi0EMCAL[nR]; // LUND plot CHARGED jets recoiling from TT pi0 EMCAL (in given hard bin)

   THnSparse* fHLundIterative_fuJet_R_pi0DAL_Sum[nR]; // LUND plot FULL jets recoiling from TT pi0 DCAL (sum over hard bins)
   THnSparse* fHLundIterative_fuJet_R_pi0DCAL[nR]; // LUND plot FULL jets recoiling from TT pi0 DCAL (in given hard bin)


   //LOOP OVER HARD BINS, READ HISTOS, NORMALIZE THEM  AND SUM THEM
   for(Int_t ih=0; ih<khb; ih++){    

      //input file for given hard bin and recoil mode
      name = Form("hb%d_%d_TT%d_%d.root", hardbins[ih],hardbins[ih+1],ttl,tth);
          
      f = new TFile(Form("%s/%s",indir.Data(), name.Data()),"READ"); //open file for given hard bin
      f->cd();
     
      hE = (TH1F*) f->Get("hNEvent");  //the number of events in given hard bin
      if(!hE) return -1;
      nevt = hE->GetBinContent(1);
     
      hXsec = (TProfile*) f->Get("Xsection"); // cross section of given hard bin
      if(!hXsec) return -2;
      xsec = hXsec->GetBinContent(1);
     
      //XsectionPbPb*fcent*lumi_mb =number of central collisions;   
      //TPbPb*xsec = number of inelastic collisions 
      scale = (XsectionPbPb*fcent*lumi_mb) * (TPbPb*xsec)/nevt;  // normalization coefficient for given hard bin
    
      //----------------
      for(Int_t ir=0; ir<nR; ir++){ //loop over jet radii
         //LUND PLOT FOR INCLUSIVE CHARGED JETS
         name = Form("fHLundIterative_chJet_R%02d",TMath::Nint(10*jetR[ir]));  
         fHLundIterative_chJet_R[ir] = (THnSparse*) f->Get(name.Data()); 
         if(!fHLundIterative_chJet_R[ir]) return -7;
         fHLundIterative_chJet_R[ir]->Scale(scale);

         //LUND PLOT FOR EMCAL INCLUSIVE FULL JETS
         name = Form("fHLundIterative_fuJet_R%02d",TMath::Nint(10*jetR[ir]));  
         fHLundIterative_fuJet_R[ir] = (THnSparse*) f->Get(name.Data()); 
         if(!fHLundIterative_fuJet_R[ir]) return -7;
         fHLundIterative_fuJet_R[ir]->Scale(scale);
   
         //LUND PLOT FOR CHARGED RECOIL JETS TT IS CHARGED HADRON
         name = Form("fHLundIterative_chJet_R%02d_h",TMath::Nint(10*jetR[ir]));  
         fHLundIterative_chJet_R_h[ir] = (THnSparse*) f->Get(name.Data()); 
         if(!fHLundIterative_chJet_R_h[ir]) return -7;
         fHLundIterative_chJet_R_h[ir]->Scale(scale);

         //LUND PLOT FOR CHARGED RECOIL JETS TT IS EMCAL PI0
         name = Form("fHLundIterative_chJet_R%02d_pi0EMCAL",TMath::Nint(10*jetR[ir]));  
         fHLundIterative_chJet_R_pi0EMCAL[ir] = (THnSparse*) f->Get(name.Data()); 
         if(!fHLundIterative_chJet_R_pi0EMCAL[ir]) return -7;
         fHLundIterative_chJet_R_pi0EMCAL[ir]->Scale(scale);

         //LUND PLOT FOR FULL RECOIL JETS TT IS DCAL PI0
         name = Form("fHLundIterative_fuJet_R%02d_pi0DCAL",TMath::Nint(10*jetR[ir]));  
         fHLundIterative_fuJet_R_pi0DCAL[ir] = (THnSparse*) f->Get(name.Data()); 
         if(!fHLundIterative_fuJet_R_pi0DCAL[ir]) return -7;
         fHLundIterative_fuJet_R_pi0DCAL[ir]->Scale(scale);
      } 


      //SUM HARD BINS
      if(ih==0){   //MAKE A COPY OF THE HISTOS IN THE FIRST HARD BIN
         for(Int_t ir=0; ir<nR; ir++){ 
            name = Form("fHLundIterative_chJet_R%02d_Sum", TMath::Nint(jetR[ir]*10));
            fHLundIterative_chJet_R_Sum[ir] = (THnSparse*) fHLundIterative_chJet_R[ir]->Clone(name.Data());
 
            name = Form("fHLundIterative_fuJet_R%02d_Sum", TMath::Nint(jetR[ir]*10));
            fHLundIterative_fuJet_R_Sum[ir] = (THnSparse*) fHLundIterative_fuJet_R[ir]->Clone(name.Data());
            //----------
            name = Form("fHLundIterative_chJet_R%02d_h_Sum", TMath::Nint(jetR[ir]*10));
            fHLundIterative_chJet_R_h_Sum[ir] = (THnSparse*) fHLundIterative_chJet_R_h[ir]->Clone(name.Data());

            name = Form("fHLundIterative_chJet_R%02d_pi0EMCAL_Sum", TMath::Nint(jetR[ir]*10));
            fHLundIterative_chJet_R_pi0EMCAL_Sum[ir] = (THnSparse*) fHLundIterative_chJet_R_pi0EMCAL[ir]->Clone(name.Data());

            name = Form("fHLundIterative_fuJet_R%02d_pi0DAL_Sum", TMath::Nint(jetR[ir]*10));
            fHLundIterative_fuJet_R_pi0DAL_Sum[ir] = (THnSparse*) fHLundIterative_fuJet_R_pi0DCAL[ir]->Clone(name.Data());
         }
      }else{ //ADD CORRESPONDING HISTOGRAMS FROM OTHER HARD BINS
         for(Int_t ir=0; ir<nR; ir++){ 
            fHLundIterative_chJet_R_Sum[ir]->Add((THnSparse*) fHLundIterative_chJet_R[ir]);
            fHLundIterative_fuJet_R_Sum[ir]->Add((THnSparse*) fHLundIterative_fuJet_R[ir]);
            fHLundIterative_chJet_R_h_Sum[ir]->Add((THnSparse*) fHLundIterative_chJet_R_h[ir]);
            fHLundIterative_chJet_R_pi0EMCAL_Sum[ir]->Add((THnSparse*) fHLundIterative_chJet_R_pi0EMCAL[ir]);
            fHLundIterative_fuJet_R_pi0DAL_Sum[ir]->Add((THnSparse*) fHLundIterative_fuJet_R_pi0DCAL[ir]);
         }
     
         f->Close();
      }
   } 

   //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //MAKE PROJECTIONS AND DRAWING ... 


   return 1;
}



