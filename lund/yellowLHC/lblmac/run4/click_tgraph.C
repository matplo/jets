TGraph *g;
Int_t i = 0;

void myexec()
{
   if (!gPad) return;

   // get event information
   int event = gPad->GetEvent();
   int px    = gPad->GetEventX();
   int py    = gPad->GetEventY();

   // some magic to get the coordinates...
   double xd = gPad->AbsPixeltoX(px);
   double yd = gPad->AbsPixeltoY(py);
   float x = gPad->PadtoX(xd);
   float y = gPad->PadtoY(yd);

   if (event==1) { // left mouse button click
      g->SetPoint(i,x,y);
      if (i==0) g->Draw("L");
      i++;
      gPad->Update();
      return;
   }

   if (event != 51 && event != 52 && event != 53)
   {
      cout <<  event << endl;
   }

   if (event==24) { // keyboard pressed
      TFile *fout = new TFile("click_tgraph.root", "recreate");
      g->Write();
      cout << "dumpding graph... to " << fout->GetName() << endl;
      fout->Close();
      delete fout;
   }
}


void click_tgraph()
{
   g = new TGraph();
   // add exec
   gPad->AddExec("myexec","myexec()");
}
