void myexec()
{
   // get event information
   int event = gPad->GetEvent();
   int px    = gPad->GetEventX();
   int py    = gPad->GetEventY();

   // some magic to get the coordinates...
   double xd = gPad->AbsPixeltoX(px);
   double yd = gPad->AbsPixeltoY(py);
   float x = gPad->PadtoX(xd);
   float y = gPad->PadtoY(yd);

   cout << "x=" << x << " y=" << y << endl;

   if (event==1) { // left mouse button click
      cout << "x=" << x << " y=" << y << endl;
   }
}
