//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#include "Physics/Resonance/XSection/MAIDHelicityAmpl.h"

using namespace genie;
using std::endl;

//____________________________________________________________________________
namespace genie {
  ostream & operator<< (ostream & stream, const MAIDHelicityAmpl & hamp)
  {
     hamp.Print(stream);
     return stream;
  }
}
//____________________________________________________________________________
MAIDHelicityAmpl::MAIDHelicityAmpl()
{
  this->Reset();
}
//____________________________________________________________________________
MAIDHelicityAmpl::MAIDHelicityAmpl(const MAIDHelicityAmpl & hamp)
{
  fValid  = hamp.IsValid();
  fA12    = hamp.AmplA12();
  fA32    = hamp.AmplA32();
  fS12    = hamp.AmplS12();
  fMinus1 = hamp.AmpMinus1();
  fPlus1  = hamp.AmpPlus1();
  fMinus3 = hamp.AmpMinus3();
  fPlus3  = hamp.AmpPlus3();
  f0Minus = hamp.Amp0Minus();
  f0Plus  = hamp.Amp0Plus();
}
//____________________________________________________________________________
void MAIDHelicityAmpl::Reset(void)
{
  fValid  = false;
  fA12    = 0.0;
  fA32    = 0.0;
  fS12    = 0.0;
  fMinus1 = 0.0;
  fPlus1  = 0.0;
  fMinus3 = 0.0;
  fPlus3  = 0.0;
  f0Minus = 0.0;
  f0Plus  = 0.0;
}
//____________________________________________________________________________
void MAIDHelicityAmpl::SetAS(double a12, double a32, double s12)
{
  fValid = true;
  fA12 = a12;
  fA32 = a32;
  fS12 = s12;
}
//____________________________________________________________________________
void MAIDHelicityAmpl::SetRSF(
  double fm1, double fp1, double fm3, double fp3, double f0m, double f0p)
{
  fValid = true;
  fMinus1 = fm1;
  fPlus1  = fp1;
  fMinus3 = fm3;
  fPlus3  = fp3;
  f0Minus = f0m;
  f0Plus  = f0p;
}
//____________________________________________________________________________
void MAIDHelicityAmpl::SetValid(bool valid)
{
  fValid = valid;
}
//____________________________________________________________________________
void MAIDHelicityAmpl::Print(ostream & stream) const
{
  stream << endl;
  stream << " valid  = " << fValid  << endl;
  stream << " A(1/2) = " << fA12    << endl;
  stream << " A(3/2) = " << fA32    << endl;
  stream << " S(1/2) = " << fS12    << endl;
  stream << " f(-1)  = " << fMinus1 << endl;
  stream << " f(+1)  = " << fPlus1  << endl;
  stream << " f(-3)  = " << fMinus3 << endl;
  stream << " f(+3)  = " << fPlus3  << endl;
  stream << " f(0-)  = " << f0Minus << endl;
  stream << " f(0+)  = " << f0Plus  << endl;
}
//____________________________________________________________________________
