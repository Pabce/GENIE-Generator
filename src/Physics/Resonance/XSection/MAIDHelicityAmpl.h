//____________________________________________________________________________
/*!
\class    genie::MAIDHelicityAmpl

\brief    Container for MAID electromagnetic resonance helicity amplitudes and
          the corresponding Rein-Sehgal-style f amplitudes.

\author   GENIE Collaboration

\created  June 25, 2026

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _MAID_HELICITY_AMPL_H_
#define _MAID_HELICITY_AMPL_H_

#include <iostream>

#include <TMath.h>

using std::ostream;

namespace genie {

class MAIDHelicityAmpl;
ostream & operator<< (ostream & stream, const MAIDHelicityAmpl & hamp);

class MAIDHelicityAmpl {

public:
  MAIDHelicityAmpl();
  MAIDHelicityAmpl(const MAIDHelicityAmpl & hamp);
  ~MAIDHelicityAmpl() { }

  void Reset(void);

  void SetAS(double a12, double a32, double s12);
  void SetRSF(double fm1, double fp1, double fm3, double fp3,
              double f0m, double f0p);
  void SetValid(bool valid);

  bool IsValid(void) const { return fValid; }

  double AmplA12(void) const { return fA12; }
  double AmplA32(void) const { return fA32; }
  double AmplS12(void) const { return fS12; }

  double Ampl2A12(void) const { return TMath::Power(fA12, 2.); }
  double Ampl2A32(void) const { return TMath::Power(fA32, 2.); }
  double Ampl2S12(void) const { return TMath::Power(fS12, 2.); }

  double AmpMinus1(void) const { return fMinus1; }
  double AmpPlus1 (void) const { return fPlus1;  }
  double AmpMinus3(void) const { return fMinus3; }
  double AmpPlus3 (void) const { return fPlus3;  }
  double Amp0Minus(void) const { return f0Minus; }
  double Amp0Plus (void) const { return f0Plus;  }

  double Amp2Minus1(void) const { return TMath::Power(fMinus1, 2.); }
  double Amp2Plus1 (void) const { return TMath::Power(fPlus1,  2.); }
  double Amp2Minus3(void) const { return TMath::Power(fMinus3, 2.); }
  double Amp2Plus3 (void) const { return TMath::Power(fPlus3,  2.); }
  double Amp20Minus(void) const { return TMath::Power(f0Minus, 2.); }
  double Amp20Plus (void) const { return TMath::Power(f0Plus,  2.); }

  friend ostream & operator<< (ostream & stream, const MAIDHelicityAmpl & hamp);

  void Print(ostream & stream) const;

private:
  bool   fValid;
  double fA12;
  double fA32;
  double fS12;
  double fMinus1;
  double fPlus1;
  double fMinus3;
  double fPlus3;
  double f0Minus;
  double f0Plus;
};

}        // genie namespace

#endif   // _MAID_HELICITY_AMPL_H_
