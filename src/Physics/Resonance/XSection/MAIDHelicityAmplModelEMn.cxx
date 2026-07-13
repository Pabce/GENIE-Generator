//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#include "Physics/Resonance/XSection/MAIDHelicityAmplModelEMn.h"

using namespace genie;

//____________________________________________________________________________
MAIDHelicityAmplModelEMn::MAIDHelicityAmplModelEMn() :
MAIDHelicityAmplModelEMp("genie::MAIDHelicityAmplModelEMn", "Default")
{
  fUseNeutronAmplitudes = true;
}
//____________________________________________________________________________
MAIDHelicityAmplModelEMn::MAIDHelicityAmplModelEMn(string config) :
MAIDHelicityAmplModelEMp("genie::MAIDHelicityAmplModelEMn", config)
{
  fUseNeutronAmplitudes = true;
}
//____________________________________________________________________________
MAIDHelicityAmplModelEMn::~MAIDHelicityAmplModelEMn()
{

}
//____________________________________________________________________________
