#include "TROOT.h"

void loadAnalysis()
{
  gROOT->ProcessLine(".L gemAnalysis.C++");
  //gROOT->ProcessLine("");
}
