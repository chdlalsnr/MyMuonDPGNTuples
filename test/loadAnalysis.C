#include "TROOT.h"

void loadAnalysis()
{
  gROOT->ProcessLine(".L gemDigiAnalysis.C++");
  //gROOT->ProcessLine("");
}
