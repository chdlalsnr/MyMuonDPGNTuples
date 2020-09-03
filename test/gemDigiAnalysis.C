#define gemDigiAnalysis_cxx
#include "gemDigiAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h> 
#include <iostream>

gemDigiAnalysis::gemDigiAnalysis(const TString & inFileName,
				 const TString & outFileName) :
  m_inFile(inFileName,"READ"),m_outFile(outFileName,"RECREATE"),fChain(0)
{

  fChain = static_cast<TTree*>(m_inFile.Get("muNtupleProducer/MuDPGTree"));
  Init(fChain);

}


void gemDigiAnalysis::Loop()
{
//   In a ROOT session, you can do:
//      root> .L gemDigiAnalysis.C
//      root> gemDigiAnalysis t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  
   book();
   
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      //std::cout<<"ok"<<std::endl;
      if (ientry < 0)
	break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      fill();
      // if (Cut(ientry) < 0) continue;
   }
   
   endJob();
   
}

void gemDigiAnalysis::book()
{

  m_outFile.cd();

  m_plots["ScatterPlotXvsY"] = new TH2F("ScatterPlotXvsY",
					"Scatter Plot; X; Y",
					100,-250.,250.,
					100,-250.,250.);

}

void gemDigiAnalysis::fill()
{

  //std::cout<<gemDigi_nDigis<<std::endl;
  for(std::size_t iDigi = 0; iDigi<gemDigi_nDigis ; ++iDigi)
    {
      Double_t x = gemDigi_g_x->at(iDigi);
      Double_t y = gemDigi_g_y->at(iDigi);
      //std::cout << "x= " << x << std::endl; 
      m_plots["ScatterPlotXvsY"]->Fill(x, y);
    }
}


void gemDigiAnalysis::endJob()
{

  m_outFile.cd();

  m_outFile.Write();
  m_outFile.Close();

}
