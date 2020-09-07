//#define gemDigiAnalysis_cxx
#include "gemDigiAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
//#include <iostream.h>

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
   if (fChain == 0) return;

   book();
   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      fill();
   }
}


void gemDigiAnalysis::book()
{
  m_outFile.cd();

  m_plots["ScatteringDigisPlotXvsY"] = new TH2F("ScatteringDigisPlotXvsY",
						"Scatter Plot GEM Digis; globalX; globalY",
						100,-250.,250.,
						100,-250.,250.);

  m_plots["ScatteringRecHitsPlotXvsY"] = new TH2F("ScatteringRecHitsPlotXvsY",
						  "Scatter Plot GEM RecHits; globalX; globalY",
						  100,-250.,250.,
						  100,-250.,250.);



}

void gemDigiAnalysis::fill()
{

  for(std::size_t iDigi = 0; iDigi < gemDigi_nDigis; ++iDigi)
    {  
      Double_t x = gemDigi_g_x->at(iDigi);
      Double_t y = gemDigi_g_y->at(iDigi);
    
      m_plots["ScatteringDigisPlotXvsY"]->Fill(x,y);
    }

  for(std::size_t iRecDigi = 0; iRecDigi < gemRecHit_nRecHits; ++ iRecDigi)
    {
      Double_t rec_x = gemRecHit_g_x->at(iRecDigi);
      Double_t rec_y = gemRecHit_g_y->at(iRecDigi);

      m_plots["ScatteringRecHitsPlotXvsY"]->Fill(rec_x,rec_y);
    }
}
    
void gemDigiAnalysis::endJob()
{

  m_outFile.cd();
  m_outFile.Write();
  m_outFile.Close();

}

