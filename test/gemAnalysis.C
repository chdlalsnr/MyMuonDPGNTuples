#define gemAnalysis_cxx
#include "gemAnalysis.h"
#include <TH2.h>
#include <TEfficiency.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <cmath>
#include <string>
#include <iostream>

gemAnalysis::gemAnalysis(const TString & inFileName,
                                 const TString & outFileName) :
  m_inFile(inFileName,"READ"),m_outFile(outFileName,"RECREATE"),fChain(0)
{

  fChain = static_cast<TTree*>(m_inFile.Get("muNtupleProducer/MuDPGTree"));
  Init(fChain);

}


void gemAnalysis::Loop()
{
//   In a ROOT session, you can do:
//      root> .L gemAnalysis.C
//      root> gemAnalysis t
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

   TCanvas *c = new TCanvas("c","canvas",800,800);

   TFile *f1 = new TFile("MuDPGNtuple_11_1_2_patch2.root","READ");
   std::cout << "File opened" << std::endl;

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

   TString var[] = {"gemDigi_g_phi","gemRecHit_g_phi"};

   TTree *tTree = (TTree*)f1->Get("muNtupleProducer/MuDPGTree");

   int n = sizeof(var)/sizeof(var[0]);
   int nBins[n];
   TH1F *histoPhi[n];
   TString binning;

   for(int i=0; i<n; i++){
     TString varname = var[i];
     TString s = std::to_string(i);

     binning = "(100,-3.15,3.15)";

     tTree->Draw(varname+">>histoPhi"+s+binning);
     
     histoPhi[i] = (TH1F*)gDirectory->Get("histoPhi"+s);

     //nBins[i] = histoPhi[i]->GetNbinsX();
     //std::cout << nBins << std::endl;
   }

   //TEfficiency *histoEff = new TEfficiency(&histoPhi[0],&histoPhi[1]);
   
   TH1F *histoEff = (TH1F*)histoPhi[1]->Clone("histoEff");
   //histoEff->SetMinimum(0.5);
   //histoEff->SetMaximum(1.4);
   histoEff->Sumw2();
   histoEff->Divide(histoPhi[0]);
   histoEff->SetMarkerStyle(21);
   histoEff->SetStats(0);
   histoEff->SetTitle("Efficiency recHits/Digis");
   histoEff->GetXaxis()->SetTitle("#phi");
   c->cd();
   histoEff->Draw("ep");
 
   m_plots["ScatteringDigisPlotXvsY"]->Draw("COLZ");
   m_plots["ScatteringRecHitsPlotXvsY"]->Draw("COLZ");


   TH2F *efficiencyMap = (TH2F*)m_plots["ScatteringRecHitsPlotXvsY"]->Clone("efficiencyMap");
   efficiencyMap->Sumw2();
   efficiencyMap->Divide(m_plots["ScatteringDigisPlotXvsY"]);
   efficiencyMap->SetMarkerStyle(21);
   efficiencyMap->SetStats(0);
   efficiencyMap->SetTitle("Efficiency Map");
   efficiencyMap->Draw("colz TEXT0");
    
   /*TCanvas *myC = new TCanvas("myC","ScatteringPlots",500,500);
   myC->cd();
   m_plots["ScatteringDigisPlotXvsY"]->Draw("COLZ");
   myC->SaveAs("./Plots/ScatteringDigisPlotXvsY.png");
   myC->Clear();
   m_plots["ScatteringRecHitsPlotXvsY"]->Draw("COLZ");
   myC->SaveAs("./Plots/ScatteringRecHitsPlotXvsY.png");
   myC->Clear();*/
   
   endJob();
}


void gemAnalysis::book()
{
  m_outFile.cd();

  m_plots["ScatteringDigisPlotXvsY"] = new TH2F("ScatteringDigisPlotXvsY",
                                                "Scatter Plot GEM Digis; globalX; globalY",
                                                100,-250.,250.,
                                                100,-250.,250.);

  m_plots["ScatteringRecHitsPlotXvsY"] = new TH2F("ScatteringRecHitsPlotXvsY",
                                                  "Scatter Plot GEM RecHits; globalX; global\
Y",
                                                  100,-250.,250.,
                                                  100,-250.,250.);
  
    
}

void gemAnalysis::fill()
{

  for(std::size_t iDigi = 0; iDigi < gemDigi_nDigis; ++iDigi)
    {
      Double_t x = gemDigi_g_x->at(iDigi);
      Double_t y = gemDigi_g_y->at(iDigi);

      m_plots["ScatteringDigisPlotXvsY"]->Fill(x,y);
      //m_plots["ScatteringDigisPlotXvsY"]->Draw("colz");
    }

  for(std::size_t iRecDigi = 0; iRecDigi < gemRecHit_nRecHits; ++ iRecDigi)
    {
      Double_t rec_x = gemRecHit_g_x->at(iRecDigi);
      Double_t rec_y = gemRecHit_g_y->at(iRecDigi);

      m_plots["ScatteringRecHitsPlotXvsY"]->Fill(rec_x,rec_y);
      //m_plots["ScatteringRecHitsPlotXvsY"]->Draw("colz");
    }
  
}

void gemAnalysis::endJob()
{

  m_outFile.cd();
  m_outFile.Write();
  m_outFile.Close();

}
