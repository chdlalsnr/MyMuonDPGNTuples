#define gemAnalysis_cxx
#include "gemAnalysis.h"
#include <TH2.h>
#include <TF1.h>
#include <TTree.h>
#include <TString.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

#include "vector"

gemAnalysis::gemAnalysis(const TString & inFileName,
                         const TString & outFileName) :
  m_inFile(inFileName,"READ"),m_outFile(outFileName,"RECREATE"),fChain(0)
{

  fChain = static_cast<TTree*>(m_inFile.Get("muNtupleProducer/MuDPGTree"));
  Init(fChain);

}

void gemAnalysis::Loop()
{

   if (fChain == 0) return;
   
   TCanvas *c = new TCanvas("c","canvas",800,800);

   TFile *f1 = new TFile("MuDPGNtuple_11_1_2_patch2.root","READ");
   std::cout << "File opened" << std::endl;

   Long64_t nentries = fChain->GetEntriesFast();
   std::size_t i=0;
   vector<float> recHit_Glb_x;
   recHit_Glb_x.clear();
   vector<float> propagated_Glb_x;
   propagated_Glb_x.clear();
   vector<float> matched_Glb_x_clsz_min2;
   matched_Glb_x_clsz_min2.clear();
   vector<float> residual_x_clsz_min2;
   residual_x_clsz_min2.clear();
   vector<float> matched_Glb_x_clsz_eq2;
   matched_Glb_x_clsz_eq2.clear();
   vector<float> residual_x_clsz_eq2;
   residual_x_clsz_eq2.clear();
   vector<float> matched_Glb_x_clsz_maj2;
   matched_Glb_x_clsz_maj2.clear();
   vector<float> residual_x_clsz_maj2;
   residual_x_clsz_maj2.clear();

   book();
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      /* if(gemRecHit_g_x->size() != 0 )
	 {
	   i=0;
	   while(i<gemRecHit_g_x->size()){
	     recHit_Glb_x.push_back(gemRecHit_g_x->at(i));
	     //std::cout << gemRecHit_g_x->at(i) << std::endl;
	     i++;
	   }
	   //std::cout << "rech" << std::endl;
	 }*/

       /*if(mu_propagatedGlb_x->size() != 0 )
         {
           i=0;
           while(i<mu_propagatedGlb_x->size()){
             propagated_Glb_x.push_back(mu_propagatedGlb_x->at(i));
	     //std::cout << mu_propagatedGlb_x->at(i) << std::endl;
             i++;
           }
	   //std::cout << "prop" << std::endl;
         }*/

      	 std::cout << "Run Number: " << event_runNumber << std::endl;
	 std::cout << "Event Number: " << event_eventNumber << std::endl;
	 /*if(mu_propagatedGlb_x->size() !=0 && gemRecHit_g_x->size() != 0)
	 {
	 	i=0;
		std::cout <<"sizeprop:"<<mu_propagatedGlb_x->size()<<std::endl;;
	 	while(i<mu_propagatedGlb_x->size())
	 	{
	 		std::cout <<"prop:"<< mu_propagatedGlb_x->at(i) << std::endl;
	 		i++;
	 	}
	 	i=0;
		std::cout <<"sizerec:"<<gemRecHit_g_x->size()<<std::endl;
	 	while(i<gemRecHit_g_x->size())
	 	{
	 		std::cout<< "rechit:"<< gemRecHit_g_x->at(i) << std::endl;
	 		i++;
	 	}
		}*/	
       
       if(mu_propagatedGlb_x->size() !=0 && gemRecHit_g_x->size() != 0)
	 {
	   if(gemRecHit_cluster_size->size() < 2 && gemRecHit_cluster_size->size() != 0)
	     {
	       vector<float> matched_x = findMatchedHit(gemRecHit_g_x,gemRecHit_g_x->size(),mu_propagatedGlb_x, mu_propagatedGlb_x->size());
	       i=0;
	       if(matched_x.size() != 0){
		 while(i< matched_x.size()){
		   matched_Glb_x_clsz_min2.push_back(matched_x.at(i));
		   i++;
		 }
	       }
	       vector<float> res_x = residual(gemRecHit_g_x,gemRecHit_g_x->size(),mu_propagatedGlb_x, mu_propagatedGlb_x->size());
	       i=0;
	       if(res_x.size() != 0){
		 while(i< res_x.size()){
		   residual_x_clsz_min2.push_back(res_x.at(i));
		   i++;
		 }
	       }
	     }
	   
	   else if(gemRecHit_cluster_size->size() == 2 && gemRecHit_cluster_size->size() != 0)
	     {
	       vector<float> matched_x = findMatchedHit(gemRecHit_g_x,gemRecHit_g_x->size(),mu_propagatedGlb_x, mu_propagatedGlb_x->size());
               i=0;
               if(matched_x.size() != 0){
                 while(i< matched_x.size()){
                   matched_Glb_x_clsz_eq2.push_back(matched_x.at(i));
                   i++;
                 }
               }
               vector<float> res_x = residual(gemRecHit_g_x,gemRecHit_g_x->size(),mu_propagatedGlb_x, mu_propagatedGlb_x->size());
               i=0;
               if(res_x.size() != 0){
                 while(i< res_x.size()){
                   residual_x_clsz_eq2.push_back(res_x.at(i));
                   i++;
                 }
               }
	     }
	   
	   else if(gemRecHit_cluster_size->size() > 2 && gemRecHit_cluster_size->size() != 0)
	     {
	       vector<float> matched_x = findMatchedHit(gemRecHit_g_x,gemRecHit_g_x->size(),mu_propagatedGlb_x, mu_propagatedGlb_x->size());
               i=0;
               if(matched_x.size() != 0){
                 while(i< matched_x.size()){
                   matched_Glb_x_clsz_maj2.push_back(matched_x.at(i));
                   i++;
                 }
               }
               vector<float> res_x = residual(gemRecHit_g_x,gemRecHit_g_x->size(),mu_propagatedGlb_x, mu_propagatedGlb_x->size());
               i=0;
               if(res_x.size() != 0){
                 while(i< res_x.size()){
                   residual_x_clsz_maj2.push_back(res_x.at(i));
                   i++;
                 }
               }
	     }
	 }	
   }

   TF1 *gauss = new TF1("gauss", "gaus", -0.5, .5);
   

   for(std::size_t iSize=0; iSize< residual_x_clsz_min2.size(); iSize++)
     {
       m_plots["ResidualPlotXCLSZMin2"]->Fill(residual_x_clsz_min2.at(iSize));
     }
   m_plots["ResidualPlotXCLSZMin2"]->Rebin(5);
   m_plots["ResidualPlotXCLSZMin2"]->Fit("gauss","R");

   for(std::size_t iSize=0; iSize< residual_x_clsz_eq2.size(); iSize++)
     {
       m_plots["ResidualPlotXCLSZEq2"]->Fill(residual_x_clsz_eq2.at(iSize));
     }
   m_plots["ResidualPlotXCLSZEq2"]->Rebin(5);
   m_plots["ResidualPlotXCLSZEq2"]->Fit("gauss","R");


   for(std::size_t iSize=0; iSize< residual_x_clsz_maj2.size(); iSize++)
     {
       m_plots["ResidualPlotXCLSZMaj2"]->Fill(residual_x_clsz_maj2.at(iSize));
     }
   m_plots["ResidualPlotXCLSZMaj2"]->Rebin(5);
   m_plots["ResidualPlotXCLSZMaj2"]->Fit("gauss","R");


   
   
   //fill_residual(residual_x_clsz_min2, residual_x_clsz_min2.size());
   //fill_matched(matched_Glb_x, matched_Glb_x.size());

      
   /*TString var[] = {"mu_propagatedGlb_x"};
   TTree *tTree = (TTree*)f1->Get("muNtupleProducer/MuDPGTree");

   int n = sizeof(var)/sizeof(var[0]);
   int nBins[n];
   int nEntries[n];
   TH1F *histo[n];
   TString binning;

   for(int i=0; i<n; i++){
     TString varname = var[i];
     TString s = std::to_string(i);

     if(varname=="mu_propagatedGlb_x") binning = "(100,-250.,250.)";
     
     tTree->Draw(varname+">>histo"+s+binning);

     histo[i] = (TH1F*)gDirectory->Get("histo"+s);

     nBins[i] = histo[i]->GetNbinsX();
     nEntries[i] = histo[i]->GetEntries();
     //std::cout << nBins << std::endl;
   }

   TH1F *histoMuGlbX = (TH1F*)histo[0]->Clone("OccupancyPlotMuGlbX");
   histoMuGlbX->SetTitle("Occupancy Plot Propagated X");
   histoMuGlbX->GetXaxis()->SetTitle("global_x");
   c->cd();
   histoMuGlbX->Draw();
   
   TH1F *histoEff = (TH1F*)m_plots["OccupancyMatchedPlotX"]->Clone("histoEff");
   //histoEff->SetMinimum(0.5);
   //histoEff->SetMaximum(1.4);
   histoEff->Sumw2();
   histoEff->Divide(histo[0]);
   histoEff->SetMarkerStyle(21);
   histoEff->SetStats(0);
   histoEff->SetTitle("Efficiency matched/propagated");
   histoEff->GetXaxis()->SetTitle("global_x");
   c->cd();
   histoEff->Draw("ep");*/

   endJob();

}

vector<float> gemAnalysis::findMatchedHit(vector<float> *recHitPositions, size_t nRecHitPos, vector<float> *muPropagatedPositions, size_t nMuPropagatedPos)
{

  float min_recHitPosition = 666; //DUMMY VALUE
  vector<float> min_recHitPositions;
  //recHitMatched.clear();
  
  for(std::size_t iMu = 0; iMu < nMuPropagatedPos; iMu++)
    {
      float min_residual = 100;
      for(std::size_t iRecHit = 0; iRecHit < nRecHitPos; iRecHit++)
	{
	  float residual = std::fabs(muPropagatedPositions->at(iMu)-recHitPositions->at(iRecHit));
	  if(residual <= min_residual)
	    {
	      min_residual = residual;
	      min_recHitPosition = recHitPositions->at(iRecHit);
	    }
	  else continue;
	}
      if(min_recHitPosition == 666) continue;
      else
	{
	  min_recHitPositions.push_back(min_recHitPosition);
	  //residuals.push_back(min_residual);
	}  
    }

  return min_recHitPositions;

}

vector<float> gemAnalysis::residual(vector<float> *recHitPositions, size_t nRecHitPos, vector<float> *muPropagatedPositions, size_t nMuPropagatedPos)
{

  float min_recHitPosition = 666; //DUMMY VALUE
  vector<float> residuals;
  //recHitMatched.clear();

  for(std::size_t iMu = 0; iMu < nMuPropagatedPos; iMu++)
    {
      float min_residual = 100;
      for(std::size_t iRecHit = 0; iRecHit < nRecHitPos; iRecHit++)
        {
          float residual = (muPropagatedPositions->at(iMu)-recHitPositions->at(iRecHit));
          if(std::fabs(residual) <= std::fabs(min_residual))
            {
              min_residual = residual;
              min_recHitPosition = recHitPositions->at(iRecHit);
            }
          else continue;
        }
      if(min_recHitPosition == 666) continue;
      else
        {
          residuals.push_back(min_residual);
        }
    }

  return residuals;

}


void gemAnalysis::book()
{
  m_outFile.cd();

  m_plots["OccupancyMatchedPlotX"] = new TH1F("OccupancyMatchedPlotX",
					      "Occupancy Matched X; globalX ; entries",
					      100,-250.,250.);

  m_plots["ResidualPlotXCLSZMin2"] = new TH1F("ResidualPlotXCLSZMin2",
                                              "Residual X ClusterSize<2; residualX ; entries",
                                              400,-2,2.);
  m_plots["ResidualPlotXCLSZEq2"] = new TH1F("ResidualPlotXCLSZEq2",
                                              "Residual X ClusterSize=2; residualX ; entries",
                                              400,-2,2.);

  m_plots["ResidualPlotXCLSZMaj2"] = new TH1F("ResidualPlotXCLSZMaj2",
                                              "Residual X ClusterSize>2; residualX ; entries",
                                              400,-2,2.);

  
}

void gemAnalysis::fill_matched(vector<float> position, size_t nSize)
{
  
  for(std::size_t size = 0; size < nSize; size++)
    {
      m_plots["OccupancyMatchedPlotX"]->Fill(position.at(size));
    }
    
}

void gemAnalysis::fill_residual(vector<float> position, size_t nSize)
{

  for(std::size_t size = 0; size < nSize; size++)
    {
      m_plots["ResidualPlotXCLSZMin2"]->Fill(position.at(size));
    }

  TF1 f1("f1","gaus",-2.,2.);
  f1.SetParameters(1,m_plots["ResidualPlotXCLSZMin2"]->GetMean(),m_plots["ResidualPlotXCLSZMin2"]->GetRMS());
  //f1.SetParameter(1,0.03);
  f1.SetParLimits(1,-.5,.5);
  f1.SetParameter(2,0.9);
  //f1.SetParLimits(2,.5,1.);
  m_plots["ResidualPlotXCLSZMin2"]->Rebin(5);
  m_plots["ResidualPlotXCLSZMin2"]->Fit("f1","R");


}

						   
void gemAnalysis::endJob()
{

  m_outFile.cd();
  m_outFile.Write();
  m_outFile.Close();

}

