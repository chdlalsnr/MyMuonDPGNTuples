#define gemAnalysis_cxx
#include "gemAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TTree.h>
#include <TString.h>
#include <iostream>
#include <string.h>

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

   Long64_t nentries = fChain->GetEntries();
   
   std::size_t i=0;
   vector<float> propagated_Glb_x;
   propagated_Glb_x.clear();
   vector<float> propagated_Glb_y;
   propagated_Glb_y.clear();

   vector<float> matched_Glb_2D;//_clsz_min2;
   matched_Glb_2D.clear();//_clsz_min2.clear();
   vector<float> residual_x;//_clsz_min2;
   residual_x.clear();//_clsz_min2.clear();

   vector<float> matched_Glb_x;
   matched_Glb_x.clear();
   vector<float> matched_Glb_y;
   matched_Glb_y.clear();
   vector<float> matched_Loc_x;
   matched_Loc_x.clear();

   book();
   
      
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = fChain->LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      runNumber = std::to_string(event_runNumber);
     
      //for(std::size_t iSize = 0; iSize < gemRecHit_nRecHits; iSize++)
      //{
      //  if(gemRecHit_region->at(iSize) == -1) 
      //    {
      if(mu_propagatedGlb_x->size() != 0)// && mu_propagatedGlb_z->at(iSize) < 0)
	{
	  vector<float> matched_g_x = findMatchedHit(gemRecHit_g_x,gemRecHit_g_y, mu_propagatedGlb_x, mu_propagatedGlb_y,0,gemRecHit_region,mu_propagatedGlb_z);
	  vector<float> matched_g_y = findMatchedHit(gemRecHit_g_x,gemRecHit_g_y, mu_propagatedGlb_x, mu_propagatedGlb_y,1,gemRecHit_region,mu_propagatedGlb_z);
	  i=0;
	  if(matched_g_x.size() != 0)
	    {
	      while(i< matched_g_x.size())
		{
		  std::cout << matched_g_x.at(i) << std::endl;
		  matched_Glb_x.push_back(matched_g_x.at(i));
		  i++;
		}
	    }
	  i=0;
	  if(matched_g_y.size() != 0)
	    {
	      while(i< matched_g_y.size())
		{
		  matched_Glb_y.push_back(matched_g_y.at(i));
		  i++;
		}
	    }
	  
	  vector<float> matched_loc_x = findLocalMatchedHit(gemRecHit_loc_x, mu_propagatedLoc_x,gemRecHit_region,mu_propagatedLoc_z);
	  i=0;
	  if(matched_loc_x.size() != 0)
	    {
	      while(i< matched_loc_x.size())
		{
		  matched_Loc_x.push_back(matched_loc_x.at(i));
		  i++;
		}
	    }
	  
	  vector<float> res_x = residual(gemRecHit_loc_x, mu_propagatedLoc_x,gemRecHit_region,mu_propagatedLoc_z);
	  i=0;
	  if(res_x.size() != 0)
	    {
	      while(i< res_x.size()){
		residual_x.push_back(res_x.at(i));
		i++;
	      }
	    }
	}
      //}
      //}
      
      for(std::size_t iRecDigi = 0; iRecDigi < gemRecHit_nRecHits; ++ iRecDigi)
	{
	  TString region = std::to_string(gemRecHit_region->at(iRecDigi));
	  if(gemRecHit_region->at(iRecDigi) == -1)
	    {
	      if(gemRecHit_layer->at(iRecDigi) == 1)
		{
		  if(gemRecHit_chamber->at(iRecDigi) == 20) continue;
		  Double_t rec_x = gemRecHit_g_x->at(iRecDigi);
		  Double_t rec_y = gemRecHit_g_y->at(iRecDigi);
		  m_plots["OccupancyRecHit2D_layer1_endcapneg"]->Fill(rec_x,rec_y);
		  m_plots["OccupancyRecHit2D_layer1_endcapneg"]->SetTitle("OccupancyRecHit2D_layer1_endcap"+region+"_"+runNumber);
		}
	      else if(gemRecHit_layer->at(iRecDigi) == 2)
		{
		  if(gemRecHit_chamber->at(iRecDigi) == 20) continue;
		  Double_t rec_x = gemRecHit_g_x->at(iRecDigi);
		  Double_t rec_y = gemRecHit_g_y->at(iRecDigi);
		  m_plots["OccupancyRecHit2D_layer2_endcapneg"]->Fill(rec_x,rec_y);
		  m_plots["OccupancyRecHit2D_layer2_endcapneg"]->SetTitle("OccupancyRecHit2D_layer2_endcap"+region+"_"+runNumber);
		}
	    }
	  else if (gemRecHit_region->at(iRecDigi) == 1)
	    {
	      if(gemRecHit_layer->at(iRecDigi) == 1)
		{
		  Double_t rec_x = gemRecHit_g_x->at(iRecDigi);
		  Double_t rec_y = gemRecHit_g_y->at(iRecDigi);
		  m_plots["OccupancyRecHit2D_layer1_endcappos"]->Fill(rec_x,rec_y);
		  m_plots["OccupancyRecHit2D_layer1_endcappos"]->SetTitle("OccupancyRecHit2D_layer1_endcap"+region+"_"+runNumber);
		}
	      else if(gemRecHit_layer->at(iRecDigi) == 2)
		{
		  Double_t rec_x = gemRecHit_g_x->at(iRecDigi);
		  Double_t rec_y = gemRecHit_g_y->at(iRecDigi);
		  m_plots["OccupancyRecHit2D_layer2_endcappos"]->Fill(rec_x,rec_y);
		  m_plots["OccupancyRecHit2D_layer2_endcappos"]->SetTitle("OccupancyRecHit2D_layer2_endcap"+region+"_"+runNumber);
		}
	    }
	}

 
      if(mu_propagatedGlb_x->size() != 0)
	{
	  for(std::size_t iSize = 0; iSize < mu_propagatedGlb_x->size(); iSize++)  
	  {
	    if(mu_propagatedGlb_z->at(iSize) < 0)
	  	{
		  m_plots["OccupancyPropagated2D"]->Fill(mu_propagatedGlb_x->at(iSize), mu_propagatedGlb_y->at(iSize));
	  	}
	  }
	  m_plots["OccupancyPropagated2D"]->GetZaxis()->SetRangeUser(0,16500000);
	  m_plots["OccupancyPropagated2D"]->SetTitle("OccupancyPropagated2D_run"+runNumber);
	}
 
    
      if(mu_propagatedLoc_x->size() != 0)
        {
          for(std::size_t iSize = 0; iSize < mu_propagatedLoc_x->size(); iSize++)
	    {
	      TString region_loc = std::to_string(mu_propagated_region->at(iSize));
	      if(mu_propagated_region->at(iSize) ==  -1)
                {
		  Double_t loc_x = mu_propagatedLoc_x->at(iSize);
                  m_plots["OccupancyPropagatedLocalX"]->Fill(loc_x);
		  m_plots["OccupancyPropagatedLocalX"]->SetTitle("OccupancyPropagatedLocalX_endcap"+region_loc+"_run"+runNumber);
		}
	    }
        }

      /*m_plots["MatchedLocalEventNumber"]->Fill(event_eventNumber, matched_loc_x.size());
	
	for(std::size_t iRec = 0; iRec < gemRecHit_nRecHits; ++iRec)
	{
	Double_t loc_x = gemRecHit_loc_x->at(iRec);
	Double_t etaPartition = gemRecHit_etaPartition->at(iRec);
	m_plots["LocalXEtaPartition"]->Fill(loc_x,etaPartition);
	}
      */
   }
      
   for(std::size_t iSize=0; iSize< residual_x.size(); iSize++)
     {
       m_plots["ResidualPlotLocX"]->Fill(residual_x.at(iSize));
     }
   m_plots["ResidualPlotLocX"]->SetTitle("ResidualPlotLocX_run"+runNumber);
   
   for(std::size_t iSize=0; iSize < matched_Loc_x.size(); iSize++)
     {
       m_plots["OccupancyMatchedLocalX"]->Fill(matched_Loc_x.at(iSize));
     }
   m_plots["OccupancyMatchedLocalX"]->SetTitle("OccupancyMatchedLocalX_run"+runNumber);
   
   for(std::size_t iSize=0; iSize < matched_Glb_y.size(); iSize++)
     {
       m_plots["OccupancyMatchedGlobal2D"]->Fill(matched_Glb_x.at(iSize),matched_Glb_y.at(iSize));
     }
   m_plots["OccupancyMatchedGlobal2D"]->GetZaxis()->SetRangeUser(0,16500000);
   m_plots["OccupancyMatchedGlobal2D"]->SetTitle("OccupancyMatchedGlobal2D_run"+runNumber);
   
      
   m_plots["OccupancyRecHit2D_layer1_endcappos"]->GetZaxis()->SetRangeUser(0,16500000);
   m_plots["OccupancyRecHit2D_layer2_endcappos"]->GetZaxis()->SetRangeUser(0,16500000);
   m_plots["OccupancyRecHit2D_layer1_endcapneg"]->GetZaxis()->SetRangeUser(0,16500000);
   m_plots["OccupancyRecHit2D_layer2_endcapneg"]->GetZaxis()->SetRangeUser(0,16500000);

   m_plots["EfficiencyMatchedLocalX_endcapneg"] = (TH1F*)m_plots["OccupancyMatchedLocalX"]->Clone("EfficiencyMatchedLocalX_endcapneg");
   m_plots["EfficiencyMatchedLocalX_endcapneg"]->SetTitle("Efficiency Matched/Propagated LocalX_"+runNumber);
   m_plots["EfficiencyMatchedLocalX_endcapneg"]->Divide(m_plots["OccupancyMatchedLocalX"]);
         
   //m__plots["MatchedLocalEventNumber"]->SetTitle("MatchedLocalEventNumber_run"+runNumber);
   //m_plots["LocalXEtaPartition"]->SetTitle("LocalXEtaPartition_run"+runNumber);

   endJob();
   
}

vector<float> gemAnalysis::findLocalMatchedHit(vector<float> *recHitPositions, vector<float> *muPropagatedPositions, vector<int> *recHit_region, vector<float> *muPropagatedPositions_z)
{
  float min_recHitPosition = 666; //DUMMY VALUE
  vector<float> min_recHitPositions;

  for(std::size_t iMu = 0; iMu < muPropagatedPositions->size(); iMu++)
    {
      if(muPropagatedPositions_z->at(iMu) < 0)
        {
	  float min_residual = 10;
	  for(std::size_t iRecHit = 0; iRecHit < recHitPositions->size(); iRecHit++)
	    {
	      if(recHit_region->at(iRecHit) == -1)
		{
		  float residual = std::fabs(muPropagatedPositions->at(iMu)-recHitPositions->at(iRecHit));
		  if(residual <= min_residual)
		    {
		      min_residual = residual;
		      {
			min_recHitPosition = recHitPositions->at(iRecHit);
		      }
		    }
		  else
		    {
		      continue;
		    }
		}
	    }
	  if(min_recHitPosition == 666) continue;
	  else
	    {
	      min_recHitPositions.push_back(min_recHitPosition);
	    }
	}
    }
  
  return min_recHitPositions;
      
}

vector<float> gemAnalysis::findMatchedHit(vector<float> *recHitPositions_x, vector<float> *recHitPositions_y, vector<float> *muPropagatedPositions_x, vector<float> *muPropagatedPositions_y, int coord, vector<int> *recHit_region, vector<float> *muPropagatedPositions_z)
{

  float min_recHitPosition = 666; //DUMMY VALUE
  vector<float> min_recHitPositions;
  //recHitMatched.clear();
    
  for(std::size_t iMu = 0; iMu < muPropagatedPositions_x->size(); iMu++)
    {
      if(muPropagatedPositions_z->at(iMu) < 0)
	{
	  float min_residual = 10;
	  for(std::size_t iRecHit = 0; iRecHit < recHitPositions_x->size(); iRecHit++)
	    {
	      if(recHit_region->at(iRecHit) == -1)
		{
		  float residualx = std::fabs(muPropagatedPositions_x->at(iMu))-std::fabs(recHitPositions_x->at(iRecHit));
		  float residualy = std::fabs(muPropagatedPositions_y->at(iMu))-std::fabs(recHitPositions_y->at(iRecHit));
		  float residual = TMath::Sqrt(TMath::Power(residualx,2)+TMath::Power(residualy,2));
		  if(residual <= min_residual)
		    {
		      min_residual = residual;
		      if (coord == 0)
			{
			  min_recHitPosition = recHitPositions_x->at(iRecHit);
			}
		      else if (coord == 1)
			{
			  min_recHitPosition = recHitPositions_y->at(iRecHit);
			}
		    }
		  else
		    { 
		      continue;
		    }
		}
	    }
	  if(min_recHitPosition == 666) continue;
	  else
	    {
	      min_recHitPositions.push_back(min_recHitPosition);
	    }
	}
    }
  
  return min_recHitPositions;
  
}


vector<float> gemAnalysis::residual2D(vector<float> *recHitPositions_x, vector<float> *recHitPositions_y, vector<float> *muPropagatedPositions_x, vector<float> *muPropagatedPositions_y)
{

  float min_recHitPosition = 666; //DUMMY VALUE
  vector<float> residuals;
  //recHitMatched.clear();

  for(std::size_t iMu = 0; iMu < muPropagatedPositions_x->size(); iMu++)
    {
      float min_residual = 100;
      for(std::size_t iRecHit = 0; iRecHit < recHitPositions_x->size(); iRecHit++)
        {
          float residualx = std::fabs(muPropagatedPositions_x->at(iMu)-recHitPositions_x->at(iRecHit));
	  float residualy = std::fabs(muPropagatedPositions_y->at(iMu)-recHitPositions_y->at(iRecHit));
          float residual = TMath::Sqrt(TMath::Power(residualx,2)+TMath::Power(residualy,2));
	  std::cout << "residual" <<std::endl;
          if(std::fabs(residual) <= std::fabs(min_residual))
            {
              min_residual = residual;
              min_recHitPosition = recHitPositions_x->at(iRecHit);
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

vector<float> gemAnalysis::residual(vector<float> *recHitPositions, vector<float> *muPropagatedPositions, vector<int> *recHit_region, vector<float> *muPropagatedPositions_z)
{

float min_recHitPosition = 666; //DUMMY VALUE
vector<float> residuals;

for(std::size_t iMu = 0; iMu < muPropagatedPositions->size(); iMu++)
  {
    if(muPropagatedPositions_z->at(iMu) < 0)
      {
	float min_residual = 10;
	for(std::size_t iRecHit = 0; iRecHit < recHitPositions->size(); iRecHit++)
	  {
	    if(recHit_region->at(iRecHit) == -1)
	      {
		float residual = std::fabs(muPropagatedPositions->at(iMu)-std::fabs(recHitPositions->at(iRecHit));
		if(residual <= min_residual)
		  {
		    min_residual = residual;
		    min_recHitPosition = recHitPositions->at(iRecHit);
		  }
		else
		  {
		    continue;
		  }
	      }
	  }
	if(min_recHitPosition == 666) continue;
	else
	  {
	    residuals.push_back(min_residual);
	  }
      }
  }
 
 return residuals;
 
}


void gemAnalysis::book()
{
  m_outFile.cd();

  m_plots["OccupancyPropagatedLocalX"] = new TH1F("OccupancyPropagatedLocalX",
					       "Propagated Local X; localX ; entries",
					       100,-25.,25.);

     
  m_plots["ResidualPlotLocX"] = new TH1F("ResidualPlotLocX",
				       "Residual LocalX; residualX ; entries",
				       400,-2,2.);
  
  
  
  m_plots["OccupancyPropagated2D"] = new TH2F("OccupancyPropagated2D",
					      "Propagated 2D; globalX ; globalY",
                                              100,-250.,250.,
                                              100,-250.,250.);

  m_plots["OccupancyMatchedGlobal2D"] = new TH2F("OccupancyMatchedGlobal2D",
					   "Matched Global 2D; globalX ; globalY",
					   100,-250.,250.,
					   100,-250.,250.);

  m_plots["OccupancyMatchedLocalX"] = new TH1F("OccupancyMatchedLocalX",
					      "Matched Local X; localX ; entries",
					      100,-25.,25.);
  

  m_plots["OccupancyRecHit2D_layer1_endcappos"] = new TH2F("OccupancyRecHit2D_layer1_endcappos",
							 "RecHit 2D GE11_layer1; globalX ; globalY",
							 100,-250.,250.,
							 100,-250.,250.);
  
  m_plots["OccupancyRecHit2D_layer2_endcappos"] = new TH2F("OccupancyRecHit2D_layer2_endcappos",
                                                         "RecHit 2D GE11_layer2; globalX ; globalY",
                                                         100,-250.,250.,
                                                         100,-250.,250.);
  
  m_plots["OccupancyRecHit2D_layer1_endcapneg"] = new TH2F("OccupancyRecHit2D_layer1_endcapneg",
                                                         "RecHit 2D GE11_layer1; globalX ; globalY",
                                                         100,-250.,250.,
                                                         100,-250.,250.);

  m_plots["OccupancyRecHit2D_layer2_endcapneg"] = new TH2F("OccupancyRecHit2D_layer2_endcapneg",
                                                         "RecHit 2D GE11_layer2; globalX ; globalY",
                                                         100,-250.,250.,
                                                         100,-250.,250.);
  
   
  m_plots["MatchedLocalEventNumber"] = new TH2F("MatchedLocalEventNumber",
					   "MatchedLocal vs EventNumber; eventNumber; # matched",
						30000000,0,30000000,
						100, 0., 100.);
 
  m_plots["LocalXEtaPartition"] = new TH2F("LocalXEtaPartition",
					   "LocalX vs EtaPartition; localX; etaPartition",
					   100,-20.,20.,
					   12,0.,12.);

  m_plots["EfficiencyMatchedLocalX_endcappos"] = new TH1F("EfficiencyMatchedLocalX_endcappos",
							  "Efficiency Matched Local X; localX ; efficiency",
							  100,-25.,25.);
  
  m_plots["EfficiencyMatchedLocalX_endcapneg"] = new TH1F("EfficiencyMatchedLocalX_endcapneg",
							  "Efficiency Matched Local X; localX ; efficiency",
							  100,-25.,25.);
  
  //m_plots["EfficiencyMatchedGlobal_endcappos"] = new TH2F("EfficiencyMatchedGlobal_endcappos",
  //							  "Efficiency Matched Global; globalX ; entries",
  //							  100,-250.,250.);


}


void gemAnalysis::endJob()
{

  m_outFile.cd();
  m_outFile.Write();
  m_outFile.Close();

}

