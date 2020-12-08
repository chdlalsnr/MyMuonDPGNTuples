#define gemAnalysis_cxx
#include "gemAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TCut.h>
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
//void gemAnalysis::Loop(const TString & inFileName,
//		       const TString & outFileName)
  
{
  //TFile *f_in = new TFile(inFileName,"READ");
  //TFile *f_out = new TFile(outFileName,"RECREATE");
  //fChain = 0 ;
  //fChain = static_cast<TTree*>(m_inFile.Get("muNtupleProducer/MuDPGTree"));
  Init(fChain);
  
  if (fChain == 0) return;
  
  //TFile *f1 = new TFile("/lustre/cms/store/user/gmilella/ExpressCosmics/CRAB3_gem_dpg_ntuple_mwgr5_run338714_express/201203_100335/0000/MuDPGNTuple_MWGR5_EXP_run338714.root","READ");
  //TFile *f1 = new TFile("MuDPGNTuple_MWGR5_EXP_run338714.root","READ");
  std::cout << "File opened" << std::endl;
  
   Long64_t nentries = fChain->GetEntries();
   
   book();
      
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = fChain->LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      runNumber = std::to_string(event_runNumber);

      	  
      findMatchedHit(gemRecHit_g_x,gemRecHit_g_y, mu_propagatedGlb_x, mu_propagatedGlb_y, gemRecHit_region, mu_propagated_region, mu_propagated_pt, mu_propagated_eta,gemRecHit_layer,gemRecHit_chamber, mu_propagated_layer, mu_propagated_chamber ,mu_isME11);
      findLocalMatchedHit(gemRecHit_loc_x, mu_propagatedLoc_x, gemRecHit_region, mu_propagated_region, mu_propagated_pt, mu_propagated_eta, gemRecHit_layer, gemRecHit_chamber, mu_propagated_layer , mu_propagated_chamber, gemRecHit_etaPartition, mu_propagated_etaP, mu_isME11);

      /*if(gemRecHit_nRecHits->size() != 0)
	{
	  for(std::size_t iMu = 0; iMu < mu_propagatedLoc_x->size(); iMu++)
	    {
	      if(mu_propagated_region->at(iMu) == -1 && !(mu_propagated_chamber->at(iMu) == 20 && mu_propagated_layer->at(iMu) == 1))
		{
		  m_plots["MuPtPropagated"]->Fill(mu_propagated_pt->at(iMu));
		  m_plots["MuEtaPropagated"]->Fill(std::fabs(mu_propagated_eta(iMu)));

		  for(std::size_t iRecHit = 0; iRecHit < gemRecHit_nRechit->size(); iRecHit++)
		    {
		      if(gemRecHit_region->at(iRecHit) == -1 && !(gemRecHit_chamber->at(iRecHit) == 20 && gemRecHit_layer->at(iRecHit) == 1))
			{
			  //matching locale
			  float residual_x = std::fabs(mu_propagatedLoc_x->at(iMu)-gemRecHit_loc_x->at(iRecHit));
			  int same_eta = mu_propagatedLoc_eta->at(iMu) - gemRecHit_etaPartition->at(iRecHit);
			  bool same_chamber = (mu_propagated_chamber->at(iMu) == gemRecHit_chamber->at(iRecHit));
			  if(residual_x < min_residualx && same_eta == 0 && same_chamber == 1)
			    {
			      
			    }
			  else continue;
			}
		      
		    }
		}
	    }
	}*/
      

   }

   m_plots["ResidualPlotGlobal_endcappos"]->SetTitle("ResidualPlotGlobal_endcap#plus_run"+runNumber);
   m_plots["ResidualPlotGlobal_endcapneg"]->SetTitle("ResidualPlotGlobal_endcap#minus_run"+runNumber);
   m_plots["ResidualPlotLocalX_endcappos"]->SetTitle("ResidualPlotLocalX_endcap#plus_run"+runNumber);
   m_plots["ResidualPlotLocalX_endcapneg"]->SetTitle("ResidualPlotLocalX_endcap#minus_run"+runNumber);
   m_plots["OccupancyMatchedLocalX_endcappos"]->SetTitle("OccupancyMatchedLocalX_endcap#plus_run"+runNumber);
   m_plots["OccupancyMatchedLocalX_endcapneg"]->SetTitle("OccupancyMatchedLocalX_endcap#minus_run"+runNumber);
   m_plots["OccupancyMatchedGlobal2D_endcapneg"]->SetTitle("OccupancyMatchedGlobal2D_endcap#minus_run"+runNumber);
   m_plots["OccupancyMatchedGlobal2D_endcappos"]->SetTitle("OccupancyMatchedGlobal2D_endcap#plus_run"+runNumber);
   m_plots["MuPtMatchedLocal_endcapneg"]->SetTitle("MuPtMatchedLocal_endcapneg_endcap#minus"+runNumber);
   m_plots["MuPtMatchedLocal_endcappos"]->SetTitle("MuPtMatchedLocal_endcappos_endcap#plus"+runNumber);
   m_plots["MuEtaMatchedLocal_endcapneg"]->SetTitle("MuEtaMatchedLocal_endcapneg_endcap#minus"+runNumber);
   m_plots["MuEtaMatchedLocal_endcappos"]->SetTitle("MuEtaMatchedLocal_endcappos_endcap#plus"+runNumber);
   m_plots["MuPtMatchedGlobal_endcapneg"]->SetTitle("MuPtMatchedGlobal_endcapneg_endcap#minus"+runNumber);
   m_plots["MuPtMatchedGlobal_endcappos"]->SetTitle("MuPtMatchedGlobal_endcappos_endcap#plus"+runNumber);
   m_plots["MuEtaMatchedGlobal_endcapneg"]->SetTitle("MuPtMatchedGlobal_endcappos_endcap#minus"+runNumber);
   m_plots["MuEtaMatchedGlobal_endcappos"]->SetTitle("MuPtMatchedGlobal_endcappos_endcap#plus"+runNumber);

      
   TString var[] = {"mu_propagated_pt", "mu_propagated_eta", "mu_propagatedLoc_x", "mu_propagatedGlb_x", "gemRecHit_g_x"};

   TTree *tTree = (TTree*)m_inFile.Get("muNtupleProducer/MuDPGTree");
   
   int n = sizeof(var)/sizeof(var[0]);
   int nBins[n];
   int nEntries[n];
   TH1F *histo_neg[n];
   TH1F *histo_pos[n];
   TH2F *histo2D_neg[n];
   TH2F *histo2D_pos[n];
   TString binning;

   for(int i=0; i<n; i++){
     TString varname = var[i];
     TString s = std::to_string(i);

     if(varname=="mu_propagated_pt") binning = "(200,0,100.)";
     if(varname=="mu_propagated_eta") binning = "(300,0.,3.)";
     if(varname=="mu_propagatedLoc_x") binning = "(100,-25.,25.,8,1.,9.)"; 
     if(varname=="mu_propagatedGlb_x") binning = "(100,-250.,250.,100,-250.,250)";
     if(varname=="gemRecHit_g_x") binning = "(100,-250.,250.,100,-250.,250)";
     
     if(varname=="mu_propagatedGlb_x")
       {
         //TCut c1 = "gemRecHit_layer && !(gemRecHit_layer == 1 && gemRecHit_chamber != 20)";
         //TCut c2 = "!gemRecHit_layer && (gemRecHit_layer == 1 && gemRecHit_chamber != 20)";
         //TCut c3 = c1 || c2;
         tTree->Draw("mu_propagatedGlb_y:"+varname+">>histo2D_neg0"+binning,"mu_propagated_region == -1 && !(mu_propagated_layer == 1 && mu_propagated_chamber == 20) && mu_isME11");
         tTree->Draw("mu_propagatedGlb_y:"+varname+">>histo2D_pos0"+binning,"mu_propagated_region == 1 && mu_isME11");
       }
     
     if(varname=="gemRecHit_g_x")
       {
	 tTree->Draw("gemRecHit_g_y:"+varname+">>histo2D_neg1"+binning,"gemRecHit_region == -1 && gemRecHit_layer == 1");//,"colz");
	 tTree->Draw("gemRecHit_g_y:"+varname+">>histo2D_pos1"+binning,"gemRecHit_region == 1");//,"colz");
	 tTree->Draw("gemRecHit_g_y:"+varname+">>histo2D_neg2"+binning,"gemRecHit_region == -1 && gemRecHit_layer == 1 && gemRecHit_chamber != 20");//,"colz");
	 tTree->Draw("gemRecHit_g_y:"+varname+">>histo2D_neg3"+binning,"gemRecHit_region == -1 && gemRecHit_layer == 2");//,"colz");
	 //tTree->Draw("gemRecHit_g_y:"+varname+">>histo2D_pos2"+binning,"gemRecHit_region == 1","colz");
       }
     
     if(varname=="mu_propagatedLoc_x")
       {
	 tTree->Draw("mu_propagated_etaP:"+varname+">>histo2D_neg4"+binning,"mu_propagated_region == -1 && !(mu_propagated_layer == 1 && mu_propagated_chamber == 20) && mu_isME11");
         tTree->Draw("mu_propagated_etaP:"+varname+">>histo2D_pos4"+binning,"mu_propagated_region == 1 && mu_isME11");
       }

     if(varname=="mu_propagated_eta")
       {
	 tTree->Draw("TMath::Abs("+varname+")>>histo_neg"+s+binning,"mu_propagated_region == -1 && !(mu_propagated_layer == 1 && mu_propagated_chamber == 20)  && mu_isME11");
         tTree->Draw("TMath::Abs("+varname+")>>histo_pos"+s+binning,"mu_propagated_region == 1  && mu_isME11");
       }

     else
       {
	 //TCut c1 = "gemRecHit_layer && !(gemRecHit_layer == 1 && gemRecHit_chamber != 20)";
	 //TCut c2 = "!gemRecHit_layer && (gemRecHit_layer == 1 && gemRecHit_chamber != 20)";
	 //TCut c3 = c1 || c2;
	 tTree->Draw(varname+">>histo_neg"+s+binning,"mu_propagated_region == -1 && !(mu_propagated_layer == 1 && mu_propagated_chamber == 20)  && mu_isME11");
	 tTree->Draw(varname+">>histo_pos"+s+binning,"mu_propagated_region == 1  && mu_isME11");
       }
     
     histo_neg[i] = (TH1F*)gDirectory->Get("histo_neg"+s);
     histo_pos[i] = (TH1F*)gDirectory->Get("histo_pos"+s);
     histo2D_neg[0] = (TH2F*)gDirectory->Get("histo2D_neg0");
     histo2D_pos[0] = (TH2F*)gDirectory->Get("histo2D_pos0");
     histo2D_neg[1] = (TH2F*)gDirectory->Get("histo2D_neg1");
     histo2D_pos[1] = (TH2F*)gDirectory->Get("histo2D_pos1");
     //histo2D_pos[2] = (TH2F*)gDirectory->Get("histo2D_pos2");
     histo2D_neg[2] = (TH2F*)gDirectory->Get("histo2D_neg2");
     histo2D_neg[3] = (TH2F*)gDirectory->Get("histo2D_neg3");
     histo2D_neg[4] = (TH2F*)gDirectory->Get("histo2D_neg4");
     histo2D_pos[4] = (TH2F*)gDirectory->Get("histo2D_pos4");
   }
  
   m_plots["OccupancyRecHit2D_endcapneg_layer1with20"] = (TH2F*)histo2D_neg[1]->Clone("OccupancyRecHit2D_endcapneg_layer1with20");
   m_plots["OccupancyRecHit2D_endcapneg_layer1with20"]->SetTitle("OccupancyRecHit2D_endcap#minus_layer1"+runNumber);
   m_plots["OccupancyRecHit2D_endcapneg_layer1with20"]->GetYaxis()->SetTitle("globalY");
   m_plots["OccupancyRecHit2D_endcapneg_layer1with20"]->GetYaxis()->SetTitleOffset(1);
   m_plots["OccupancyRecHit2D_endcapneg_layer1with20"]->GetXaxis()->SetTitle("globalX");
   m_plots["OccupancyRecHit2D_endcapneg_layer1"] = (TH2F*)histo2D_neg[2]->Clone("OccupancyRecHit2D_endcapneg_layer1");
   m_plots["OccupancyRecHit2D_endcapneg_layer1"]->SetTitle("OccupancyRecHit2D_endcap#minus_layer1"+runNumber);
   m_plots["OccupancyRecHit2D_endcapneg_layer1"]->GetYaxis()->SetTitle("globalY");
   m_plots["OccupancyRecHit2D_endcapneg_layer1"]->GetYaxis()->SetTitleOffset(1);
   m_plots["OccupancyRecHit2D_endcapneg_layer1"]->GetXaxis()->SetTitle("globalX");

   m_plots["OccupancyRecHit2D_endcapneg_layer2"] = (TH2F*)histo2D_neg[3]->Clone("OccupancyRecHit2D_endcapneg_layer2");
   m_plots["OccupancyRecHit2D_endcapneg_layer2"]->SetTitle("OccupancyRecHit2D_endcap#minus_layer2"+runNumber);
   m_plots["OccupancyRecHit2D_endcapneg_layer2"]->GetYaxis()->SetTitle("globalY");
   m_plots["OccupancyRecHit2D_endcapneg_layer2"]->GetYaxis()->SetTitleOffset(1);
   m_plots["OccupancyRecHit2D_endcapneg_layer2"]->GetXaxis()->SetTitle("globalX");

   m_plots["OccupancyRecHit2D_endcappos"] = (TH2F*)histo2D_pos[1]->Clone("OccupancyRecHit2D_endcappos");
   m_plots["OccupancyRecHit2D_endcappos"]->SetTitle("OccupancyRecHit2D_endcap#plus"+runNumber);
   m_plots["OccupancyRecHit2D_endcappos"]->GetYaxis()->SetTitle("globalY");
   m_plots["OccupancyRecHit2D_endcappos"]->GetYaxis()->SetTitleOffset(1);
   m_plots["OccupancyRecHit2D_endcappos"]->GetXaxis()->SetTitle("globalX");
  
   m_plots["OccupancyPropagatedLocalX_endcapneg"] = (TH2F*)histo2D_neg[4]->Clone("OccupancyPropagatedLocalX_endcapneg");
   m_plots["OccupancyPropagatedLocalX_endcapneg"]->SetTitle("OccupancyPropagatedLocalX_endcap#minus"+runNumber);
   m_plots["OccupancyPropagatedLocalX_endcapneg"]->GetYaxis()->SetTitle("globalY");
   m_plots["OccupancyPropagatedLocalX_endcapneg"]->GetYaxis()->SetTitleOffset(1);
   m_plots["OccupancyPropagatedLocalX_endcapneg"]->GetXaxis()->SetTitle("globalX");
   m_plots["OccupancyPropagatedLocalX_endcappos"] = (TH2F*)histo2D_pos[4]->Clone("OccupancyPropagatedLocalX_endcappos");
   m_plots["OccupancyPropagatedLocalX_endcappos"]->SetTitle("OccupancyPropagatedLocalX_endcap#plus"+runNumber);
   m_plots["OccupancyPropagatedLocalX_endcappos"]->GetYaxis()->SetTitle("globalY");
   m_plots["OccupancyPropagatedLocalX_endcappos"]->GetYaxis()->SetTitleOffset(1);
   m_plots["OccupancyPropagatedLocalX_endcappos"]->GetXaxis()->SetTitle("globalX");
   
   m_plots["OccupancyPropagated2D_endcapneg"] = (TH2F*)histo2D_neg[0]->Clone("OccupancyPropagated2D_endcapneg");
   m_plots["OccupancyPropagated2D_endcapneg"]->SetTitle("OccupancyPropagated2D_endcap#minus"+runNumber);
   m_plots["OccupancyPropagated2D_endcapneg"]->GetYaxis()->SetTitle("globalY");
   m_plots["OccupancyPropagated2D_endcapneg"]->GetYaxis()->SetTitleOffset(1);
   m_plots["OccupancyPropagated2D_endcapneg"]->GetXaxis()->SetTitle("globalX");
   m_plots["OccupancyPropagated2D_endcappos"] = (TH2F*)histo2D_pos[0]->Clone("OccupancyPropagated2D_endcappos");
   m_plots["OccupancyPropagated2D_endcappos"]->SetTitle("OccupancyPropagated2D_endcap#plus"+runNumber);
   m_plots["OccupancyPropagated2D_endcappos"]->GetYaxis()->SetTitle("globalY");
   m_plots["OccupancyPropagated2D_endcappos"]->GetYaxis()->SetTitleOffset(1);
   m_plots["OccupancyPropagated2D_endcappos"]->GetXaxis()->SetTitle("globalX");


   /*m_plots["EfficiencyMatchedLocalX_endcapneg"] = (TH1F*)m_plots["OccupancyMatchedLocalX_endcapneg"]->Clone("EfficiencyMatchedLocalX_endcapneg");
   m_plots["EfficiencyMatchedLocalX_endcapneg"]->SetTitle("Efficiency Matched/Propagated LocalX_endcap#minus_"+runNumber);
   m_plots["EfficiencyMatchedLocalX_endcapneg"]->Sumw2();
   m_plots["EfficiencyMatchedLocalX_endcapneg"]->Divide(m_plots["OccupancyPropagatedLocalX_endcapneg"]);*/
   /*m_plots["EfficiencyMatchedLocalX_endcappos"] = (TH1F*)m_plots["OccupancyMatchedLocalX_endcappos"]->Clone("EfficiencyMatchedLocalX_endcappos");
   m_plots["EfficiencyMatchedLocalX_endcappos"]->SetTitle("Efficiency Matched/Propagated LocalX_endcap#plus"+runNumber);
   m_plots["EfficiencyMatchedLocalX_endcappos"]->Sumw2();
   m_plots["EfficiencyMatchedLocalX_endcappos"]->Divide(m_plots["OccupancyPropagatedLocalX_endcappos"]);*/
 
   /*m_plots["EfficiencyMatchedGlobal2D_endcapneg"] = (TH1F*)m_plots["OccupancyMatchedGlobal2D_endcapneg"]->Clone("EfficiencyMatchedGlobal2D_endcapneg");
   m_plots["EfficiencyMatchedGlobal2D_endcapneg"]->SetTitle("Efficiency Matched/Propagated Global2D_endcap#minus_"+runNumber);
   m_plots["EfficiencyMatchedGlobal2D_endcapneg"]->Sumw2();
   m_plots["EfficiencyMatchedGlobal2D_endcapneg"]->Divide(m_plots["OccupancyPropagated2D_endcapneg"]);*/
   /*m_plots["EfficiencyMatchedGlobal2D_endcappos"] = (TH1F*)m_plots["OccupancyMatchedGlobal2D_endcappos"]->Clone("EfficiencyMatchedGlobal2D_endcappos");
   m_plots["EfficiencyMatchedGlobal2D_endcappos"]->SetTitle("Efficiency Matched/Propagated Global2D_endcap#plus"+runNumber);
   m_plots["EfficiencyMatchedGlobal2D_endcappos"]->Sumw2();
   m_plots["EfficiencyMatchedGlobal2D_endcappos"]->Divide(m_plots["OccupancyPropagated2D_endcappos"]);*/
     
   //m_plots["EfficiencyPtMatchedLocalX_endcapneg"] = (TH1F*)m_plots["MuPtMatchedLocal_endcapneg"]->Clone("EfficiencyPtMatchedLocalX_endcapneg");
   m_plots["EfficiencyPtMatchedLocalX_endcapneg"] = (TH1F*)m_plots["MuPtMatchedLocal_endcapneg"]->Clone("EfficiencyPtMatchedLocalX_endcapneg");
   m_plots["EfficiencyPtMatchedLocalX_endcapneg"]->SetTitle("Efficiency Pt Local_endcap#minus_"+runNumber);
   m_plots["EfficiencyPtMatchedLocalX_endcapneg"]->Sumw2();
   m_plots["EfficiencyPtMatchedLocalX_endcapneg"]->Divide(histo_neg[0]);
   /*m_plots["EfficiencyPtMatchedLocalX_endcappos"] = (TH1F*)m_plots["MuPtMatchedLocal_endcappos"]->Clone("EfficiencyPtMatchedLocalX_endcappos");
   m_plots["EfficiencyPtMatchedLocalX_endcappos"]->SetTitle("Efficiency Pt Local_endcap#plus_"+runNumber);
   m_plots["EfficiencyPtMatchedLocalX_endcappos"]->Sumw2();
   m_plots["EfficiencyPtMatchedLocalX_endcappos"]->Divide(histo_pos[0]);*/
   
   m_plots["EfficiencyEtaMatchedLocalX_endcapneg"] = (TH1F*)m_plots["MuEtaMatchedLocal_endcapneg"]->Clone("EfficiencyEtaMatchedLocalX_endcapneg");
   m_plots["EfficiencyEtaMatchedLocalX_endcapneg"]->SetTitle("Efficiency Eta Local_endcap#minus_"+runNumber);
   m_plots["EfficiencyEtaMatchedLocalX_endcapneg"]->Sumw2();
   m_plots["EfficiencyEtaMatchedLocalX_endcapneg"]->Divide(histo_neg[1]);
   m_plots["EfficiencyEtaMatchedLocalX_endcapneg"]->GetXaxis()->SetRangeUser(1.50,2.20);
   /*m_plots["EfficiencyEtaMatchedLocalX_endcappos"] = (TH1F*)m_plots["MuEtaMatchedLocal_endcappos"]->Clone("EfficiencyEtaMatchedLocalX_endcappos");
   m_plots["EfficiencyEtaMatchedLocalX_endcappos"]->SetTitle("Efficiency Eta Local_endcap#plus_"+runNumber);
   m_plots["EfficiencyEtaMatchedLocalX_endcappos"]->Sumw2();
   m_plots["EfficiencyEtaMatchedLocalX_endcappos"]->Divide(histo_pos[1]);*/

   m_plots["EfficiencyPtMatchedGlobal2D_endcapneg"] = (TH1F*)m_plots["MuPtMatchedLocal_endcapneg"]->Clone("EfficiencyPtMatchedGlobal2D_endcapneg");
   m_plots["EfficiencyPtMatchedGlobal2D_endcapneg"]->SetTitle("Efficiency Pt Global_endcap#minus_"+runNumber);
   m_plots["EfficiencyPtMatchedGlobal2D_endcapneg"]->Sumw2();
   m_plots["EfficiencyPtMatchedGlobal2D_endcapneg"]->Divide(histo_neg[0]);
   /*m_plots["EfficiencyPtMatchedGlobal2D_endcappos"] = (TH1F*)m_plots["MuPtMatchedGlobal_endcappos"]->Clone("EfficiencyPtMatchedGlobal2D_endcappos");
   m_plots["EfficiencyPtMatchedGlobal2D_endcappos"]->SetTitle("Efficiency Pt Global_endcap#plus_"+runNumber);
   m_plots["EfficiencyPtMatchedGlobal2D_endcappos"]->Sumw2();
   m_plots["EfficiencyPtMatchedGlobal2D_endcappos"]->Divide(histo_pos[0]);*/

   m_plots["EfficiencyEtaMatchedGlobal2D_endcapneg"] = (TH1F*)m_plots["MuEtaMatchedLocal_endcapneg"]->Clone("EfficiencyEtaMatchedGlobal2D_endcapneg");
   m_plots["EfficiencyEtaMatchedGlobal2D_endcapneg"]->SetTitle("Efficiency Eta Global_endcap#minus_"+runNumber);
   m_plots["EfficiencyEtaMatchedGlobal2D_endcapneg"]->Sumw2();
   m_plots["EfficiencyEtaMatchedGlobal2D_endcapneg"]->Divide(histo_neg[1]);
   m_plots["EfficiencyEtaMatchedGlobal2D_endcapneg"]->GetXaxis()->SetRangeUser(1.50,2.20);
   /*m_plots["EfficiencyEtaMatchedGlobal2D_endcappos"] = (TH1F*)m_plots["MuEtaMatchedGlobal_endcappos"]->Clone("EfficiencyEtaMatchedGlobal2D_endcappos");
   m_plots["EfficiencyEtaMatchedGlobal2D_endcappos"]->SetTitle("Efficiency Eta Global_endcap#plus_"+runNumber);
   m_plots["EfficiencyEtaMatchedGlobal2D_endcappos"]->Sumw2();
   m_plots["EfficiencyEtaMatchedGlobal2D_endcappos"]->Divide(histo_pos[1]);*/

         
   //m__plots["MatchedLocalEventNumber"]->SetTitle("MatchedLocalEventNumber_run"+runNumber);
   //m_plots["LocalXEtaPartition"]->SetTitle("LocalXEtaPartition_run"+runNumber);

   endJob();
   
}

void gemAnalysis::findLocalMatchedHit(vector<float> *recHitPositions, vector<float> *muPropagatedPositions, vector<int> *recHitregion, vector<int> *muPropagatedRegion, vector<float> *muPt, vector<float> *muEta, vector<int> *recHitlayer, vector<int> *recHitchamber, vector<int> *muLayer, vector<int> *muChamber, vector<int> *recHitEta, vector<int> *muPropagatedEta, vector<bool> *muME11)
{
  float min_recHitPosition = 666; //DUMMY VALUE
  vector<float> min_recHitPositions;
  float mupt;
  float mueta;
  float etaP;
  float min_residual = 10;

  for(std::size_t iMu = 0; iMu < muPropagatedPositions->size(); iMu++)
    {
      if( muPropagatedRegion->at(iMu) == -1 && !(muChamber->at(iMu) == 20 && muLayer->at(iMu) == 1) && muME11->at(iMu) == 1)
	{
	  for(std::size_t iRecHit = 0; iRecHit < recHitPositions->size(); iRecHit++)
		{
		  if( recHitregion->at(iRecHit) == -1)
		    {
		      
		      if(recHitlayer->at(iRecHit) == 1 && recHitchamber->at(iRecHit) == 20) continue;
		      float residual_eta = std::fabs(muPropagatedEta->at(iMu)-recHitEta->at(iRecHit));
		      float residual = (muPropagatedPositions->at(iMu)-recHitPositions->at(iRecHit));
		      bool same_chamber = muChamber->at(iMu) == recHitchamber->at(iRecHit);
		      if(std::fabs(residual) < min_residual && residual_eta == 0 && same_chamber)
			{
			  min_residual = residual;
			  min_recHitPosition = recHitPositions->at(iRecHit);
			  //if(recHitchamber->(atiRecHit) == 8 || recHitchamber->(atiRecHit) == 9) continue;
			  etaP = recHitEta->at(iRecHit);
			  mupt = muPt->at(iMu);
			  mueta = TMath::Abs(muEta->at(iMu));
			}
		      else continue;
		    }
		}
	  if(min_recHitPosition == 666) continue;
	  else
	    {
	      m_plots["ResidualPlotLocalX_endcapneg"]->Fill(min_residual);
	      m_plots["OccupancyMatchedLocalX_endcapneg"]->Fill(min_recHitPosition,etaP);
	      m_plots["MuPtMatchedLocal_endcapneg"]->Fill(mupt);
	      m_plots["MuEtaMatchedLocal_endcapneg"]->Fill(mueta);
	    }
	 
	}
      
      else if( muPropagatedRegion->at(iMu) == 1 && muME11->at(iMu) == 1)
	{
	  for(std::size_t iRecHit = 0; iRecHit < recHitPositions->size(); iRecHit++)
	    {
	      if( recHitregion->at(iRecHit) == 1)
		{
		    float residual_eta = std::fabs(muPropagatedEta->at(iMu)-recHitEta->at(iRecHit));
		    float residual = (muPropagatedPositions->at(iMu)-recHitPositions->at(iRecHit));
		    if(std::fabs(residual) < min_residual && residual_eta == 0)
		      {
			min_residual = residual;
			min_recHitPosition = recHitPositions->at(iRecHit);
			etaP = recHitEta->at(iRecHit);
			mupt = muPt->at(iMu);
			mueta = muEta->at(iMu);
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
	      m_plots["ResidualPlotLocalX_endcappos"]->Fill(min_residual);
	      m_plots["OccupancyMatchedLocalX_endcappos"]->Fill(min_recHitPosition,etaP);
	      m_plots["MuPtMatchedLocal_endcappos"]->Fill(mupt);
	      m_plots["MuEtaMatchedLocal_endcappos"]->Fill(mueta);
	    }
	}
    }
}

void gemAnalysis::findMatchedHit(vector<float> *recHitPositions_x, vector<float> *recHitPositions_y, vector<float> *muPropagatedPositions_x, vector<float> *muPropagatedPositions_y, vector<int> *recHitregion, vector<int> *muPropagatedRegion, vector<float> *muPt, vector<float> *muEta, vector<int> *recHitlayer, vector<int> *recHitchamber, vector<int> *muLayer, vector<int> *muChamber , vector<bool> *isME11)
{

  float min_recHitPosition_x = 666; //DUMMY VALUE
  float min_recHitPosition_y = 666;
  float mupt,mueta;
  float min_residual = 10;
  vector<float> min_recHitPositions;
  //recHitMatched.clear();
  
  for(std::size_t iMu = 0; iMu < muPropagatedPositions_x->size(); iMu++)
    {
      if(muPropagatedRegion->at(iMu) == -1 && !(muChamber->at(iMu) == 20 && muLayer->at(iMu) == 1) && isME11->at(iMu) == 1 )
	{
	  for(std::size_t iRecHit = 0; iRecHit < recHitPositions_x->size(); iRecHit++)
	    {
	      if(recHitregion->at(iRecHit) == -1)
		{
		  if(recHitlayer->at(iRecHit)==1 && recHitchamber->at(iRecHit)==20) continue;
		  float residualx = std::fabs(muPropagatedPositions_x->at(iMu)-recHitPositions_x->at(iRecHit));
		  float residualy = std::fabs(muPropagatedPositions_y->at(iMu)-recHitPositions_y->at(iRecHit));
		  float residual = TMath::Sqrt(TMath::Power(residualx,2)+TMath::Power(residualy,2));
		  if(residual <= min_residual)
		    {
		      min_residual = residual;
		      min_recHitPosition_x = recHitPositions_x->at(iRecHit);
		      min_recHitPosition_y = recHitPositions_y->at(iRecHit);
		      //if(recHitchamber->(atiRecHit) == 8 || recHitchamber->(atiRecHit) == 9) continue;
		      mupt = muPt->at(iMu);
		      mueta = TMath::Abs(muEta->at(iMu));
		    }
		  else continue;
		}
	    }
	  
	  if(min_recHitPosition_x == 666 && min_recHitPosition_y == 666) continue;
	  else
	    {
	      m_plots["ResidualPlotGlobal_endcapneg"]->Fill(min_residual);
	      m_plots["OccupancyMatchedGlobal2D_endcapneg"]->Fill(min_recHitPosition_x,min_recHitPosition_y);
	      m_plots["MuPtMatchedGlobal_endcapneg"]->Fill(mupt);
	      m_plots["MuEtaMatchedGlobal_endcapneg"]->Fill(mueta);
	    }
	}
    
  
      if(muPropagatedRegion->at(iMu) == 1 && isME11->at(iMu) == 1)
	{
	  for(std::size_t iRecHit = 0; iRecHit < recHitPositions_x->size(); iRecHit++)
	    {
	      if(recHitregion->at(iRecHit) == 1)
		{
		  float residualx = std::fabs(muPropagatedPositions_x->at(iMu)-recHitPositions_x->at(iRecHit));
		  float residualy = std::fabs(muPropagatedPositions_y->at(iMu)-recHitPositions_y->at(iRecHit));
		  float residual = TMath::Sqrt(TMath::Power(residualx,2)+TMath::Power(residualy,2));
		  if(residual <= min_residual)
		    {
		      min_residual = residual;
		      min_recHitPosition_x = recHitPositions_x->at(iRecHit);
		      min_recHitPosition_y = recHitPositions_y->at(iRecHit);
		      mupt = muPt->at(iMu);
		      mueta = muEta->at(iMu);
		    }
		  else
		    {
		      continue;
		    }
		}
	    }
	  if(min_recHitPosition_x == 666 && min_recHitPosition_y == 666) continue;
	  else
	    {
	      m_plots["ResidualPlotGlobal_endcappos"]->Fill(min_residual);
	      m_plots["OccupancyMatchedGlobal2D_endcappos"]->Fill(min_recHitPosition_x,min_recHitPosition_y);
	      m_plots["MuPtMatchedGlobal_endcappos"]->Fill(mupt);
	      m_plots["MuEtaMatchedGlobal_endcappos"]->Fill(mueta);
	    }
	}
    }
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

void gemAnalysis::residual(vector<float> *recHitPositions, vector<float> *muPropagatedPositions, vector<int> *recHitregion, vector<int> *muPropagatedRegion, vector<int> *recHitlayer, vector<int> *recHitchamber)
{

float min_recHitPosition = 666; //DUMMY VALUE
vector<float> residuals;

for(std::size_t iMu = 0; iMu < muPropagatedPositions->size(); iMu++)
  {
    if(muPropagatedRegion->at(iMu) == -1)
      {
	float min_residual = 10;
	for(std::size_t iRecHit = 0; iRecHit < recHitPositions->size(); iRecHit++)
	  {
	    if(recHitregion->at(iRecHit) == -1)
	      {
		if(recHitlayer->at(iRecHit)==1 && recHitchamber->at(iRecHit)==20) continue;
		float residual = muPropagatedPositions->at(iMu)-recHitPositions->at(iRecHit);
		if(std::fabs(residual) <= min_residual)
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
	    m_plots["ResidualPlotLocX_endcapneg"]->Fill(min_residual);
	  }
      }
    
    if(muPropagatedRegion->at(iMu) == 1)
      {
        float min_residual = 10;
        for(std::size_t iRecHit = 0; iRecHit < recHitPositions->size(); iRecHit++)
          {
            if(recHitregion->at(iRecHit) == 1)
              {
                float residual = muPropagatedPositions->at(iMu)-recHitPositions->at(iRecHit);
                if(std::fabs(residual) <= min_residual)
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
            m_plots["ResidualPlotLocX_endcappos"]->Fill(min_residual);
          }
      }
  }
}


void gemAnalysis::book()
{
  m_outFile.cd();

       
  m_plots["ResidualPlotLocalX_endcapneg"] = new TH1F("ResidualPlotLocalX_endcapneg",
						   "Residual LocalX; residualX ; entries",
						   400,-2,2.);
  
  m_plots["ResidualPlotLocalX_endcappos"] = new TH1F("ResidualPlotLocalX_endcappos",
						   "Residual LocalX; residualX ; entries",
						   400,-2,2.);

  m_plots["ResidualPlotGlobal_endcapneg"] = new TH1F("ResidualPlotGlobal_endcapneg",
						     "Residual Global; residual; entries",
						     200,0,10.);
  
  m_plots["ResidualPlotGlobal_endcappos"] = new TH1F("ResidualPlotGlobal_endcappos",
                                                     "Residual Global; residual; entries",
                                                     200,0,10.);
    
  m_plots["OccupancyMatchedGlobal2D_endcapneg"] = new TH2F("OccupancyMatchedGlobal2D_endcapneg",
							   "Matched Global 2D; globalX ; globalY",
							   100,-250.,250.,
							   100,-250.,250.);
  
  m_plots["OccupancyMatchedGlobal2D_endcappos"] = new TH2F("OccupancyMatchedGlobal2D_endcappos",
							   "Matched Global 2D; globalX ; globalY",
							   100,-250.,250.,
							   100,-250.,250.);
  
  m_plots["OccupancyMatchedLocalX_endcapneg"] = new TH2F("OccupancyMatchedLocalX_endcapneg",
							 "Matched Local X; localX ; etaPartition",
							 100,-25.,25.,
							 8,1.,9.);
  
  m_plots["OccupancyMatchedLocalX_endcappos"] = new TH2F("OccupancyMatchedLocalX_endcappos",
							 "Matched Local X; localX ; etaPartition",
							 100,-25.,25.,
							 8,1.,9.);
  

  /*m_plots["OccupancyRecHit2D_layer1_endcappos"] = new TH2F("OccupancyRecHit2D_layer1_endcappos",
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
                                                         100,-250.,250.);*/
  
   
  /*m_plots["MatchedLocalEventNumber"] = new TH2F("MatchedLocalEventNumber",
					   "MatchedLocal vs EventNumber; eventNumber; # matched",
						30000000,0,30000000,
						100, 0., 100.);
 
  m_plots["LocalXEtaPartition"] = new TH2F("LocalXEtaPartition",
					   "LocalX vs EtaPartition; localX; etaPartition",
					   100,-20.,20.,
					   12,0.,12.);*/

  m_plots["MuPtMatchedLocal_endcapneg"] = new TH1F("MuPtMatchedLocal_endcapneg",
					      "MuPtMatched Local; p_{t}; entries",
					      200, 0., 100.);
    
  m_plots["MuPtMatchedLocal_endcappos"] = new TH1F("MuPtMatchedLocal_endcappos",
					     "MuPtMatched Local; p_{t}; entries",
					      200, 0., 100.);
  
  m_plots["MuEtaMatchedLocal_endcapneg"] = new TH1F("MuEtaMatchedLocal_endcapneg",
					       "MuEtaMatched Local; #eta; entries",
					       300, 0., 3.);

  m_plots["MuEtaMatchedLocal_endcappos"] = new TH1F("MuEtaMatchedLocal_endcappos",
                                               "MuEtaMatched Local; #eta; entries",
                                               300, 0., 3.);

  m_plots["MuPtMatchedGlobal_endcapneg"] = new TH1F("MuPtMatchedGlobal_endcapneg",
						   "MuPtMatched Global; p_{t}; entries",
						   200, 0., 100.);

  m_plots["MuPtMatchedGlobal_endcappos"] = new TH1F("MuPtMatchedGlobal_endcappos",
						   "MuPtMatched Global; p_{t}; entries",
						   200, 0., 100.);

  m_plots["MuEtaMatchedGlobal_endcapneg"] = new TH1F("MuEtaMatchedGlobal_endcapneg",
						    "MuEtaMatched Global; #eta; entries",
						    300, 0., 3.);

  m_plots["MuEtaMatchedGlobal_endcappos"] = new TH1F("MuEtaMatchedGlobal_endcappos",
						    "MuEtaMatched Global; #eta; entries",
						    600, -3., 3.);
  
}


void gemAnalysis::endJob()
{

  m_outFile.cd();
  m_outFile.Write();
  m_outFile.Close();

}

