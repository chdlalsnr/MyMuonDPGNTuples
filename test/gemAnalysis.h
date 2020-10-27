//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 13 11:45:56 2020 by ROOT version 6.22/02
// from TTree MuDPGTree/Mu DPG Tree
// found on file: MuDPGNtuple_11_1_2_patch2.root
//////////////////////////////////////////////////////////

#ifndef gemAnalysis_h
#define gemAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TString.h>
#include <TTree.h>
#include <iostream>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "TClonesArray.h"
#include "vector"
#include "vector"
#include "vector"



class gemAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   gemAnalysis(const TString & inFileName,
               const TString & outFileName);
   
   // Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types

   Int_t           event_runNumber;
   Int_t           event_lumiBlock;
   Long64_t        event_eventNumber;
   UInt_t          dtDigi_nDigis;
   vector<short>   *dtDigi_wheel;
   vector<short>   *dtDigi_sector;
   vector<short>   *dtDigi_station;
   vector<short>   *dtDigi_superLayer;
   vector<short>   *dtDigi_layer;
   vector<short>   *dtDigi_wire;
   vector<float>   *dtDigi_time;
   UInt_t          ph2DtDigi_nDigis;
   vector<short>   *ph2DtDigi_wheel;
   vector<short>   *ph2DtDigi_sector;
   vector<short>   *ph2DtDigi_station;
   vector<short>   *ph2DtDigi_superLayer;
   vector<short>   *ph2DtDigi_layer;
   vector<short>   *ph2DtDigi_wire;
   vector<float>   *ph2DtDigi_time;
   UInt_t          dtSeg_nSegments;
   vector<short>   *dtSeg_wheel;
   vector<short>   *dtSeg_sector;
   vector<short>   *dtSeg_station;
   vector<short>   *dtSeg_hasPhi;
   vector<short>   *dtSeg_hasZed;
   vector<float>   *dtSeg_posLoc_x;
   vector<float>   *dtSeg_posLoc_y;
   vector<float>   *dtSeg_posLoc_z;
   vector<float>   *dtSeg_dirLoc_x;
   vector<float>   *dtSeg_dirLoc_y;
   vector<float>   *dtSeg_dirLoc_z;
   vector<float>   *dtSeg_posLoc_x_SL1;
   vector<float>   *dtSeg_posLoc_x_SL3;
   vector<float>   *dtSeg_posLoc_x_midPlane;
   vector<float>   *dtSeg_posGlb_phi;
   vector<float>   *dtSeg_posGlb_eta;
   vector<float>   *dtSeg_dirGlb_phi;
   vector<float>   *dtSeg_dirGlb_eta;
   TClonesArray    *dtSeg_hitsExpPos;
   TClonesArray    *dtSeg_hitsExpPosCh;
   TClonesArray    *dtSeg_hitsExpWire;
   vector<float>   *dtSeg_phi_t0;
   vector<float>   *dtSeg_phi_vDrift;
   vector<float>   *dtSeg_phi_normChi2;
   vector<short>   *dtSeg_phi_nHits;
   TClonesArray    *dtSeg_phiHits_pos;
   TClonesArray    *dtSeg_phiHits_posCh;
   TClonesArray    *dtSeg_phiHits_posErr;
   TClonesArray    *dtSeg_phiHits_side;
   TClonesArray    *dtSeg_phiHits_wire;
   TClonesArray    *dtSeg_phiHits_wirePos;
   TClonesArray    *dtSeg_phiHits_layer;
   TClonesArray    *dtSeg_phiHits_superLayer;
   TClonesArray    *dtSeg_phiHits_time;
   TClonesArray    *dtSeg_phiHits_timeCali;
   vector<float>   *dtSeg_z_normChi2;
   vector<short>   *dtSeg_z_nHits;
   TClonesArray    *dtSeg_zHits_pos;
   TClonesArray    *dtSeg_zHits_posCh;
   TClonesArray    *dtSeg_zHits_posErr;
   TClonesArray    *dtSeg_zHits_side;
   TClonesArray    *dtSeg_zHits_wire;
   TClonesArray    *dtSeg_zHits_wirePos;
   TClonesArray    *dtSeg_zHits_layer;
   TClonesArray    *dtSeg_zHits_time;
   TClonesArray    *dtSeg_zHits_timeCali;
   UInt_t          ph2DtSeg_nSegments;
   vector<short>   *ph2DtSeg_wheel;
   vector<short>   *ph2DtSeg_sector;
   vector<short>   *ph2DtSeg_station;
   vector<short>   *ph2DtSeg_hasPhi;
   vector<short>   *ph2DtSeg_hasZed;
   vector<float>   *ph2DtSeg_posLoc_x;
   vector<float>   *ph2DtSeg_posLoc_y;
   vector<float>   *ph2DtSeg_posLoc_z;
   vector<float>   *ph2DtSeg_dirLoc_x;
   vector<float>   *ph2DtSeg_dirLoc_y;
   vector<float>   *ph2DtSeg_dirLoc_z;
   vector<float>   *ph2DtSeg_posLoc_x_SL1;
   vector<float>   *ph2DtSeg_posLoc_x_SL3;
   vector<float>   *ph2DtSeg_posLoc_x_midPlane;
   vector<float>   *ph2DtSeg_posGlb_phi;
   vector<float>   *ph2DtSeg_posGlb_eta;
   vector<float>   *ph2DtSeg_dirGlb_phi;
   vector<float>   *ph2DtSeg_dirGlb_eta;
   TClonesArray    *ph2DtSeg_hitsExpPos;
   TClonesArray    *ph2DtSeg_hitsExpPosCh;
   TClonesArray    *ph2DtSeg_hitsExpWire;
   vector<float>   *ph2DtSeg_phi_t0;
   vector<float>   *ph2DtSeg_phi_vDrift;
   vector<float>   *ph2DtSeg_phi_normChi2;
   vector<short>   *ph2DtSeg_phi_nHits;
   TClonesArray    *ph2DtSeg_phiHits_pos;
   TClonesArray    *ph2DtSeg_phiHits_posCh;
   TClonesArray    *ph2DtSeg_phiHits_posErr;
   TClonesArray    *ph2DtSeg_phiHits_side;
   TClonesArray    *ph2DtSeg_phiHits_wire;
   TClonesArray    *ph2DtSeg_phiHits_wirePos;
   TClonesArray    *ph2DtSeg_phiHits_layer;
   TClonesArray    *ph2DtSeg_phiHits_superLayer;
   TClonesArray    *ph2DtSeg_phiHits_time;
   TClonesArray    *ph2DtSeg_phiHits_timeCali;
   vector<float>   *ph2DtSeg_z_normChi2;
   vector<short>   *ph2DtSeg_z_nHits;
   TClonesArray    *ph2DtSeg_zHits_pos;
   TClonesArray    *ph2DtSeg_zHits_posCh;
   TClonesArray    *ph2DtSeg_zHits_posErr;
   TClonesArray    *ph2DtSeg_zHits_side;
   TClonesArray    *ph2DtSeg_zHits_wire;
   TClonesArray    *ph2DtSeg_zHits_wirePos;
   TClonesArray    *ph2DtSeg_zHits_layer;
   TClonesArray    *ph2DtSeg_zHits_time;
   TClonesArray    *ph2DtSeg_zHits_timeCali;
   UInt_t          gemDigi_nDigis;
   vector<short>   *gemDigi_station;
   vector<short>   *gemDigi_region;
   vector<short>   *gemDigi_roll;
   vector<short>   *gemDigi_bx;
   vector<short>   *gemDigi_strip;
   vector<short>   *gemDigi_pad;
   vector<float>   *gemDigi_g_r;
   vector<float>   *gemDigi_g_phi;
   vector<float>   *gemDigig_eta;
   vector<float>   *gemDigi_g_x;
   vector<float>   *gemDigi_g_y;
   vector<float>   *gemDigi_g_z;
   UInt_t          gemRecHit_nRecHits;
   vector<int>     *gemRecHit_cluster_size;
   vector<float>   *gemRecHit_g_r;
   vector<float>   *gemRecHit_g_phi;
   vector<float>   *gemRecHit_g_x;
   vector<float>   *gemRecHit_g_y;
   vector<float>   *gemRecHit_g_z;
   UInt_t          gemSegment_nSegments;
   vector<short>   *gemSegment_region;
   vector<short>   *gemSegment_ring;
   vector<short>   *gemSegment_station;
   vector<float>   *gemSegment_posLoc_x;
   vector<float>   *gemSegment_posLoc_y;
   vector<float>   *gemSegment_posLoc_z;
   vector<float>   *gemSegment_dirLoc_x;
   vector<float>   *gemSegment_dirLoc_y;
   vector<float>   *gemSegment_dirLoc_z;
   vector<float>   *gemSegment_posGlb_x;
   vector<float>   *gemSegment_posGlb_y;
   vector<float>   *gemSegment_posGlb_z;
   vector<float>   *gemSegment_time;
   vector<float>   *gemSegment_time_err;
   vector<double>  *gemSegment_chi2;
   vector<float>   *gemSegment_posGlb_phi;
   vector<float>   *gemSegment_posGlb_eta;
   vector<float>   *gemSegment_dirGlb_phi;
   vector<float>   *gemSegment_dirGlb_eta;
   UInt_t          mu_nMuons;
   vector<float>   *mu_pt;
   vector<float>   *mu_phi;
   vector<float>   *mu_eta;
   vector<short>   *mu_charge;
   vector<bool>    *mu_isGlobal;
   vector<bool>    *mu_isStandalone;
   vector<bool>    *mu_isTracker;
   vector<bool>    *mu_isGEM;
   vector<bool>    *mu_isLoose;
   vector<bool>    *mu_isMedium;
   vector<bool>    *mu_isTight;
   vector<float>   *mu_propagatedLoc_x;
   vector<float>   *mu_propagatedLoc_y;
   vector<float>   *mu_propagatedLoc_z;
   vector<float>   *mu_propagatedLoc_r;
   vector<float>   *mu_propagatedGlb_x;
   vector<float>   *mu_propagatedGlb_y;
   vector<float>   *mu_propagatedGlb_z;
   vector<float>   *mu_propagatedGlb_r;

   // List of branches
   TBranch        *b_dtDigi_nDigis;   //!
   TBranch        *b_dtDigi_wheel;   //!
   TBranch        *b_dtDigi_sector;   //!
   TBranch        *b_dtDigi_station;   //!
   TBranch        *b_dtDigi_superLayer;   //!
   TBranch        *b_dtDigi_layer;   //!
   TBranch        *b_dtDigi_wire;   //!
   TBranch        *b_dtDigi_time;   //!
   TBranch        *b_ph2DtDigi_nDigis;   //!
   TBranch        *b_ph2DtDigi_wheel;   //!
   TBranch        *b_ph2DtDigi_sector;   //!
   TBranch        *b_ph2DtDigi_station;   //!
   TBranch        *b_ph2DtDigi_superLayer;   //!
   TBranch        *b_ph2DtDigi_layer;   //!
   TBranch        *b_ph2DtDigi_wire;   //!
   TBranch        *b_ph2DtDigi_time;   //!
   TBranch        *b_dtSeg_nSegments;   //!
   TBranch        *b_dtSeg_wheel;   //!
   TBranch        *b_dtSeg_sector;   //!
   TBranch        *b_dtSeg_station;   //!
   TBranch        *b_dtSeg_hasPhi;   //!
   TBranch        *b_dtSeg_hasZed;   //!
   TBranch        *b_dtSeg_posLoc_x;   //!
   TBranch        *b_dtSeg_posLoc_y;   //!
   TBranch        *b_dtSeg_posLoc_z;   //!
   TBranch        *b_dtSeg_dirLoc_x;   //!
   TBranch        *b_dtSeg_dirLoc_y;   //!
   TBranch        *b_dtSeg_dirLoc_z;   //!
   TBranch        *b_dtSeg_posLoc_x_SL1;   //!
   TBranch        *b_dtSeg_posLoc_x_SL3;   //!
   TBranch        *b_dtSeg_posLoc_x_midPlane;   //!
   TBranch        *b_dtSeg_posGlb_phi;   //!
   TBranch        *b_dtSeg_posGlb_eta;   //!
   TBranch        *b_dtSeg_dirGlb_phi;   //!
   TBranch        *b_dtSeg_dirGlb_eta;   //!
   TBranch        *b_dtSeg_hitsExpPos;   //!
   TBranch        *b_dtSeg_hitsExpPosCh;   //!
   TBranch        *b_dtSeg_hitsExpWire;   //!
   TBranch        *b_dtSeg_phi_t0;   //!
   TBranch        *b_dtSeg_phi_vDrift;   //!
   TBranch        *b_dtSeg_phi_normChi2;   //!
   TBranch        *b_dtSeg_phi_nHits;   //!
   TBranch        *b_dtSeg_phiHits_pos;   //!
   TBranch        *b_dtSeg_phiHits_posCh;   //!
   TBranch        *b_dtSeg_phiHits_posErr;   //!
   TBranch        *b_dtSeg_phiHits_side;   //!
   TBranch        *b_dtSeg_phiHits_wire;   //!
   TBranch        *b_dtSeg_phiHits_wirePos;   //!
   TBranch        *b_dtSeg_phiHits_layer;   //!
   TBranch        *b_dtSeg_phiHits_superLayer;   //!
   TBranch        *b_dtSeg_phiHits_time;   //!
   TBranch        *b_dtSeg_phiHits_timeCali;   //!
   TBranch        *b_dtSeg_z_normChi2;   //!
   TBranch        *b_dtSeg_z_nHits;   //!
   TBranch        *b_dtSeg_zHits_pos;   //!
   TBranch        *b_dtSeg_zHits_posCh;   //!
   TBranch        *b_dtSeg_zHits_posErr;   //!
   TBranch        *b_dtSeg_zHits_side;   //!
   TBranch        *b_dtSeg_zHits_wire;   //!
   TBranch        *b_dtSeg_zHits_wirePos;   //!
   TBranch        *b_dtSeg_zHits_layer;   //!
   TBranch        *b_dtSeg_zHits_time;   //!
   TBranch        *b_dtSeg_zHits_timeCali;   //!
   TBranch        *b_ph2DtSeg_nSegments;   //!
   TBranch        *b_ph2DtSeg_wheel;   //!
   TBranch        *b_ph2DtSeg_sector;   //!
   TBranch        *b_ph2DtSeg_station;   //!
   TBranch        *b_ph2DtSeg_hasPhi;   //!
   TBranch        *b_ph2DtSeg_hasZed;   //!
   TBranch        *b_ph2DtSeg_posLoc_x;   //!
   TBranch        *b_ph2DtSeg_posLoc_y;   //!
   TBranch        *b_ph2DtSeg_posLoc_z;   //!
   TBranch        *b_ph2DtSeg_dirLoc_x;   //!
   TBranch        *b_ph2DtSeg_dirLoc_y;   //!
   TBranch        *b_ph2DtSeg_dirLoc_z;   //!
   TBranch        *b_ph2DtSeg_posLoc_x_SL1;   //!
   TBranch        *b_ph2DtSeg_posLoc_x_SL3;   //!
   TBranch        *b_ph2DtSeg_posLoc_x_midPlane;   //!
   TBranch        *b_ph2DtSeg_posGlb_phi;   //!
   TBranch        *b_ph2DtSeg_posGlb_eta;   //!
   TBranch        *b_ph2DtSeg_dirGlb_phi;   //!
   TBranch        *b_ph2DtSeg_dirGlb_eta;   //!
   TBranch        *b_ph2DtSeg_hitsExpPos;   //!
   TBranch        *b_ph2DtSeg_hitsExpPosCh;   //!
   TBranch        *b_ph2DtSeg_hitsExpWire;   //!
   TBranch        *b_ph2DtSeg_phi_t0;   //!
   TBranch        *b_ph2DtSeg_phi_vDrift;   //!
   TBranch        *b_ph2DtSeg_phi_normChi2;   //!
   TBranch        *b_ph2DtSeg_phi_nHits;   //!
   TBranch        *b_ph2DtSeg_phiHits_pos;   //!
   TBranch        *b_ph2DtSeg_phiHits_posCh;   //!
   TBranch        *b_ph2DtSeg_phiHits_posErr;   //!
   TBranch        *b_ph2DtSeg_phiHits_side;   //!
   TBranch        *b_ph2DtSeg_phiHits_wire;   //!
   TBranch        *b_ph2DtSeg_phiHits_wirePos;   //!
   TBranch        *b_ph2DtSeg_phiHits_layer;   //!
   TBranch        *b_ph2DtSeg_phiHits_superLayer;   //!
   TBranch        *b_ph2DtSeg_phiHits_time;   //!
   TBranch        *b_ph2DtSeg_phiHits_timeCali;   //!
   TBranch        *b_ph2DtSeg_z_normChi2;   //!
   TBranch        *b_ph2DtSeg_z_nHits;   //!
   TBranch        *b_ph2DtSeg_zHits_pos;   //!
   TBranch        *b_ph2DtSeg_zHits_posCh;   //!
   TBranch        *b_ph2DtSeg_zHits_posErr;   //!
   TBranch        *b_ph2DtSeg_zHits_side;   //!
   TBranch        *b_ph2DtSeg_zHits_wire;   //!
   TBranch        *b_ph2DtSeg_zHits_wirePos;   //!
   TBranch        *b_ph2DtSeg_zHits_layer;   //!
   TBranch        *b_ph2DtSeg_zHits_time;   //!
   TBranch        *b_ph2DtSeg_zHits_timeCali;   //!
   TBranch        *b_gemDigi_nDigis;   //!
   TBranch        *b_gemDigi_station;   //!
   TBranch        *b_gemDigi_region;   //!
   TBranch        *b_gemDigi_roll;   //!
   TBranch        *b_gemDigi_bx;   //!
   TBranch        *b_gemDigi_strip;   //!
   TBranch        *b_gemDigi_pad;   //!
   TBranch        *b_gemDigi_g_r;   //!
   TBranch        *b_gemDigi_g_phi;   //!
   TBranch        *b_gemDigig_eta;   //!
   TBranch        *b_gemDigi_g_x;   //!
   TBranch        *b_gemDigi_g_y;   //!
   TBranch        *b_gemDigi_g_z;   //!
   TBranch        *b_gemRecHit_nRecHits;   //!
   TBranch        *b_gemRecHit_cluster_size;   //!
   TBranch        *b_gemRecHit_g_r;   //!
   TBranch        *b_gemRecHit_g_phi;   //!
   TBranch        *b_gemRecHit_g_x;   //!
   TBranch        *b_gemRecHit_g_y;   //!
   TBranch        *b_gemRecHit_g_z;   //!
   TBranch        *b_gemSegment_nSegments;   //!
   TBranch        *b_gemSegment_region;   //!
   TBranch        *b_gemSegment_ring;   //!
   TBranch        *b_gemSegment_station;   //!
   TBranch        *b_gemSegment_posLoc_x;   //!
   TBranch        *b_gemSegment_posLoc_y;   //!
   TBranch        *b_gemSegment_posLoc_z;   //!
   TBranch        *b_gemSegment_dirLoc_x;   //!
   TBranch        *b_gemSegment_dirLoc_y;   //!
   TBranch        *b_gemSegment_dirLoc_z;   //!
   TBranch        *b_gemSegment_posGlb_x;   //!
   TBranch        *b_gemSegment_posGlb_y;   //!
   TBranch        *b_gemSegment_posGlb_z;   //!
   TBranch        *b_gemSegment_time;   //!
   TBranch        *b_gemSegment_time_err;   //!
   TBranch        *b_gemSegment_chi2;   //!
   TBranch        *b_gemSegment_posGlb_phi;   //!
   TBranch        *b_gemSegment_posGlb_eta;   //!
   TBranch        *b_gemSegment_dirGlb_phi;   //!
   TBranch        *b_gemSegment_dirGlb_eta;   //!
   TBranch        *b_mu_nMuons;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_charge;   //!
   TBranch        *b_mu_isGlobal;   //!
   TBranch        *b_mu_isStandalone;   //!
   TBranch        *b_mu_isTracker;   //!
   TBranch        *b_mu_isGEM;   //!
   TBranch        *b_mu_isLoose;   //!
   TBranch        *b_mu_isMedium;   //!
   TBranch        *b_mu_isTight;   //!
   TBranch        *b_mu_propagatedLoc_x;   //!
   TBranch        *b_mu_propagatedLoc_y;   //!
   TBranch        *b_mu_propagatedLoc_z;   //!
   TBranch        *b_mu_propagatedLoc_r;   //!
   TBranch        *b_mu_propagatedGlb_x;   //!
   TBranch        *b_mu_propagatedGlb_y;   //!
   TBranch        *b_mu_propagatedGlb_z;   //!
   TBranch        *b_mu_propagatedGlb_r;   //!

   gemAnalysis(TTree *tree=0);
   virtual ~gemAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

 protected:

   TFile m_inFile;
   TFile m_outFile;

   std::map<std::string, TH1*> m_plots;

   virtual vector<float> findMatchedHit(vector<float> *recHitPositions, size_t nRecHitPos, vector<float> *muPropagatedPositions, size_t nMuPropagatedPos);
   virtual vector<float> residual(vector<float> *recHitPositions, size_t nRecHitPos, vector<float> *muPropagatedPositions, size_t nMuPropagatedPos);
   virtual void book();
   virtual void fill_matched(vector<float> position, size_t nSize);
   virtual void fill_residual(vector<float> position, size_t nSize);
   virtual void endJob();

};

#endif

#ifdef gemAnalysis_cxx
gemAnalysis::gemAnalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("MuDPGNtuple_11_1_2_patch2.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("MuDPGNtuple_11_1_2_patch2.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("MuDPGNtuple_11_1_2_patch2.root:/muNtupleProducer");
      dir->GetObject("MuDPGTree",tree);

   }
   Init(tree);
}

gemAnalysis::~gemAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gemAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gemAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void gemAnalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   dtDigi_wheel = 0;
   dtDigi_sector = 0;
   dtDigi_station = 0;
   dtDigi_superLayer = 0;
   dtDigi_layer = 0;
   dtDigi_wire = 0;
   dtDigi_time = 0;
   ph2DtDigi_wheel = 0;
   ph2DtDigi_sector = 0;
   ph2DtDigi_station = 0;
   ph2DtDigi_superLayer = 0;
   ph2DtDigi_layer = 0;
   ph2DtDigi_wire = 0;
   ph2DtDigi_time = 0;
   dtSeg_wheel = 0;
   dtSeg_sector = 0;
   dtSeg_station = 0;
   dtSeg_hasPhi = 0;
   dtSeg_hasZed = 0;
   dtSeg_posLoc_x = 0;
   dtSeg_posLoc_y = 0;
   dtSeg_posLoc_z = 0;
   dtSeg_dirLoc_x = 0;
   dtSeg_dirLoc_y = 0;
   dtSeg_dirLoc_z = 0;
   dtSeg_posLoc_x_SL1 = 0;
   dtSeg_posLoc_x_SL3 = 0;
   dtSeg_posLoc_x_midPlane = 0;
   dtSeg_posGlb_phi = 0;
   dtSeg_posGlb_eta = 0;
   dtSeg_dirGlb_phi = 0;
   dtSeg_dirGlb_eta = 0;
   dtSeg_hitsExpPos = 0;
   dtSeg_hitsExpPosCh = 0;
   dtSeg_hitsExpWire = 0;
   dtSeg_phi_t0 = 0;
   dtSeg_phi_vDrift = 0;
   dtSeg_phi_normChi2 = 0;
   dtSeg_phi_nHits = 0;
   dtSeg_phiHits_pos = 0;
   dtSeg_phiHits_posCh = 0;
   dtSeg_phiHits_posErr = 0;
   dtSeg_phiHits_side = 0;
   dtSeg_phiHits_wire = 0;
   dtSeg_phiHits_wirePos = 0;
   dtSeg_phiHits_layer = 0;
   dtSeg_phiHits_superLayer = 0;
   dtSeg_phiHits_time = 0;
   dtSeg_phiHits_timeCali = 0;
   dtSeg_z_normChi2 = 0;
   dtSeg_z_nHits = 0;
   dtSeg_zHits_pos = 0;
   dtSeg_zHits_posCh = 0;
   dtSeg_zHits_posErr = 0;
   dtSeg_zHits_side = 0;
   dtSeg_zHits_wire = 0;
   dtSeg_zHits_wirePos = 0;
   dtSeg_zHits_layer = 0;
   dtSeg_zHits_time = 0;
   dtSeg_zHits_timeCali = 0;
   ph2DtSeg_wheel = 0;
   ph2DtSeg_sector = 0;
   ph2DtSeg_station = 0;
   ph2DtSeg_hasPhi = 0;
   ph2DtSeg_hasZed = 0;
   ph2DtSeg_posLoc_x = 0;
   ph2DtSeg_posLoc_y = 0;
   ph2DtSeg_posLoc_z = 0;
   ph2DtSeg_dirLoc_x = 0;
   ph2DtSeg_dirLoc_y = 0;
   ph2DtSeg_dirLoc_z = 0;
   ph2DtSeg_posLoc_x_SL1 = 0;
   ph2DtSeg_posLoc_x_SL3 = 0;
   ph2DtSeg_posLoc_x_midPlane = 0;
   ph2DtSeg_posGlb_phi = 0;
   ph2DtSeg_posGlb_eta = 0;
   ph2DtSeg_dirGlb_phi = 0;
   ph2DtSeg_dirGlb_eta = 0;
   ph2DtSeg_hitsExpPos = 0;
   ph2DtSeg_hitsExpPosCh = 0;
   ph2DtSeg_hitsExpWire = 0;
   ph2DtSeg_phi_t0 = 0;
   ph2DtSeg_phi_vDrift = 0;
   ph2DtSeg_phi_normChi2 = 0;
   ph2DtSeg_phi_nHits = 0;
   ph2DtSeg_phiHits_pos = 0;
   ph2DtSeg_phiHits_posCh = 0;
   ph2DtSeg_phiHits_posErr = 0;
   ph2DtSeg_phiHits_side = 0;
   ph2DtSeg_phiHits_wire = 0;
   ph2DtSeg_phiHits_wirePos = 0;
   ph2DtSeg_phiHits_layer = 0;
   ph2DtSeg_phiHits_superLayer = 0;
   ph2DtSeg_phiHits_time = 0;
   ph2DtSeg_phiHits_timeCali = 0;
   ph2DtSeg_z_normChi2 = 0;
   ph2DtSeg_z_nHits = 0;
   ph2DtSeg_zHits_pos = 0;
   ph2DtSeg_zHits_posCh = 0;
   ph2DtSeg_zHits_posErr = 0;
   ph2DtSeg_zHits_side = 0;
   ph2DtSeg_zHits_wire = 0;
   ph2DtSeg_zHits_wirePos = 0;
   ph2DtSeg_zHits_layer = 0;
   ph2DtSeg_zHits_time = 0;
   ph2DtSeg_zHits_timeCali = 0;
   gemDigi_station = 0;
   gemDigi_region = 0;
   gemDigi_roll = 0;
   gemDigi_bx = 0;
   gemDigi_strip = 0;
   gemDigi_pad = 0;
   gemDigi_g_r = 0;
   gemDigi_g_phi = 0;
   gemDigig_eta = 0;
   gemDigi_g_x = 0;
   gemDigi_g_y = 0;
   gemDigi_g_z = 0;
   gemRecHit_cluster_size = 0;
   gemRecHit_g_r = 0;
   gemRecHit_g_phi = 0;
   gemRecHit_g_x = 0;
   gemRecHit_g_y = 0;
   gemRecHit_g_z = 0;
   gemSegment_region = 0;
   gemSegment_ring = 0;
   gemSegment_station = 0;
   gemSegment_posLoc_x = 0;
   gemSegment_posLoc_y = 0;
   gemSegment_posLoc_z = 0;
   gemSegment_dirLoc_x = 0;
   gemSegment_dirLoc_y = 0;
   gemSegment_dirLoc_z = 0;
   gemSegment_posGlb_x = 0;
   gemSegment_posGlb_y = 0;
   gemSegment_posGlb_z = 0;
   gemSegment_time = 0;
   gemSegment_time_err = 0;
   gemSegment_chi2 = 0;
   gemSegment_posGlb_phi = 0;
   gemSegment_posGlb_eta = 0;
   gemSegment_dirGlb_phi = 0;
   gemSegment_dirGlb_eta = 0;
   mu_pt = 0;
   mu_phi = 0;
   mu_eta = 0;
   mu_charge = 0;
   mu_isGlobal = 0;
   mu_isStandalone = 0;
   mu_isTracker = 0;
   mu_isGEM = 0;
   mu_isLoose = 0;
   mu_isMedium = 0;
   mu_isTight = 0;
   mu_propagatedLoc_x = 0;
   mu_propagatedLoc_y = 0;
   mu_propagatedLoc_z = 0;
   mu_propagatedLoc_r = 0;
   mu_propagatedGlb_x = 0;
   mu_propagatedGlb_y = 0;
   mu_propagatedGlb_z = 0;
   mu_propagatedGlb_r = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("dtDigi_nDigis", &dtDigi_nDigis, &b_dtDigi_nDigis);
   fChain->SetBranchAddress("dtDigi_wheel", &dtDigi_wheel, &b_dtDigi_wheel);
   fChain->SetBranchAddress("dtDigi_sector", &dtDigi_sector, &b_dtDigi_sector);
   fChain->SetBranchAddress("dtDigi_station", &dtDigi_station, &b_dtDigi_station);
   fChain->SetBranchAddress("dtDigi_superLayer", &dtDigi_superLayer, &b_dtDigi_superLayer);
   fChain->SetBranchAddress("dtDigi_layer", &dtDigi_layer, &b_dtDigi_layer);
   fChain->SetBranchAddress("dtDigi_wire", &dtDigi_wire, &b_dtDigi_wire);
   fChain->SetBranchAddress("dtDigi_time", &dtDigi_time, &b_dtDigi_time);
   fChain->SetBranchAddress("ph2DtDigi_nDigis", &ph2DtDigi_nDigis, &b_ph2DtDigi_nDigis);
   fChain->SetBranchAddress("ph2DtDigi_wheel", &ph2DtDigi_wheel, &b_ph2DtDigi_wheel);
   fChain->SetBranchAddress("ph2DtDigi_sector", &ph2DtDigi_sector, &b_ph2DtDigi_sector);
   fChain->SetBranchAddress("ph2DtDigi_station", &ph2DtDigi_station, &b_ph2DtDigi_station);
   fChain->SetBranchAddress("ph2DtDigi_superLayer", &ph2DtDigi_superLayer, &b_ph2DtDigi_superLayer);
   fChain->SetBranchAddress("ph2DtDigi_layer", &ph2DtDigi_layer, &b_ph2DtDigi_layer);
   fChain->SetBranchAddress("ph2DtDigi_wire", &ph2DtDigi_wire, &b_ph2DtDigi_wire);
   fChain->SetBranchAddress("ph2DtDigi_time", &ph2DtDigi_time, &b_ph2DtDigi_time);
   fChain->SetBranchAddress("dtSeg_nSegments", &dtSeg_nSegments, &b_dtSeg_nSegments);
   fChain->SetBranchAddress("dtSeg_wheel", &dtSeg_wheel, &b_dtSeg_wheel);
   fChain->SetBranchAddress("dtSeg_sector", &dtSeg_sector, &b_dtSeg_sector);
   fChain->SetBranchAddress("dtSeg_station", &dtSeg_station, &b_dtSeg_station);
   fChain->SetBranchAddress("dtSeg_hasPhi", &dtSeg_hasPhi, &b_dtSeg_hasPhi);
   fChain->SetBranchAddress("dtSeg_hasZed", &dtSeg_hasZed, &b_dtSeg_hasZed);
   fChain->SetBranchAddress("dtSeg_posLoc_x", &dtSeg_posLoc_x, &b_dtSeg_posLoc_x);
   fChain->SetBranchAddress("dtSeg_posLoc_y", &dtSeg_posLoc_y, &b_dtSeg_posLoc_y);
   fChain->SetBranchAddress("dtSeg_posLoc_z", &dtSeg_posLoc_z, &b_dtSeg_posLoc_z);
   fChain->SetBranchAddress("dtSeg_dirLoc_x", &dtSeg_dirLoc_x, &b_dtSeg_dirLoc_x);
   fChain->SetBranchAddress("dtSeg_dirLoc_y", &dtSeg_dirLoc_y, &b_dtSeg_dirLoc_y);
   fChain->SetBranchAddress("dtSeg_dirLoc_z", &dtSeg_dirLoc_z, &b_dtSeg_dirLoc_z);
   fChain->SetBranchAddress("dtSeg_posLoc_x_SL1", &dtSeg_posLoc_x_SL1, &b_dtSeg_posLoc_x_SL1);
   fChain->SetBranchAddress("dtSeg_posLoc_x_SL3", &dtSeg_posLoc_x_SL3, &b_dtSeg_posLoc_x_SL3);
   fChain->SetBranchAddress("dtSeg_posLoc_x_midPlane", &dtSeg_posLoc_x_midPlane, &b_dtSeg_posLoc_x_midPlane);
   fChain->SetBranchAddress("dtSeg_posGlb_phi", &dtSeg_posGlb_phi, &b_dtSeg_posGlb_phi);
   fChain->SetBranchAddress("dtSeg_posGlb_eta", &dtSeg_posGlb_eta, &b_dtSeg_posGlb_eta);
   fChain->SetBranchAddress("dtSeg_dirGlb_phi", &dtSeg_dirGlb_phi, &b_dtSeg_dirGlb_phi);
   fChain->SetBranchAddress("dtSeg_dirGlb_eta", &dtSeg_dirGlb_eta, &b_dtSeg_dirGlb_eta);
   fChain->SetBranchAddress("dtSeg_hitsExpPos", &dtSeg_hitsExpPos, &b_dtSeg_hitsExpPos);
   fChain->SetBranchAddress("dtSeg_hitsExpPosCh", &dtSeg_hitsExpPosCh, &b_dtSeg_hitsExpPosCh);
   fChain->SetBranchAddress("dtSeg_hitsExpWire", &dtSeg_hitsExpWire, &b_dtSeg_hitsExpWire);
   fChain->SetBranchAddress("dtSeg_phi_t0", &dtSeg_phi_t0, &b_dtSeg_phi_t0);
   fChain->SetBranchAddress("dtSeg_phi_vDrift", &dtSeg_phi_vDrift, &b_dtSeg_phi_vDrift);
   fChain->SetBranchAddress("dtSeg_phi_normChi2", &dtSeg_phi_normChi2, &b_dtSeg_phi_normChi2);
   fChain->SetBranchAddress("dtSeg_phi_nHits", &dtSeg_phi_nHits, &b_dtSeg_phi_nHits);
   fChain->SetBranchAddress("dtSeg_phiHits_pos", &dtSeg_phiHits_pos, &b_dtSeg_phiHits_pos);
   fChain->SetBranchAddress("dtSeg_phiHits_posCh", &dtSeg_phiHits_posCh, &b_dtSeg_phiHits_posCh);
   fChain->SetBranchAddress("dtSeg_phiHits_posErr", &dtSeg_phiHits_posErr, &b_dtSeg_phiHits_posErr);
   fChain->SetBranchAddress("dtSeg_phiHits_side", &dtSeg_phiHits_side, &b_dtSeg_phiHits_side);
   fChain->SetBranchAddress("dtSeg_phiHits_wire", &dtSeg_phiHits_wire, &b_dtSeg_phiHits_wire);
   fChain->SetBranchAddress("dtSeg_phiHits_wirePos", &dtSeg_phiHits_wirePos, &b_dtSeg_phiHits_wirePos);
   fChain->SetBranchAddress("dtSeg_phiHits_layer", &dtSeg_phiHits_layer, &b_dtSeg_phiHits_layer);
   fChain->SetBranchAddress("dtSeg_phiHits_superLayer", &dtSeg_phiHits_superLayer, &b_dtSeg_phiHits_superLayer);
   fChain->SetBranchAddress("dtSeg_phiHits_time", &dtSeg_phiHits_time, &b_dtSeg_phiHits_time);
   fChain->SetBranchAddress("dtSeg_phiHits_timeCali", &dtSeg_phiHits_timeCali, &b_dtSeg_phiHits_timeCali);
   fChain->SetBranchAddress("dtSeg_z_normChi2", &dtSeg_z_normChi2, &b_dtSeg_z_normChi2);
   fChain->SetBranchAddress("dtSeg_z_nHits", &dtSeg_z_nHits, &b_dtSeg_z_nHits);
   fChain->SetBranchAddress("dtSeg_zHits_pos", &dtSeg_zHits_pos, &b_dtSeg_zHits_pos);
   fChain->SetBranchAddress("dtSeg_zHits_posCh", &dtSeg_zHits_posCh, &b_dtSeg_zHits_posCh);
   fChain->SetBranchAddress("dtSeg_zHits_posErr", &dtSeg_zHits_posErr, &b_dtSeg_zHits_posErr);
   fChain->SetBranchAddress("dtSeg_zHits_side", &dtSeg_zHits_side, &b_dtSeg_zHits_side);
   fChain->SetBranchAddress("dtSeg_zHits_wire", &dtSeg_zHits_wire, &b_dtSeg_zHits_wire);
   fChain->SetBranchAddress("dtSeg_zHits_wirePos", &dtSeg_zHits_wirePos, &b_dtSeg_zHits_wirePos);
   fChain->SetBranchAddress("dtSeg_zHits_layer", &dtSeg_zHits_layer, &b_dtSeg_zHits_layer);
   fChain->SetBranchAddress("dtSeg_zHits_time", &dtSeg_zHits_time, &b_dtSeg_zHits_time);
   fChain->SetBranchAddress("dtSeg_zHits_timeCali", &dtSeg_zHits_timeCali, &b_dtSeg_zHits_timeCali);
   fChain->SetBranchAddress("ph2DtSeg_nSegments", &ph2DtSeg_nSegments, &b_ph2DtSeg_nSegments);
   fChain->SetBranchAddress("ph2DtSeg_wheel", &ph2DtSeg_wheel, &b_ph2DtSeg_wheel);
   fChain->SetBranchAddress("ph2DtSeg_sector", &ph2DtSeg_sector, &b_ph2DtSeg_sector);
   fChain->SetBranchAddress("ph2DtSeg_station", &ph2DtSeg_station, &b_ph2DtSeg_station);
   fChain->SetBranchAddress("ph2DtSeg_hasPhi", &ph2DtSeg_hasPhi, &b_ph2DtSeg_hasPhi);
   fChain->SetBranchAddress("ph2DtSeg_hasZed", &ph2DtSeg_hasZed, &b_ph2DtSeg_hasZed);
   fChain->SetBranchAddress("ph2DtSeg_posLoc_x", &ph2DtSeg_posLoc_x, &b_ph2DtSeg_posLoc_x);
   fChain->SetBranchAddress("ph2DtSeg_posLoc_y", &ph2DtSeg_posLoc_y, &b_ph2DtSeg_posLoc_y);
   fChain->SetBranchAddress("ph2DtSeg_posLoc_z", &ph2DtSeg_posLoc_z, &b_ph2DtSeg_posLoc_z);
   fChain->SetBranchAddress("ph2DtSeg_dirLoc_x", &ph2DtSeg_dirLoc_x, &b_ph2DtSeg_dirLoc_x);
   fChain->SetBranchAddress("ph2DtSeg_dirLoc_y", &ph2DtSeg_dirLoc_y, &b_ph2DtSeg_dirLoc_y);
   fChain->SetBranchAddress("ph2DtSeg_dirLoc_z", &ph2DtSeg_dirLoc_z, &b_ph2DtSeg_dirLoc_z);
   fChain->SetBranchAddress("ph2DtSeg_posLoc_x_SL1", &ph2DtSeg_posLoc_x_SL1, &b_ph2DtSeg_posLoc_x_SL1);
   fChain->SetBranchAddress("ph2DtSeg_posLoc_x_SL3", &ph2DtSeg_posLoc_x_SL3, &b_ph2DtSeg_posLoc_x_SL3);
   fChain->SetBranchAddress("ph2DtSeg_posLoc_x_midPlane", &ph2DtSeg_posLoc_x_midPlane, &b_ph2DtSeg_posLoc_x_midPlane);
   fChain->SetBranchAddress("ph2DtSeg_posGlb_phi", &ph2DtSeg_posGlb_phi, &b_ph2DtSeg_posGlb_phi);
   fChain->SetBranchAddress("ph2DtSeg_posGlb_eta", &ph2DtSeg_posGlb_eta, &b_ph2DtSeg_posGlb_eta);
   fChain->SetBranchAddress("ph2DtSeg_dirGlb_phi", &ph2DtSeg_dirGlb_phi, &b_ph2DtSeg_dirGlb_phi);
   fChain->SetBranchAddress("ph2DtSeg_dirGlb_eta", &ph2DtSeg_dirGlb_eta, &b_ph2DtSeg_dirGlb_eta);
   fChain->SetBranchAddress("ph2DtSeg_hitsExpPos", &ph2DtSeg_hitsExpPos, &b_ph2DtSeg_hitsExpPos);
   fChain->SetBranchAddress("ph2DtSeg_hitsExpPosCh", &ph2DtSeg_hitsExpPosCh, &b_ph2DtSeg_hitsExpPosCh);
   fChain->SetBranchAddress("ph2DtSeg_hitsExpWire", &ph2DtSeg_hitsExpWire, &b_ph2DtSeg_hitsExpWire);
   fChain->SetBranchAddress("ph2DtSeg_phi_t0", &ph2DtSeg_phi_t0, &b_ph2DtSeg_phi_t0);
   fChain->SetBranchAddress("ph2DtSeg_phi_vDrift", &ph2DtSeg_phi_vDrift, &b_ph2DtSeg_phi_vDrift);
   fChain->SetBranchAddress("ph2DtSeg_phi_normChi2", &ph2DtSeg_phi_normChi2, &b_ph2DtSeg_phi_normChi2);
   fChain->SetBranchAddress("ph2DtSeg_phi_nHits", &ph2DtSeg_phi_nHits, &b_ph2DtSeg_phi_nHits);
   fChain->SetBranchAddress("ph2DtSeg_phiHits_pos", &ph2DtSeg_phiHits_pos, &b_ph2DtSeg_phiHits_pos);
   fChain->SetBranchAddress("ph2DtSeg_phiHits_posCh", &ph2DtSeg_phiHits_posCh, &b_ph2DtSeg_phiHits_posCh);
   fChain->SetBranchAddress("ph2DtSeg_phiHits_posErr", &ph2DtSeg_phiHits_posErr, &b_ph2DtSeg_phiHits_posErr);
   fChain->SetBranchAddress("ph2DtSeg_phiHits_side", &ph2DtSeg_phiHits_side, &b_ph2DtSeg_phiHits_side);
   fChain->SetBranchAddress("ph2DtSeg_phiHits_wire", &ph2DtSeg_phiHits_wire, &b_ph2DtSeg_phiHits_wire);
   fChain->SetBranchAddress("ph2DtSeg_phiHits_wirePos", &ph2DtSeg_phiHits_wirePos, &b_ph2DtSeg_phiHits_wirePos);
   fChain->SetBranchAddress("ph2DtSeg_phiHits_layer", &ph2DtSeg_phiHits_layer, &b_ph2DtSeg_phiHits_layer);
   fChain->SetBranchAddress("ph2DtSeg_phiHits_superLayer", &ph2DtSeg_phiHits_superLayer, &b_ph2DtSeg_phiHits_superLayer);
   fChain->SetBranchAddress("ph2DtSeg_phiHits_time", &ph2DtSeg_phiHits_time, &b_ph2DtSeg_phiHits_time);
   fChain->SetBranchAddress("ph2DtSeg_phiHits_timeCali", &ph2DtSeg_phiHits_timeCali, &b_ph2DtSeg_phiHits_timeCali);
   fChain->SetBranchAddress("ph2DtSeg_z_normChi2", &ph2DtSeg_z_normChi2, &b_ph2DtSeg_z_normChi2);
   fChain->SetBranchAddress("ph2DtSeg_z_nHits", &ph2DtSeg_z_nHits, &b_ph2DtSeg_z_nHits);
   fChain->SetBranchAddress("ph2DtSeg_zHits_pos", &ph2DtSeg_zHits_pos, &b_ph2DtSeg_zHits_pos);
   fChain->SetBranchAddress("ph2DtSeg_zHits_posCh", &ph2DtSeg_zHits_posCh, &b_ph2DtSeg_zHits_posCh);
   fChain->SetBranchAddress("ph2DtSeg_zHits_posErr", &ph2DtSeg_zHits_posErr, &b_ph2DtSeg_zHits_posErr);
   fChain->SetBranchAddress("ph2DtSeg_zHits_side", &ph2DtSeg_zHits_side, &b_ph2DtSeg_zHits_side);
   fChain->SetBranchAddress("ph2DtSeg_zHits_wire", &ph2DtSeg_zHits_wire, &b_ph2DtSeg_zHits_wire);
   fChain->SetBranchAddress("ph2DtSeg_zHits_wirePos", &ph2DtSeg_zHits_wirePos, &b_ph2DtSeg_zHits_wirePos);
   fChain->SetBranchAddress("ph2DtSeg_zHits_layer", &ph2DtSeg_zHits_layer, &b_ph2DtSeg_zHits_layer);
   fChain->SetBranchAddress("ph2DtSeg_zHits_time", &ph2DtSeg_zHits_time, &b_ph2DtSeg_zHits_time);
   fChain->SetBranchAddress("ph2DtSeg_zHits_timeCali", &ph2DtSeg_zHits_timeCali, &b_ph2DtSeg_zHits_timeCali);
   fChain->SetBranchAddress("gemDigi_nDigis", &gemDigi_nDigis, &b_gemDigi_nDigis);
   fChain->SetBranchAddress("gemDigi_station", &gemDigi_station, &b_gemDigi_station);
   fChain->SetBranchAddress("gemDigi_region", &gemDigi_region, &b_gemDigi_region);
   fChain->SetBranchAddress("gemDigi_roll", &gemDigi_roll, &b_gemDigi_roll);
   fChain->SetBranchAddress("gemDigi_bx", &gemDigi_bx, &b_gemDigi_bx);
   fChain->SetBranchAddress("gemDigi_strip", &gemDigi_strip, &b_gemDigi_strip);
   fChain->SetBranchAddress("gemDigi_pad", &gemDigi_pad, &b_gemDigi_pad);
   fChain->SetBranchAddress("gemDigi_g_r", &gemDigi_g_r, &b_gemDigi_g_r);
   fChain->SetBranchAddress("gemDigi_g_phi", &gemDigi_g_phi, &b_gemDigi_g_phi);
   fChain->SetBranchAddress("gemDigig_eta", &gemDigig_eta, &b_gemDigig_eta);
   fChain->SetBranchAddress("gemDigi_g_x", &gemDigi_g_x, &b_gemDigi_g_x);
   fChain->SetBranchAddress("gemDigi_g_y", &gemDigi_g_y, &b_gemDigi_g_y);
   fChain->SetBranchAddress("gemDigi_g_z", &gemDigi_g_z, &b_gemDigi_g_z);
   fChain->SetBranchAddress("gemRecHit_nRecHits", &gemRecHit_nRecHits, &b_gemRecHit_nRecHits);
   fChain->SetBranchAddress("gemRecHit_cluster_size", &gemRecHit_cluster_size, &b_gemRecHit_cluster_size);
   fChain->SetBranchAddress("gemRecHit_g_r", &gemRecHit_g_r, &b_gemRecHit_g_r);
   fChain->SetBranchAddress("gemRecHit_g_phi", &gemRecHit_g_phi, &b_gemRecHit_g_phi);
   fChain->SetBranchAddress("gemRecHit_g_x", &gemRecHit_g_x, &b_gemRecHit_g_x);
   fChain->SetBranchAddress("gemRecHit_g_y", &gemRecHit_g_y, &b_gemRecHit_g_y);
   fChain->SetBranchAddress("gemRecHit_g_z", &gemRecHit_g_z, &b_gemRecHit_g_z);
   fChain->SetBranchAddress("gemSegment_nSegments", &gemSegment_nSegments, &b_gemSegment_nSegments);
   fChain->SetBranchAddress("gemSegment_region", &gemSegment_region, &b_gemSegment_region);
   fChain->SetBranchAddress("gemSegment_ring", &gemSegment_ring, &b_gemSegment_ring);
   fChain->SetBranchAddress("gemSegment_station", &gemSegment_station, &b_gemSegment_station);
   fChain->SetBranchAddress("gemSegment_posLoc_x", &gemSegment_posLoc_x, &b_gemSegment_posLoc_x);
   fChain->SetBranchAddress("gemSegment_posLoc_y", &gemSegment_posLoc_y, &b_gemSegment_posLoc_y);
   fChain->SetBranchAddress("gemSegment_posLoc_z", &gemSegment_posLoc_z, &b_gemSegment_posLoc_z);
   fChain->SetBranchAddress("gemSegment_dirLoc_x", &gemSegment_dirLoc_x, &b_gemSegment_dirLoc_x);
   fChain->SetBranchAddress("gemSegment_dirLoc_y", &gemSegment_dirLoc_y, &b_gemSegment_dirLoc_y);
   fChain->SetBranchAddress("gemSegment_dirLoc_z", &gemSegment_dirLoc_z, &b_gemSegment_dirLoc_z);
   fChain->SetBranchAddress("gemSegment_posGlb_x", &gemSegment_posGlb_x, &b_gemSegment_posGlb_x);
   fChain->SetBranchAddress("gemSegment_posGlb_y", &gemSegment_posGlb_y, &b_gemSegment_posGlb_y);
   fChain->SetBranchAddress("gemSegment_posGlb_z", &gemSegment_posGlb_z, &b_gemSegment_posGlb_z);
   fChain->SetBranchAddress("gemSegment_time", &gemSegment_time, &b_gemSegment_time);
   fChain->SetBranchAddress("gemSegment_time_err", &gemSegment_time_err, &b_gemSegment_time_err);
   fChain->SetBranchAddress("gemSegment_chi2", &gemSegment_chi2, &b_gemSegment_chi2);
   fChain->SetBranchAddress("gemSegment_posGlb_phi", &gemSegment_posGlb_phi, &b_gemSegment_posGlb_phi);
   fChain->SetBranchAddress("gemSegment_posGlb_eta", &gemSegment_posGlb_eta, &b_gemSegment_posGlb_eta);
   fChain->SetBranchAddress("gemSegment_dirGlb_phi", &gemSegment_dirGlb_phi, &b_gemSegment_dirGlb_phi);
   fChain->SetBranchAddress("gemSegment_dirGlb_eta", &gemSegment_dirGlb_eta, &b_gemSegment_dirGlb_eta);
   fChain->SetBranchAddress("mu_nMuons", &mu_nMuons, &b_mu_nMuons);
   fChain->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
   fChain->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
   fChain->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
   fChain->SetBranchAddress("mu_charge", &mu_charge, &b_mu_charge);
   fChain->SetBranchAddress("mu_isGlobal", &mu_isGlobal, &b_mu_isGlobal);
   fChain->SetBranchAddress("mu_isStandalone", &mu_isStandalone, &b_mu_isStandalone);
   fChain->SetBranchAddress("mu_isTracker", &mu_isTracker, &b_mu_isTracker);
   fChain->SetBranchAddress("mu_isGEM", &mu_isGEM, &b_mu_isGEM);
   fChain->SetBranchAddress("mu_isLoose", &mu_isLoose, &b_mu_isLoose);
   fChain->SetBranchAddress("mu_isMedium", &mu_isMedium, &b_mu_isMedium);
   fChain->SetBranchAddress("mu_isTight", &mu_isTight, &b_mu_isTight);
   fChain->SetBranchAddress("mu_propagatedLoc_x", &mu_propagatedLoc_x, &b_mu_propagatedLoc_x);
   fChain->SetBranchAddress("mu_propagatedLoc_y", &mu_propagatedLoc_y, &b_mu_propagatedLoc_y);
   fChain->SetBranchAddress("mu_propagatedLoc_z", &mu_propagatedLoc_z, &b_mu_propagatedLoc_z);
   fChain->SetBranchAddress("mu_propagatedLoc_r", &mu_propagatedLoc_r, &b_mu_propagatedLoc_r);
   fChain->SetBranchAddress("mu_propagatedGlb_x", &mu_propagatedGlb_x, &b_mu_propagatedGlb_x);
   fChain->SetBranchAddress("mu_propagatedGlb_y", &mu_propagatedGlb_y, &b_mu_propagatedGlb_y);
   fChain->SetBranchAddress("mu_propagatedGlb_z", &mu_propagatedGlb_z, &b_mu_propagatedGlb_z);
   fChain->SetBranchAddress("mu_propagatedGlb_r", &mu_propagatedGlb_r, &b_mu_propagatedGlb_r);
   Notify();
}

Bool_t gemAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gemAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gemAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gemAnalysis_cxx
