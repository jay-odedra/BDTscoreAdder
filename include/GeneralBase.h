//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 24 00:10:28 2023 by ROOT version 6.26/11
// from TChain tree/
//////////////////////////////////////////////////////////

#ifndef GeneralBase_h
#define GeneralBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class GeneralBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           luminosityBlock;
   Long64_t        event;
   Int_t           bunchCrossing;
   Int_t           nBToKEE;
   Int_t           nElectron;
   Int_t           nProbeTracks;
   Bool_t          HLT_DoubleEle10_eta1p22_mMax6;
   Bool_t          HLT_DoubleEle9p5_eta1p22_mMax6;
   Bool_t          HLT_DoubleEle9_eta1p22_mMax6;
   Bool_t          HLT_DoubleEle8p5_eta1p22_mMax6;
   Bool_t          HLT_DoubleEle8_eta1p22_mMax6;
   Bool_t          HLT_DoubleEle7p5_eta1p22_mMax6;
   Bool_t          HLT_DoubleEle7_eta1p22_mMax6;
   Bool_t          HLT_DoubleEle6p5_eta1p22_mMax6;
   Bool_t          HLT_DoubleEle6_eta1p22_mMax6;
   Bool_t          HLT_DoubleEle5p5_eta1p22_mMax6;
   Bool_t          HLT_DoubleEle5_eta1p22_mMax6;
   Bool_t          HLT_DoubleEle4p5_eta1p22_mMax6;
   Bool_t          HLT_DoubleEle4_eta1p22_mMax6;
   Bool_t          L1_DoubleEG11_er1p2_dR_Max0p6;
   UChar_t         L1_DoubleEG10p5_er1p2_dR_Max0p6;
   Bool_t          L1_DoubleEG10_er1p2_dR_Max0p6;
   Bool_t          L1_DoubleEG9p5_er1p2_dR_Max0p6;
   Bool_t          L1_DoubleEG9_er1p2_dR_Max0p7;
   Bool_t          L1_DoubleEG8p5_er1p2_dR_Max0p7;
   Bool_t          L1_DoubleEG8_er1p2_dR_Max0p7;
   Bool_t          L1_DoubleEG7p5_er1p2_dR_Max0p7;
   Bool_t          L1_DoubleEG7_er1p2_dR_Max0p8;
   Bool_t          L1_DoubleEG6p5_er1p2_dR_Max0p8;
   Bool_t          L1_DoubleEG6_er1p2_dR_Max0p8;
   Bool_t          L1_DoubleEG5p5_er1p2_dR_Max0p8;
   Bool_t          L1_DoubleEG5_er1p2_dR_Max0p9;
   Bool_t          L1_DoubleEG4p5_er1p2_dR_Max0p9;
   Bool_t          L1_DoubleEG4_er1p2_dR_Max0p9;
   Bool_t          L1_4p5_HLT_4p0;
   Bool_t          L1_5p0_HLT_4p0;
   Bool_t          L1_5p5_HLT_4p0;
   Bool_t          L1_5p5_HLT_6p0;
   Bool_t          L1_6p0_HLT_4p0;
   Bool_t          L1_6p5_HLT_4p5;
   Bool_t          L1_7p0_HLT_5p0;
   Bool_t          L1_7p5_HLT_5p0;
   Bool_t          L1_8p0_HLT_5p0;
   Bool_t          L1_8p5_HLT_5p0;
   Bool_t          L1_8p5_HLT_5p5;
   Bool_t          L1_9p0_HLT_6p0;
   Bool_t          L1_10p5_HLT_5p0;
   Bool_t          L1_10p5_HLT_6p5;
   Bool_t          L1_11p0_HLT_6p5;
   Bool_t          trigger_OR;
   Int_t           nTrigObj;
   Int_t           nOtherPV;
   Float_t         PV_ndof;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_chi2;
   Float_t         PV_score;
   Int_t           PV_npvs;
   Int_t           PV_npvsGood;
   Int_t           nSV;
   Int_t           nSkimBToKEE;
   Float_t         BToKEE_fit_pt;
   Float_t         BToKEE_fit_eta;
   Float_t         BToKEE_fit_phi;
   Float_t         BToKEE_D0_mass_LepToK_KToPi;
   Float_t         BToKEE_D0_mass_LepToPi_KToK;
   Float_t         BToKEE_fit_mass;
   Float_t         BToKEE_l_xy;
   Float_t         BToKEE_l_xy_unc;
   Float_t         BToKEE_fit_cos2D;
   Float_t         BToKEE_svprob;
   Float_t         BToKEE_fit_massErr;
   Float_t         BToKEE_b_iso04;
   Float_t         BToKEE_mll_fullfit;
   Float_t         BToKEE_vtx_x;
   Float_t         BToKEE_vtx_y;
   Float_t         BToKEE_vtx_z;
   Float_t         BToKEE_l1Idx;
   Float_t         BToKEE_l2Idx;
   Float_t         BToKEE_kIdx;
   Float_t         BToKEE_fit_k_pt;
   Float_t         BToKEE_fit_k_eta;
   Float_t         BToKEE_fit_k_phi;
   Float_t         BToKEE_fit_l1_pt;
   Float_t         BToKEE_fit_l1_eta;
   Float_t         BToKEE_fit_l1_phi;
   Float_t         BToKEE_fit_l2_pt;
   Float_t         BToKEE_fit_l2_eta;
   Float_t         BToKEE_fit_l2_phi;
   Float_t         BToKEE_l1_iso04;
   Float_t         BToKEE_l2_iso04;
   Float_t         BToKEE_l1_isPF;
   Float_t         BToKEE_l2_isPF;
   Float_t         BToKEE_k_iso04;
   Float_t         BToKEE_l1_isPFoverlap;
   Float_t         BToKEE_l2_isPFoverlap;
   Float_t         BToKEE_l1_LooseId;
   Float_t         BToKEE_l2_LooseId;
   Float_t         BToKEE_l1_MediumId;
   Float_t         BToKEE_l2_MediumId;
   Float_t         BToKEE_l1_TightId;
   Float_t         BToKEE_l2_TightId;
   Float_t         BToKEE_l1_ConvVeto;
   Float_t         BToKEE_l2_ConvVeto;
   Float_t         BToKEE_l1_PFMvaID_Fall17;
   Float_t         BToKEE_l2_PFMvaID_Fall17;
   Float_t         BToKEE_l1_PFMvaID_retrained;
   Float_t         BToKEE_l2_PFMvaID_retrained;
   Float_t         BToKEE_l1_iso04_dca;
   Float_t         BToKEE_l2_iso04_dca;
   Float_t         BToKEE_b_iso04_dca;
   Float_t         BToKEE_k_iso04_dca;
   Float_t         BToKEE_k_svip3d;
   Float_t         BToKEE_k_svip3d_err;
   Float_t         BToKEE_l1_n_isotrk_dca;
   Float_t         BToKEE_l2_n_isotrk_dca;
   Float_t         BToKEE_k_n_isotrk_dca;
   Float_t         BToKEE_l1_vx;
   Float_t         BToKEE_l1_vy;
   Float_t         BToKEE_l1_vz;
   Float_t         BToKEE_l1_charge;
   Float_t         BToKEE_l2_vx;
   Float_t         BToKEE_l2_vy;
   Float_t         BToKEE_l2_vz;
   Float_t         BToKEE_l2_charge;
   Float_t         BToKEE_k_vx;
   Float_t         BToKEE_k_vy;
   Float_t         BToKEE_k_vz;
   Float_t         BToKEE_k_dca_sig;
   Float_t         BToKEE_k_dz;
   Float_t         BToKEE_kMu_matched;
   Float_t         BToKEE_k_charge;
   Float_t         BToKEE_l_xy_sig;
   Float_t         BToKEE_l1l2Dz;
   Float_t         BToKEE_lKDz;
   Float_t         BToKEE_l1l2Dr;
   Float_t         BToKEE_lKDr;
   Float_t         BToKEE_kl_massKPi;
   Float_t         BToKEE_p_assymetry;
   Float_t         Presel_BDT;
   Float_t         BDTSCORE_noah;
   Float_t         BDTSCORE_analysis_llkdr;
   Float_t         BDTSCORE_analysis_llkdr_withlow;
   Float_t         BDTSCORE_analysis_llkdr_withlow_withnomllcutmc;
   Float_t         BToKEE_k_svipd2d;
   Float_t         BToKEE_k_svipd2d_err;
   Float_t         BToKEE_l1_dzTrg;
   Float_t         BToKEE_l2_dzTrg;
   Float_t         BToKEE_l1_dxy;
   Float_t         BToKEE_l2_dxy;
   Float_t         BToKEE_l1_dxy_err;
   Float_t         BToKEE_l2_dxy_err;
   Float_t         trig_wgt;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_nBToKEE;   //!
   TBranch        *b_nElectron;   //!
   TBranch        *b_nProbeTracks;   //!
   TBranch        *b_HLT_DoubleEle10_eta1p22_mMax6;   //!
   TBranch        *b_HLT_DoubleEle9p5_eta1p22_mMax6;   //!
   TBranch        *b_HLT_DoubleEle9_eta1p22_mMax6;   //!
   TBranch        *b_HLT_DoubleEle8p5_eta1p22_mMax6;   //!
   TBranch        *b_HLT_DoubleEle8_eta1p22_mMax6;   //!
   TBranch        *b_HLT_DoubleEle7p5_eta1p22_mMax6;   //!
   TBranch        *b_HLT_DoubleEle7_eta1p22_mMax6;   //!
   TBranch        *b_HLT_DoubleEle6p5_eta1p22_mMax6;   //!
   TBranch        *b_HLT_DoubleEle6_eta1p22_mMax6;   //!
   TBranch        *b_HLT_DoubleEle5p5_eta1p22_mMax6;   //!
   TBranch        *b_HLT_DoubleEle5_eta1p22_mMax6;   //!
   TBranch        *b_HLT_DoubleEle4p5_eta1p22_mMax6;   //!
   TBranch        *b_HLT_DoubleEle4_eta1p22_mMax6;   //!
   TBranch        *b_L1_DoubleEG11_er1p2_dR_Max0p6;   //!
   TBranch        *b_L1_DoubleEG10p5_er1p2_dR_Max0p6;   //!
   TBranch        *b_L1_DoubleEG10_er1p2_dR_Max0p6;   //!
   TBranch        *b_L1_DoubleEG9p5_er1p2_dR_Max0p6;   //!
   TBranch        *b_L1_DoubleEG9_er1p2_dR_Max0p7;   //!
   TBranch        *b_L1_DoubleEG8p5_er1p2_dR_Max0p7;   //!
   TBranch        *b_L1_DoubleEG8_er1p2_dR_Max0p7;   //!
   TBranch        *b_L1_DoubleEG7p5_er1p2_dR_Max0p7;   //!
   TBranch        *b_L1_DoubleEG7_er1p2_dR_Max0p8;   //!
   TBranch        *b_L1_DoubleEG6p5_er1p2_dR_Max0p8;   //!
   TBranch        *b_L1_DoubleEG6_er1p2_dR_Max0p8;   //!
   TBranch        *b_L1_DoubleEG5p5_er1p2_dR_Max0p8;   //!
   TBranch        *b_L1_DoubleEG5_er1p2_dR_Max0p9;   //!
   TBranch        *b_L1_DoubleEG4p5_er1p2_dR_Max0p9;   //!
   TBranch        *b_L1_DoubleEG4_er1p2_dR_Max0p9;   //!
   TBranch        *b_L1_4p5_HLT_4p0;   //!
   TBranch        *b_L1_5p0_HLT_4p0;   //!
   TBranch        *b_L1_5p5_HLT_4p0;   //!
   TBranch        *b_L1_5p5_HLT_6p0;   //!
   TBranch        *b_L1_6p0_HLT_4p0;   //!
   TBranch        *b_L1_6p5_HLT_4p5;   //!
   TBranch        *b_L1_7p0_HLT_5p0;   //!
   TBranch        *b_L1_7p5_HLT_5p0;   //!
   TBranch        *b_L1_8p0_HLT_5p0;   //!
   TBranch        *b_L1_8p5_HLT_5p0;   //!
   TBranch        *b_L1_8p5_HLT_5p5;   //!
   TBranch        *b_L1_9p0_HLT_6p0;   //!
   TBranch        *b_L1_10p5_HLT_5p0;   //!
   TBranch        *b_L1_10p5_HLT_6p5;   //!
   TBranch        *b_L1_11p0_HLT_6p5;   //!
   TBranch        *b_trigger_OR;   //!
   TBranch        *b_nTrigObj;   //!
   TBranch        *b_nOtherPV;   //!
   TBranch        *b_PV_ndof;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_chi2;   //!
   TBranch        *b_PV_score;   //!
   TBranch        *b_PV_npvs;   //!
   TBranch        *b_PV_npvsGood;   //!
   TBranch        *b_nSV;   //!
   TBranch        *b_nSkimBToKEE;   //!
   TBranch        *b_BToKEE_fit_pt;   //!
   TBranch        *b_BToKEE_fit_eta;   //!
   TBranch        *b_BToKEE_fit_phi;   //!
   TBranch        *b_BToKEE_D0_mass_LepToK_KToPi;   //!
   TBranch        *b_BToKEE_D0_mass_LepToPi_KToK;   //!
   TBranch        *b_BToKEE_fit_mass;   //!
   TBranch        *b_BToKEE_l_xy;   //!
   TBranch        *b_BToKEE_l_xy_unc;   //!
   TBranch        *b_BToKEE_fit_cos2D;   //!
   TBranch        *b_BToKEE_svprob;   //!
   TBranch        *b_BToKEE_fit_massErr;   //!
   TBranch        *b_BToKEE_b_iso04;   //!
   TBranch        *b_BToKEE_mll_fullfit;   //!
   TBranch        *b_BToKEE_vtx_x;   //!
   TBranch        *b_BToKEE_vtx_y;   //!
   TBranch        *b_BToKEE_vtx_z;   //!
   TBranch        *b_BToKEE_l1Idx;   //!
   TBranch        *b_BToKEE_l2Idx;   //!
   TBranch        *b_BToKEE_kIdx;   //!
   TBranch        *b_BToKEE_fit_k_pt;   //!
   TBranch        *b_BToKEE_fit_k_eta;   //!
   TBranch        *b_BToKEE_fit_k_phi;   //!
   TBranch        *b_BToKEE_fit_l1_pt;   //!
   TBranch        *b_BToKEE_fit_l1_eta;   //!
   TBranch        *b_BToKEE_fit_l1_phi;   //!
   TBranch        *b_BToKEE_fit_l2_pt;   //!
   TBranch        *b_BToKEE_fit_l2_eta;   //!
   TBranch        *b_BToKEE_fit_l2_phi;   //!
   TBranch        *b_BToKEE_l1_iso04;   //!
   TBranch        *b_BToKEE_l2_iso04;   //!
   TBranch        *b_BToKEE_l1_isPF;   //!
   TBranch        *b_BToKEE_l2_isPF;   //!
   TBranch        *b_BToKEE_k_iso04;   //!
   TBranch        *b_BToKEE_l1_isPFoverlap;   //!
   TBranch        *b_BToKEE_l2_isPFoverlap;   //!
   TBranch        *b_BToKEE_l1_LooseId;   //!
   TBranch        *b_BToKEE_l2_LooseId;   //!
   TBranch        *b_BToKEE_l1_MediumId;   //!
   TBranch        *b_BToKEE_l2_MediumId;   //!
   TBranch        *b_BToKEE_l1_TightId;   //!
   TBranch        *b_BToKEE_l2_TightId;   //!
   TBranch        *b_BToKEE_l1_ConvVeto;   //!
   TBranch        *b_BToKEE_l2_ConvVeto;   //!
   TBranch        *b_BToKEE_l1_PFMvaID_Fall17;   //!
   TBranch        *b_BToKEE_l2_PFMvaID_Fall17;   //!
   TBranch        *b_BToKEE_l1_PFMvaID_retrained;   //!
   TBranch        *b_BToKEE_l2_PFMvaID_retrained;   //!
   TBranch        *b_BToKEE_l1_iso04_dca;   //!
   TBranch        *b_BToKEE_l2_iso04_dca;   //!
   TBranch        *b_BToKEE_b_iso04_dca;   //!
   TBranch        *b_BToKEE_k_iso04_dca;   //!
   TBranch        *b_BToKEE_k_svip3d;   //!
   TBranch        *b_BToKEE_k_svip3d_err;   //!
   TBranch        *b_BToKEE_l1_n_isotrk_dca;   //!
   TBranch        *b_BToKEE_l2_n_isotrk_dca;   //!
   TBranch        *b_BToKEE_k_n_isotrk_dca;   //!
   TBranch        *b_BToKEE_l1_vx;   //!
   TBranch        *b_BToKEE_l1_vy;   //!
   TBranch        *b_BToKEE_l1_vz;   //!
   TBranch        *b_BToKEE_l1_charge;   //!
   TBranch        *b_BToKEE_l2_vx;   //!
   TBranch        *b_BToKEE_l2_vy;   //!
   TBranch        *b_BToKEE_l2_vz;   //!
   TBranch        *b_BToKEE_l2_charge;   //!
   TBranch        *b_BToKEE_k_vx;   //!
   TBranch        *b_BToKEE_k_vy;   //!
   TBranch        *b_BToKEE_k_vz;   //!
   TBranch        *b_BToKEE_k_dca_sig;   //!
   TBranch        *b_BToKEE_k_dz;   //!
   TBranch        *b_BToKEE_kMu_matched;   //!
   TBranch        *b_BToKEE_k_charge;   //!
   TBranch        *b_BToKEE_l_xy_sig;   //!
   TBranch        *b_BToKEE_l1l2Dz;   //!
   TBranch        *b_BToKEE_lKDz;   //!
   TBranch        *b_BToKEE_l1l2Dr;   //!
   TBranch        *b_BToKEE_lKDr;   //!
   TBranch        *b_BToKEE_kl_massKPi;   //!
   TBranch        *b_BToKEE_p_assymetry;   //!
   TBranch        *b_Presel_BDT;   //!
   TBranch        *b_BDTSCORE_noah;   //!
   TBranch        *b_BDTSCORE_analysis_llkdr;   //!
   TBranch        *b_BDTSCORE_analysis_llkdr_withlow;   //!
   TBranch        *b_BDTSCORE_analysis_llkdr_withlow_withnomllcutmc;   //!
   TBranch        *b_BToKEE_k_svipd2d;   //!
   TBranch        *b_BToKEE_k_svipd2d_err;   //!
   TBranch        *b_BToKEE_l1_dzTrg;   //!
   TBranch        *b_BToKEE_l2_dzTrg;   //!
   TBranch        *b_BToKEE_l1_dxy;   //!
   TBranch        *b_BToKEE_l2_dxy;   //!
   TBranch        *b_BToKEE_l1_dxy_err;   //!
   TBranch        *b_BToKEE_l2_dxy_err;   //!
   TBranch        *b_trig_wgt;   //!

   GeneralBase(TTree *tree=0);
   virtual ~GeneralBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef GeneralBase_cxx
GeneralBase::GeneralBase(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("tree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("tree","");
      chain->Add("outputty_mc_MC_.root/tree");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

GeneralBase::~GeneralBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t GeneralBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t GeneralBase::LoadTree(Long64_t entry)
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

void GeneralBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("nBToKEE", &nBToKEE, &b_nBToKEE);
   fChain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
   fChain->SetBranchAddress("nProbeTracks", &nProbeTracks, &b_nProbeTracks);
   fChain->SetBranchAddress("HLT_DoubleEle10_eta1p22_mMax6", &HLT_DoubleEle10_eta1p22_mMax6, &b_HLT_DoubleEle10_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle9p5_eta1p22_mMax6", &HLT_DoubleEle9p5_eta1p22_mMax6, &b_HLT_DoubleEle9p5_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle9_eta1p22_mMax6", &HLT_DoubleEle9_eta1p22_mMax6, &b_HLT_DoubleEle9_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle8p5_eta1p22_mMax6", &HLT_DoubleEle8p5_eta1p22_mMax6, &b_HLT_DoubleEle8p5_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle8_eta1p22_mMax6", &HLT_DoubleEle8_eta1p22_mMax6, &b_HLT_DoubleEle8_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle7p5_eta1p22_mMax6", &HLT_DoubleEle7p5_eta1p22_mMax6, &b_HLT_DoubleEle7p5_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle7_eta1p22_mMax6", &HLT_DoubleEle7_eta1p22_mMax6, &b_HLT_DoubleEle7_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle6p5_eta1p22_mMax6", &HLT_DoubleEle6p5_eta1p22_mMax6, &b_HLT_DoubleEle6p5_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle6_eta1p22_mMax6", &HLT_DoubleEle6_eta1p22_mMax6, &b_HLT_DoubleEle6_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle5p5_eta1p22_mMax6", &HLT_DoubleEle5p5_eta1p22_mMax6, &b_HLT_DoubleEle5p5_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle5_eta1p22_mMax6", &HLT_DoubleEle5_eta1p22_mMax6, &b_HLT_DoubleEle5_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle4p5_eta1p22_mMax6", &HLT_DoubleEle4p5_eta1p22_mMax6, &b_HLT_DoubleEle4p5_eta1p22_mMax6);
   fChain->SetBranchAddress("HLT_DoubleEle4_eta1p22_mMax6", &HLT_DoubleEle4_eta1p22_mMax6, &b_HLT_DoubleEle4_eta1p22_mMax6);
   fChain->SetBranchAddress("L1_DoubleEG11_er1p2_dR_Max0p6", &L1_DoubleEG11_er1p2_dR_Max0p6, &b_L1_DoubleEG11_er1p2_dR_Max0p6);
   fChain->SetBranchAddress("L1_DoubleEG10p5_er1p2_dR_Max0p6", &L1_DoubleEG10p5_er1p2_dR_Max0p6, &b_L1_DoubleEG10p5_er1p2_dR_Max0p6);
   fChain->SetBranchAddress("L1_DoubleEG10_er1p2_dR_Max0p6", &L1_DoubleEG10_er1p2_dR_Max0p6, &b_L1_DoubleEG10_er1p2_dR_Max0p6);
   fChain->SetBranchAddress("L1_DoubleEG9p5_er1p2_dR_Max0p6", &L1_DoubleEG9p5_er1p2_dR_Max0p6, &b_L1_DoubleEG9p5_er1p2_dR_Max0p6);
   fChain->SetBranchAddress("L1_DoubleEG9_er1p2_dR_Max0p7", &L1_DoubleEG9_er1p2_dR_Max0p7, &b_L1_DoubleEG9_er1p2_dR_Max0p7);
   fChain->SetBranchAddress("L1_DoubleEG8p5_er1p2_dR_Max0p7", &L1_DoubleEG8p5_er1p2_dR_Max0p7, &b_L1_DoubleEG8p5_er1p2_dR_Max0p7);
   fChain->SetBranchAddress("L1_DoubleEG8_er1p2_dR_Max0p7", &L1_DoubleEG8_er1p2_dR_Max0p7, &b_L1_DoubleEG8_er1p2_dR_Max0p7);
   fChain->SetBranchAddress("L1_DoubleEG7p5_er1p2_dR_Max0p7", &L1_DoubleEG7p5_er1p2_dR_Max0p7, &b_L1_DoubleEG7p5_er1p2_dR_Max0p7);
   fChain->SetBranchAddress("L1_DoubleEG7_er1p2_dR_Max0p8", &L1_DoubleEG7_er1p2_dR_Max0p8, &b_L1_DoubleEG7_er1p2_dR_Max0p8);
   fChain->SetBranchAddress("L1_DoubleEG6p5_er1p2_dR_Max0p8", &L1_DoubleEG6p5_er1p2_dR_Max0p8, &b_L1_DoubleEG6p5_er1p2_dR_Max0p8);
   fChain->SetBranchAddress("L1_DoubleEG6_er1p2_dR_Max0p8", &L1_DoubleEG6_er1p2_dR_Max0p8, &b_L1_DoubleEG6_er1p2_dR_Max0p8);
   fChain->SetBranchAddress("L1_DoubleEG5p5_er1p2_dR_Max0p8", &L1_DoubleEG5p5_er1p2_dR_Max0p8, &b_L1_DoubleEG5p5_er1p2_dR_Max0p8);
   fChain->SetBranchAddress("L1_DoubleEG5_er1p2_dR_Max0p9", &L1_DoubleEG5_er1p2_dR_Max0p9, &b_L1_DoubleEG5_er1p2_dR_Max0p9);
   fChain->SetBranchAddress("L1_DoubleEG4p5_er1p2_dR_Max0p9", &L1_DoubleEG4p5_er1p2_dR_Max0p9, &b_L1_DoubleEG4p5_er1p2_dR_Max0p9);
   fChain->SetBranchAddress("L1_DoubleEG4_er1p2_dR_Max0p9", &L1_DoubleEG4_er1p2_dR_Max0p9, &b_L1_DoubleEG4_er1p2_dR_Max0p9);
   fChain->SetBranchAddress("L1_4p5_HLT_4p0", &L1_4p5_HLT_4p0, &b_L1_4p5_HLT_4p0);
   fChain->SetBranchAddress("L1_5p0_HLT_4p0", &L1_5p0_HLT_4p0, &b_L1_5p0_HLT_4p0);
   fChain->SetBranchAddress("L1_5p5_HLT_4p0", &L1_5p5_HLT_4p0, &b_L1_5p5_HLT_4p0);
   fChain->SetBranchAddress("L1_5p5_HLT_6p0", &L1_5p5_HLT_6p0, &b_L1_5p5_HLT_6p0);
   fChain->SetBranchAddress("L1_6p0_HLT_4p0", &L1_6p0_HLT_4p0, &b_L1_6p0_HLT_4p0);
   fChain->SetBranchAddress("L1_6p5_HLT_4p5", &L1_6p5_HLT_4p5, &b_L1_6p5_HLT_4p5);
   fChain->SetBranchAddress("L1_7p0_HLT_5p0", &L1_7p0_HLT_5p0, &b_L1_7p0_HLT_5p0);
   fChain->SetBranchAddress("L1_7p5_HLT_5p0", &L1_7p5_HLT_5p0, &b_L1_7p5_HLT_5p0);
   fChain->SetBranchAddress("L1_8p0_HLT_5p0", &L1_8p0_HLT_5p0, &b_L1_8p0_HLT_5p0);
   fChain->SetBranchAddress("L1_8p5_HLT_5p0", &L1_8p5_HLT_5p0, &b_L1_8p5_HLT_5p0);
   fChain->SetBranchAddress("L1_8p5_HLT_5p5", &L1_8p5_HLT_5p5, &b_L1_8p5_HLT_5p5);
   fChain->SetBranchAddress("L1_9p0_HLT_6p0", &L1_9p0_HLT_6p0, &b_L1_9p0_HLT_6p0);
   fChain->SetBranchAddress("L1_10p5_HLT_5p0", &L1_10p5_HLT_5p0, &b_L1_10p5_HLT_5p0);
   fChain->SetBranchAddress("L1_10p5_HLT_6p5", &L1_10p5_HLT_6p5, &b_L1_10p5_HLT_6p5);
   fChain->SetBranchAddress("L1_11p0_HLT_6p5", &L1_11p0_HLT_6p5, &b_L1_11p0_HLT_6p5);
   fChain->SetBranchAddress("trigger_OR", &trigger_OR, &b_trigger_OR);
   fChain->SetBranchAddress("nTrigObj", &nTrigObj, &b_nTrigObj);
   fChain->SetBranchAddress("nOtherPV", &nOtherPV, &b_nOtherPV);
   fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   fChain->SetBranchAddress("nSV", &nSV, &b_nSV);
   fChain->SetBranchAddress("nSkimBToKEE", &nSkimBToKEE, &b_nSkimBToKEE);
   fChain->SetBranchAddress("BToKEE_fit_pt", &BToKEE_fit_pt, &b_BToKEE_fit_pt);
   fChain->SetBranchAddress("BToKEE_fit_eta", &BToKEE_fit_eta, &b_BToKEE_fit_eta);
   fChain->SetBranchAddress("BToKEE_fit_phi", &BToKEE_fit_phi, &b_BToKEE_fit_phi);
   fChain->SetBranchAddress("BToKEE_D0_mass_LepToK_KToPi", &BToKEE_D0_mass_LepToK_KToPi, &b_BToKEE_D0_mass_LepToK_KToPi);
   fChain->SetBranchAddress("BToKEE_D0_mass_LepToPi_KToK", &BToKEE_D0_mass_LepToPi_KToK, &b_BToKEE_D0_mass_LepToPi_KToK);
   fChain->SetBranchAddress("BToKEE_fit_mass", &BToKEE_fit_mass, &b_BToKEE_fit_mass);
   fChain->SetBranchAddress("BToKEE_l_xy", &BToKEE_l_xy, &b_BToKEE_l_xy);
   fChain->SetBranchAddress("BToKEE_l_xy_unc", &BToKEE_l_xy_unc, &b_BToKEE_l_xy_unc);
   fChain->SetBranchAddress("BToKEE_fit_cos2D", &BToKEE_fit_cos2D, &b_BToKEE_fit_cos2D);
   fChain->SetBranchAddress("BToKEE_svprob", &BToKEE_svprob, &b_BToKEE_svprob);
   fChain->SetBranchAddress("BToKEE_fit_massErr", &BToKEE_fit_massErr, &b_BToKEE_fit_massErr);
   fChain->SetBranchAddress("BToKEE_b_iso04", &BToKEE_b_iso04, &b_BToKEE_b_iso04);
   fChain->SetBranchAddress("BToKEE_mll_fullfit", &BToKEE_mll_fullfit, &b_BToKEE_mll_fullfit);
   fChain->SetBranchAddress("BToKEE_vtx_x", &BToKEE_vtx_x, &b_BToKEE_vtx_x);
   fChain->SetBranchAddress("BToKEE_vtx_y", &BToKEE_vtx_y, &b_BToKEE_vtx_y);
   fChain->SetBranchAddress("BToKEE_vtx_z", &BToKEE_vtx_z, &b_BToKEE_vtx_z);
   fChain->SetBranchAddress("BToKEE_l1Idx", &BToKEE_l1Idx, &b_BToKEE_l1Idx);
   fChain->SetBranchAddress("BToKEE_l2Idx", &BToKEE_l2Idx, &b_BToKEE_l2Idx);
   fChain->SetBranchAddress("BToKEE_kIdx", &BToKEE_kIdx, &b_BToKEE_kIdx);
   fChain->SetBranchAddress("BToKEE_fit_k_pt", &BToKEE_fit_k_pt, &b_BToKEE_fit_k_pt);
   fChain->SetBranchAddress("BToKEE_fit_k_eta", &BToKEE_fit_k_eta, &b_BToKEE_fit_k_eta);
   fChain->SetBranchAddress("BToKEE_fit_k_phi", &BToKEE_fit_k_phi, &b_BToKEE_fit_k_phi);
   fChain->SetBranchAddress("BToKEE_fit_l1_pt", &BToKEE_fit_l1_pt, &b_BToKEE_fit_l1_pt);
   fChain->SetBranchAddress("BToKEE_fit_l1_eta", &BToKEE_fit_l1_eta, &b_BToKEE_fit_l1_eta);
   fChain->SetBranchAddress("BToKEE_fit_l1_phi", &BToKEE_fit_l1_phi, &b_BToKEE_fit_l1_phi);
   fChain->SetBranchAddress("BToKEE_fit_l2_pt", &BToKEE_fit_l2_pt, &b_BToKEE_fit_l2_pt);
   fChain->SetBranchAddress("BToKEE_fit_l2_eta", &BToKEE_fit_l2_eta, &b_BToKEE_fit_l2_eta);
   fChain->SetBranchAddress("BToKEE_fit_l2_phi", &BToKEE_fit_l2_phi, &b_BToKEE_fit_l2_phi);
   fChain->SetBranchAddress("BToKEE_l1_iso04", &BToKEE_l1_iso04, &b_BToKEE_l1_iso04);
   fChain->SetBranchAddress("BToKEE_l2_iso04", &BToKEE_l2_iso04, &b_BToKEE_l2_iso04);
   fChain->SetBranchAddress("BToKEE_l1_isPF", &BToKEE_l1_isPF, &b_BToKEE_l1_isPF);
   fChain->SetBranchAddress("BToKEE_l2_isPF", &BToKEE_l2_isPF, &b_BToKEE_l2_isPF);
   fChain->SetBranchAddress("BToKEE_k_iso04", &BToKEE_k_iso04, &b_BToKEE_k_iso04);
   fChain->SetBranchAddress("BToKEE_l1_isPFoverlap", &BToKEE_l1_isPFoverlap, &b_BToKEE_l1_isPFoverlap);
   fChain->SetBranchAddress("BToKEE_l2_isPFoverlap", &BToKEE_l2_isPFoverlap, &b_BToKEE_l2_isPFoverlap);
   fChain->SetBranchAddress("BToKEE_l1_LooseId", &BToKEE_l1_LooseId, &b_BToKEE_l1_LooseId);
   fChain->SetBranchAddress("BToKEE_l2_LooseId", &BToKEE_l2_LooseId, &b_BToKEE_l2_LooseId);
   fChain->SetBranchAddress("BToKEE_l1_MediumId", &BToKEE_l1_MediumId, &b_BToKEE_l1_MediumId);
   fChain->SetBranchAddress("BToKEE_l2_MediumId", &BToKEE_l2_MediumId, &b_BToKEE_l2_MediumId);
   fChain->SetBranchAddress("BToKEE_l1_TightId", &BToKEE_l1_TightId, &b_BToKEE_l1_TightId);
   fChain->SetBranchAddress("BToKEE_l2_TightId", &BToKEE_l2_TightId, &b_BToKEE_l2_TightId);
   fChain->SetBranchAddress("BToKEE_l1_ConvVeto", &BToKEE_l1_ConvVeto, &b_BToKEE_l1_ConvVeto);
   fChain->SetBranchAddress("BToKEE_l2_ConvVeto", &BToKEE_l2_ConvVeto, &b_BToKEE_l2_ConvVeto);
   fChain->SetBranchAddress("BToKEE_l1_PFMvaID_Fall17", &BToKEE_l1_PFMvaID_Fall17, &b_BToKEE_l1_PFMvaID_Fall17);
   fChain->SetBranchAddress("BToKEE_l2_PFMvaID_Fall17", &BToKEE_l2_PFMvaID_Fall17, &b_BToKEE_l2_PFMvaID_Fall17);
   fChain->SetBranchAddress("BToKEE_l1_PFMvaID_retrained", &BToKEE_l1_PFMvaID_retrained, &b_BToKEE_l1_PFMvaID_retrained);
   fChain->SetBranchAddress("BToKEE_l2_PFMvaID_retrained", &BToKEE_l2_PFMvaID_retrained, &b_BToKEE_l2_PFMvaID_retrained);
   fChain->SetBranchAddress("BToKEE_l1_iso04_dca", &BToKEE_l1_iso04_dca, &b_BToKEE_l1_iso04_dca);
   fChain->SetBranchAddress("BToKEE_l2_iso04_dca", &BToKEE_l2_iso04_dca, &b_BToKEE_l2_iso04_dca);
   fChain->SetBranchAddress("BToKEE_b_iso04_dca", &BToKEE_b_iso04_dca, &b_BToKEE_b_iso04_dca);
   fChain->SetBranchAddress("BToKEE_k_iso04_dca", &BToKEE_k_iso04_dca, &b_BToKEE_k_iso04_dca);
   fChain->SetBranchAddress("BToKEE_k_svip3d", &BToKEE_k_svip3d, &b_BToKEE_k_svip3d);
   fChain->SetBranchAddress("BToKEE_k_svip3d_err", &BToKEE_k_svip3d_err, &b_BToKEE_k_svip3d_err);
   fChain->SetBranchAddress("BToKEE_l1_n_isotrk_dca", &BToKEE_l1_n_isotrk_dca, &b_BToKEE_l1_n_isotrk_dca);
   fChain->SetBranchAddress("BToKEE_l2_n_isotrk_dca", &BToKEE_l2_n_isotrk_dca, &b_BToKEE_l2_n_isotrk_dca);
   fChain->SetBranchAddress("BToKEE_k_n_isotrk_dca", &BToKEE_k_n_isotrk_dca, &b_BToKEE_k_n_isotrk_dca);
   fChain->SetBranchAddress("BToKEE_l1_vx", &BToKEE_l1_vx, &b_BToKEE_l1_vx);
   fChain->SetBranchAddress("BToKEE_l1_vy", &BToKEE_l1_vy, &b_BToKEE_l1_vy);
   fChain->SetBranchAddress("BToKEE_l1_vz", &BToKEE_l1_vz, &b_BToKEE_l1_vz);
   fChain->SetBranchAddress("BToKEE_l1_charge", &BToKEE_l1_charge, &b_BToKEE_l1_charge);
   fChain->SetBranchAddress("BToKEE_l2_vx", &BToKEE_l2_vx, &b_BToKEE_l2_vx);
   fChain->SetBranchAddress("BToKEE_l2_vy", &BToKEE_l2_vy, &b_BToKEE_l2_vy);
   fChain->SetBranchAddress("BToKEE_l2_vz", &BToKEE_l2_vz, &b_BToKEE_l2_vz);
   fChain->SetBranchAddress("BToKEE_l2_charge", &BToKEE_l2_charge, &b_BToKEE_l2_charge);
   fChain->SetBranchAddress("BToKEE_k_vx", &BToKEE_k_vx, &b_BToKEE_k_vx);
   fChain->SetBranchAddress("BToKEE_k_vy", &BToKEE_k_vy, &b_BToKEE_k_vy);
   fChain->SetBranchAddress("BToKEE_k_vz", &BToKEE_k_vz, &b_BToKEE_k_vz);
   fChain->SetBranchAddress("BToKEE_k_dca_sig", &BToKEE_k_dca_sig, &b_BToKEE_k_dca_sig);
   fChain->SetBranchAddress("BToKEE_k_dz", &BToKEE_k_dz, &b_BToKEE_k_dz);
   fChain->SetBranchAddress("BToKEE_kMu_matched", &BToKEE_kMu_matched, &b_BToKEE_kMu_matched);
   fChain->SetBranchAddress("BToKEE_k_charge", &BToKEE_k_charge, &b_BToKEE_k_charge);
   fChain->SetBranchAddress("BToKEE_l_xy_sig", &BToKEE_l_xy_sig, &b_BToKEE_l_xy_sig);
   fChain->SetBranchAddress("BToKEE_l1l2Dz", &BToKEE_l1l2Dz, &b_BToKEE_l1l2Dz);
   fChain->SetBranchAddress("BToKEE_lKDz", &BToKEE_lKDz, &b_BToKEE_lKDz);
   fChain->SetBranchAddress("BToKEE_l1l2Dr", &BToKEE_l1l2Dr, &b_BToKEE_l1l2Dr);
   fChain->SetBranchAddress("BToKEE_lKDr", &BToKEE_lKDr, &b_BToKEE_lKDr);
   fChain->SetBranchAddress("BToKEE_kl_massKPi", &BToKEE_kl_massKPi, &b_BToKEE_kl_massKPi);
   fChain->SetBranchAddress("BToKEE_p_assymetry", &BToKEE_p_assymetry, &b_BToKEE_p_assymetry);
   fChain->SetBranchAddress("Presel_BDT", &Presel_BDT, &b_Presel_BDT);
   fChain->SetBranchAddress("BDTSCORE_noah", &BDTSCORE_noah, &b_BDTSCORE_noah);
   fChain->SetBranchAddress("BDTSCORE_analysis_llkdr", &BDTSCORE_analysis_llkdr, &b_BDTSCORE_analysis_llkdr);
   fChain->SetBranchAddress("BDTSCORE_analysis_llkdr_withlow", &BDTSCORE_analysis_llkdr_withlow, &b_BDTSCORE_analysis_llkdr_withlow);
   fChain->SetBranchAddress("BDTSCORE_analysis_llkdr_withlow_withnomllcutmc", &BDTSCORE_analysis_llkdr_withlow_withnomllcutmc, &b_BDTSCORE_analysis_llkdr_withlow_withnomllcutmc);
   fChain->SetBranchAddress("BToKEE_k_svipd2d", &BToKEE_k_svipd2d, &b_BToKEE_k_svipd2d);
   fChain->SetBranchAddress("BToKEE_k_svipd2d_err", &BToKEE_k_svipd2d_err, &b_BToKEE_k_svipd2d_err);
   fChain->SetBranchAddress("BToKEE_l1_dzTrg", &BToKEE_l1_dzTrg, &b_BToKEE_l1_dzTrg);
   fChain->SetBranchAddress("BToKEE_l2_dzTrg", &BToKEE_l2_dzTrg, &b_BToKEE_l2_dzTrg);
   fChain->SetBranchAddress("BToKEE_l1_dxy", &BToKEE_l1_dxy, &b_BToKEE_l1_dxy);
   fChain->SetBranchAddress("BToKEE_l2_dxy", &BToKEE_l2_dxy, &b_BToKEE_l2_dxy);
   fChain->SetBranchAddress("BToKEE_l1_dxy_err", &BToKEE_l1_dxy_err, &b_BToKEE_l1_dxy_err);
   fChain->SetBranchAddress("BToKEE_l2_dxy_err", &BToKEE_l2_dxy_err, &b_BToKEE_l2_dxy_err);   
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("trig_wgt", &trig_wgt, &b_trig_wgt);

   Notify();
}

Bool_t GeneralBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GeneralBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GeneralBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Class_cxx
