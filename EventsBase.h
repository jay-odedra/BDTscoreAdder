//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul  4 15:23:15 2023 by ROOT version 6.26/11
// from TChain Events/
//////////////////////////////////////////////////////////

#ifndef EventsBase_h
#define EventsBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string>

// Header file for the classes stored in the TTree if any.

class EventsBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   UInt_t          bunchCrossing;
   UInt_t          nBToKEE;
   UInt_t          nElectron;
   UInt_t          nProbeTracks;
   UChar_t         HLT_DoubleEle10_eta1p22_mMax6;
   UChar_t         HLT_DoubleEle9p5_eta1p22_mMax6;
   UChar_t         HLT_DoubleEle9_eta1p22_mMax6;
   UChar_t         HLT_DoubleEle8p5_eta1p22_mMax6;
   UChar_t         HLT_DoubleEle8_eta1p22_mMax6;
   UChar_t         HLT_DoubleEle7p5_eta1p22_mMax6;
   UChar_t         HLT_DoubleEle7_eta1p22_mMax6;
   UChar_t         HLT_DoubleEle6p5_eta1p22_mMax6;
   UChar_t         HLT_DoubleEle6_eta1p22_mMax6;
   UChar_t         HLT_DoubleEle5p5_eta1p22_mMax6;
   UChar_t         HLT_DoubleEle5_eta1p22_mMax6;
   UChar_t         HLT_DoubleEle4p5_eta1p22_mMax6;
   UChar_t         HLT_DoubleEle4_eta1p22_mMax6;
   UChar_t         L1_DoubleEG11_er1p2_dR_Max0p6;
   UChar_t         L1_DoubleEG10p5_er1p2_dR_Max0p6;
   UChar_t         L1_DoubleEG10_er1p2_dR_Max0p6;
   UChar_t         L1_DoubleEG9p5_er1p2_dR_Max0p6;
   UChar_t         L1_DoubleEG9_er1p2_dR_Max0p7;
   UChar_t         L1_DoubleEG8p5_er1p2_dR_Max0p7;
   UChar_t         L1_DoubleEG8_er1p2_dR_Max0p7;
   UChar_t         L1_DoubleEG7p5_er1p2_dR_Max0p7;
   UChar_t         L1_DoubleEG7_er1p2_dR_Max0p8;
   UChar_t         L1_DoubleEG6p5_er1p2_dR_Max0p8;
   UChar_t         L1_DoubleEG6_er1p2_dR_Max0p8;
   UChar_t         L1_DoubleEG5p5_er1p2_dR_Max0p8;
   UChar_t         L1_DoubleEG5_er1p2_dR_Max0p9;
   UChar_t         L1_DoubleEG4p5_er1p2_dR_Max0p9;
   UChar_t         L1_DoubleEG4_er1p2_dR_Max0p9;
   UInt_t          nTrigObj;
   UInt_t          nOtherPV;
   Float_t         PV_ndof;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_chi2;
   Float_t         PV_score;
   Int_t           PV_npvs;
   Int_t           PV_npvsGood;
   UInt_t          nSV;
   UInt_t          nSkimBToKEE;
   Double_t        SkimBToKEE_fit_pt[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_fit_eta[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_fit_phi[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_D0_mass_LepToK_KToPi[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_D0_mass_LepToPi_KToK[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_fit_mass[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l_xy[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l_xy_unc[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_fit_cos2D[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_svprob[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_fit_massErr[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_b_iso04[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_mll_fullfit[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_vtx_x[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_vtx_y[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_vtx_z[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l1Idx[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l2Idx[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_kIdx[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_fit_k_pt[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_fit_k_eta[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_fit_k_phi[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_fit_l1_pt[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_fit_l1_eta[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_fit_l1_phi[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_fit_l2_pt[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_fit_l2_eta[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_fit_l2_phi[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l1_iso04[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l2_iso04[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l1_isPF[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l2_isPF[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_k_iso04[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l1_isPFoverlap[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l2_isPFoverlap[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l1_LooseId[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l2_LooseId[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l1_MediumId[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l2_MediumId[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l1_TightId[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l2_TightId[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l1_ConvVeto[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l2_ConvVeto[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l1_PFMvaID_Fall17[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l2_PFMvaID_Fall17[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l1_PFMvaID_retrained[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l2_PFMvaID_retrained[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l1_iso04_dca[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l2_iso04_dca[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_b_iso04_dca[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_k_iso04_dca[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_k_svip3d[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_k_svip3d_err[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l1_n_isotrk_dca[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l2_n_isotrk_dca[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_k_n_isotrk_dca[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l1_vx[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l1_vy[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l1_vz[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l1_charge[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l2_vx[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l2_vy[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l2_vz[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l2_charge[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_k_vx[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_k_vy[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_k_vz[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_k_dca_sig[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_k_dz[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_kMu_matched[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_k_charge[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l_xy_sig[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l1l2Dz[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_lKDz[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_l1l2Dr[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_lKDr[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_kl_massKPi[570];   //[nSkimBToKEE]
   Double_t        SkimBToKEE_p_assymetry[570];   //[nSkimBToKEE]

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
   TBranch        *b_SkimBToKEE_fit_pt;   //!
   TBranch        *b_SkimBToKEE_fit_eta;   //!
   TBranch        *b_SkimBToKEE_fit_phi;   //!
   TBranch        *b_SkimBToKEE_D0_mass_LepToK_KToPi;   //!
   TBranch        *b_SkimBToKEE_D0_mass_LepToPi_KToK;   //!
   TBranch        *b_SkimBToKEE_fit_mass;   //!
   TBranch        *b_SkimBToKEE_l_xy;   //!
   TBranch        *b_SkimBToKEE_l_xy_unc;   //!
   TBranch        *b_SkimBToKEE_fit_cos2D;   //!
   TBranch        *b_SkimBToKEE_svprob;   //!
   TBranch        *b_SkimBToKEE_fit_massErr;   //!
   TBranch        *b_SkimBToKEE_b_iso04;   //!
   TBranch        *b_SkimBToKEE_mll_fullfit;   //!
   TBranch        *b_SkimBToKEE_vtx_x;   //!
   TBranch        *b_SkimBToKEE_vtx_y;   //!
   TBranch        *b_SkimBToKEE_vtx_z;   //!
   TBranch        *b_SkimBToKEE_l1Idx;   //!
   TBranch        *b_SkimBToKEE_l2Idx;   //!
   TBranch        *b_SkimBToKEE_kIdx;   //!
   TBranch        *b_SkimBToKEE_fit_k_pt;   //!
   TBranch        *b_SkimBToKEE_fit_k_eta;   //!
   TBranch        *b_SkimBToKEE_fit_k_phi;   //!
   TBranch        *b_SkimBToKEE_fit_l1_pt;   //!
   TBranch        *b_SkimBToKEE_fit_l1_eta;   //!
   TBranch        *b_SkimBToKEE_fit_l1_phi;   //!
   TBranch        *b_SkimBToKEE_fit_l2_pt;   //!
   TBranch        *b_SkimBToKEE_fit_l2_eta;   //!
   TBranch        *b_SkimBToKEE_fit_l2_phi;   //!
   TBranch        *b_SkimBToKEE_l1_iso04;   //!
   TBranch        *b_SkimBToKEE_l2_iso04;   //!
   TBranch        *b_SkimBToKEE_l1_isPF;   //!
   TBranch        *b_SkimBToKEE_l2_isPF;   //!
   TBranch        *b_SkimBToKEE_k_iso04;   //!
   TBranch        *b_SkimBToKEE_l1_isPFoverlap;   //!
   TBranch        *b_SkimBToKEE_l2_isPFoverlap;   //!
   TBranch        *b_SkimBToKEE_l1_LooseId;   //!
   TBranch        *b_SkimBToKEE_l2_LooseId;   //!
   TBranch        *b_SkimBToKEE_l1_MediumId;   //!
   TBranch        *b_SkimBToKEE_l2_MediumId;   //!
   TBranch        *b_SkimBToKEE_l1_TightId;   //!
   TBranch        *b_SkimBToKEE_l2_TightId;   //!
   TBranch        *b_SkimBToKEE_l1_ConvVeto;   //!
   TBranch        *b_SkimBToKEE_l2_ConvVeto;   //!
   TBranch        *b_SkimBToKEE_l1_PFMvaID_Fall17;   //!
   TBranch        *b_SkimBToKEE_l2_PFMvaID_Fall17;   //!
   TBranch        *b_SkimBToKEE_l1_PFMvaID_retrained;   //!
   TBranch        *b_SkimBToKEE_l2_PFMvaID_retrained;   //!
   TBranch        *b_SkimBToKEE_l1_iso04_dca;   //!
   TBranch        *b_SkimBToKEE_l2_iso04_dca;   //!
   TBranch        *b_SkimBToKEE_b_iso04_dca;   //!
   TBranch        *b_SkimBToKEE_k_iso04_dca;   //!
   TBranch        *b_SkimBToKEE_k_svip3d;   //!
   TBranch        *b_SkimBToKEE_k_svip3d_err;   //!
   TBranch        *b_SkimBToKEE_l1_n_isotrk_dca;   //!
   TBranch        *b_SkimBToKEE_l2_n_isotrk_dca;   //!
   TBranch        *b_SkimBToKEE_k_n_isotrk_dca;   //!
   TBranch        *b_SkimBToKEE_l1_vx;   //!
   TBranch        *b_SkimBToKEE_l1_vy;   //!
   TBranch        *b_SkimBToKEE_l1_vz;   //!
   TBranch        *b_SkimBToKEE_l1_charge;   //!
   TBranch        *b_SkimBToKEE_l2_vx;   //!
   TBranch        *b_SkimBToKEE_l2_vy;   //!
   TBranch        *b_SkimBToKEE_l2_vz;   //!
   TBranch        *b_SkimBToKEE_l2_charge;   //!
   TBranch        *b_SkimBToKEE_k_vx;   //!
   TBranch        *b_SkimBToKEE_k_vy;   //!
   TBranch        *b_SkimBToKEE_k_vz;   //!
   TBranch        *b_SkimBToKEE_k_dca_sig;   //!
   TBranch        *b_SkimBToKEE_k_dz;   //!
   TBranch        *b_SkimBToKEE_kMu_matched;   //!
   TBranch        *b_SkimBToKEE_k_charge;   //!
   TBranch        *b_SkimBToKEE_l_xy_sig;   //!
   TBranch        *b_SkimBToKEE_l1l2Dz;   //!
   TBranch        *b_SkimBToKEE_lKDz;   //!
   TBranch        *b_SkimBToKEE_l1l2Dr;   //!
   TBranch        *b_SkimBToKEE_lKDr;   //!
   TBranch        *b_SkimBToKEE_kl_massKPi;   //!
   TBranch        *b_SkimBToKEE_p_assymetry;   //!

   EventsBase(TTree *tree=0);
   virtual ~EventsBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::string outname, std::string outdir);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Events_cxx
EventsBase::EventsBase(TTree *tree) : fChain(0) 
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
      f->GetObject("Events",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("Events","");
      chain->Add("/eos/user/j/jodedra/fastforesttester/CMSSW_13_1_0/src/Run2022F_part1_0000.root_Run2022F_part1_0000_all.root/Events");
      chain->Add("/eos/user/j/jodedra/fastforesttester/CMSSW_13_1_0/src/smallsample_Chunk0.root/Events");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

EventsBase::~EventsBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EventsBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EventsBase::LoadTree(Long64_t entry)
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

void EventsBase::Init(TTree *tree)
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
   fChain->SetBranchAddress("SkimBToKEE_fit_pt", SkimBToKEE_fit_pt, &b_SkimBToKEE_fit_pt);
   fChain->SetBranchAddress("SkimBToKEE_fit_eta", SkimBToKEE_fit_eta, &b_SkimBToKEE_fit_eta);
   fChain->SetBranchAddress("SkimBToKEE_fit_phi", SkimBToKEE_fit_phi, &b_SkimBToKEE_fit_phi);
   fChain->SetBranchAddress("SkimBToKEE_D0_mass_LepToK_KToPi", SkimBToKEE_D0_mass_LepToK_KToPi, &b_SkimBToKEE_D0_mass_LepToK_KToPi);
   fChain->SetBranchAddress("SkimBToKEE_D0_mass_LepToPi_KToK", SkimBToKEE_D0_mass_LepToPi_KToK, &b_SkimBToKEE_D0_mass_LepToPi_KToK);
   fChain->SetBranchAddress("SkimBToKEE_fit_mass", SkimBToKEE_fit_mass, &b_SkimBToKEE_fit_mass);
   fChain->SetBranchAddress("SkimBToKEE_l_xy", SkimBToKEE_l_xy, &b_SkimBToKEE_l_xy);
   fChain->SetBranchAddress("SkimBToKEE_l_xy_unc", SkimBToKEE_l_xy_unc, &b_SkimBToKEE_l_xy_unc);
   fChain->SetBranchAddress("SkimBToKEE_fit_cos2D", SkimBToKEE_fit_cos2D, &b_SkimBToKEE_fit_cos2D);
   fChain->SetBranchAddress("SkimBToKEE_svprob", SkimBToKEE_svprob, &b_SkimBToKEE_svprob);
   fChain->SetBranchAddress("SkimBToKEE_fit_massErr", SkimBToKEE_fit_massErr, &b_SkimBToKEE_fit_massErr);
   fChain->SetBranchAddress("SkimBToKEE_b_iso04", SkimBToKEE_b_iso04, &b_SkimBToKEE_b_iso04);
   fChain->SetBranchAddress("SkimBToKEE_mll_fullfit", SkimBToKEE_mll_fullfit, &b_SkimBToKEE_mll_fullfit);
   fChain->SetBranchAddress("SkimBToKEE_vtx_x", SkimBToKEE_vtx_x, &b_SkimBToKEE_vtx_x);
   fChain->SetBranchAddress("SkimBToKEE_vtx_y", SkimBToKEE_vtx_y, &b_SkimBToKEE_vtx_y);
   fChain->SetBranchAddress("SkimBToKEE_vtx_z", SkimBToKEE_vtx_z, &b_SkimBToKEE_vtx_z);
   fChain->SetBranchAddress("SkimBToKEE_l1Idx", SkimBToKEE_l1Idx, &b_SkimBToKEE_l1Idx);
   fChain->SetBranchAddress("SkimBToKEE_l2Idx", SkimBToKEE_l2Idx, &b_SkimBToKEE_l2Idx);
   fChain->SetBranchAddress("SkimBToKEE_kIdx", SkimBToKEE_kIdx, &b_SkimBToKEE_kIdx);
   fChain->SetBranchAddress("SkimBToKEE_fit_k_pt", SkimBToKEE_fit_k_pt, &b_SkimBToKEE_fit_k_pt);
   fChain->SetBranchAddress("SkimBToKEE_fit_k_eta", SkimBToKEE_fit_k_eta, &b_SkimBToKEE_fit_k_eta);
   fChain->SetBranchAddress("SkimBToKEE_fit_k_phi", SkimBToKEE_fit_k_phi, &b_SkimBToKEE_fit_k_phi);
   fChain->SetBranchAddress("SkimBToKEE_fit_l1_pt", SkimBToKEE_fit_l1_pt, &b_SkimBToKEE_fit_l1_pt);
   fChain->SetBranchAddress("SkimBToKEE_fit_l1_eta", SkimBToKEE_fit_l1_eta, &b_SkimBToKEE_fit_l1_eta);
   fChain->SetBranchAddress("SkimBToKEE_fit_l1_phi", SkimBToKEE_fit_l1_phi, &b_SkimBToKEE_fit_l1_phi);
   fChain->SetBranchAddress("SkimBToKEE_fit_l2_pt", SkimBToKEE_fit_l2_pt, &b_SkimBToKEE_fit_l2_pt);
   fChain->SetBranchAddress("SkimBToKEE_fit_l2_eta", SkimBToKEE_fit_l2_eta, &b_SkimBToKEE_fit_l2_eta);
   fChain->SetBranchAddress("SkimBToKEE_fit_l2_phi", SkimBToKEE_fit_l2_phi, &b_SkimBToKEE_fit_l2_phi);
   fChain->SetBranchAddress("SkimBToKEE_l1_iso04", SkimBToKEE_l1_iso04, &b_SkimBToKEE_l1_iso04);
   fChain->SetBranchAddress("SkimBToKEE_l2_iso04", SkimBToKEE_l2_iso04, &b_SkimBToKEE_l2_iso04);
   fChain->SetBranchAddress("SkimBToKEE_l1_isPF", SkimBToKEE_l1_isPF, &b_SkimBToKEE_l1_isPF);
   fChain->SetBranchAddress("SkimBToKEE_l2_isPF", SkimBToKEE_l2_isPF, &b_SkimBToKEE_l2_isPF);
   fChain->SetBranchAddress("SkimBToKEE_k_iso04", SkimBToKEE_k_iso04, &b_SkimBToKEE_k_iso04);
   fChain->SetBranchAddress("SkimBToKEE_l1_isPFoverlap", SkimBToKEE_l1_isPFoverlap, &b_SkimBToKEE_l1_isPFoverlap);
   fChain->SetBranchAddress("SkimBToKEE_l2_isPFoverlap", SkimBToKEE_l2_isPFoverlap, &b_SkimBToKEE_l2_isPFoverlap);
   fChain->SetBranchAddress("SkimBToKEE_l1_LooseId", SkimBToKEE_l1_LooseId, &b_SkimBToKEE_l1_LooseId);
   fChain->SetBranchAddress("SkimBToKEE_l2_LooseId", SkimBToKEE_l2_LooseId, &b_SkimBToKEE_l2_LooseId);
   fChain->SetBranchAddress("SkimBToKEE_l1_MediumId", SkimBToKEE_l1_MediumId, &b_SkimBToKEE_l1_MediumId);
   fChain->SetBranchAddress("SkimBToKEE_l2_MediumId", SkimBToKEE_l2_MediumId, &b_SkimBToKEE_l2_MediumId);
   fChain->SetBranchAddress("SkimBToKEE_l1_TightId", SkimBToKEE_l1_TightId, &b_SkimBToKEE_l1_TightId);
   fChain->SetBranchAddress("SkimBToKEE_l2_TightId", SkimBToKEE_l2_TightId, &b_SkimBToKEE_l2_TightId);
   fChain->SetBranchAddress("SkimBToKEE_l1_ConvVeto", SkimBToKEE_l1_ConvVeto, &b_SkimBToKEE_l1_ConvVeto);
   fChain->SetBranchAddress("SkimBToKEE_l2_ConvVeto", SkimBToKEE_l2_ConvVeto, &b_SkimBToKEE_l2_ConvVeto);
   fChain->SetBranchAddress("SkimBToKEE_l1_PFMvaID_Fall17", SkimBToKEE_l1_PFMvaID_Fall17, &b_SkimBToKEE_l1_PFMvaID_Fall17);
   fChain->SetBranchAddress("SkimBToKEE_l2_PFMvaID_Fall17", SkimBToKEE_l2_PFMvaID_Fall17, &b_SkimBToKEE_l2_PFMvaID_Fall17);
   fChain->SetBranchAddress("SkimBToKEE_l1_PFMvaID_retrained", SkimBToKEE_l1_PFMvaID_retrained, &b_SkimBToKEE_l1_PFMvaID_retrained);
   fChain->SetBranchAddress("SkimBToKEE_l2_PFMvaID_retrained", SkimBToKEE_l2_PFMvaID_retrained, &b_SkimBToKEE_l2_PFMvaID_retrained);
   fChain->SetBranchAddress("SkimBToKEE_l1_iso04_dca", SkimBToKEE_l1_iso04_dca, &b_SkimBToKEE_l1_iso04_dca);
   fChain->SetBranchAddress("SkimBToKEE_l2_iso04_dca", SkimBToKEE_l2_iso04_dca, &b_SkimBToKEE_l2_iso04_dca);
   fChain->SetBranchAddress("SkimBToKEE_b_iso04_dca", SkimBToKEE_b_iso04_dca, &b_SkimBToKEE_b_iso04_dca);
   fChain->SetBranchAddress("SkimBToKEE_k_iso04_dca", SkimBToKEE_k_iso04_dca, &b_SkimBToKEE_k_iso04_dca);
   fChain->SetBranchAddress("SkimBToKEE_k_svip3d", SkimBToKEE_k_svip3d, &b_SkimBToKEE_k_svip3d);
   fChain->SetBranchAddress("SkimBToKEE_k_svip3d_err", SkimBToKEE_k_svip3d_err, &b_SkimBToKEE_k_svip3d_err);
   fChain->SetBranchAddress("SkimBToKEE_l1_n_isotrk_dca", SkimBToKEE_l1_n_isotrk_dca, &b_SkimBToKEE_l1_n_isotrk_dca);
   fChain->SetBranchAddress("SkimBToKEE_l2_n_isotrk_dca", SkimBToKEE_l2_n_isotrk_dca, &b_SkimBToKEE_l2_n_isotrk_dca);
   fChain->SetBranchAddress("SkimBToKEE_k_n_isotrk_dca", SkimBToKEE_k_n_isotrk_dca, &b_SkimBToKEE_k_n_isotrk_dca);
   fChain->SetBranchAddress("SkimBToKEE_l1_vx", SkimBToKEE_l1_vx, &b_SkimBToKEE_l1_vx);
   fChain->SetBranchAddress("SkimBToKEE_l1_vy", SkimBToKEE_l1_vy, &b_SkimBToKEE_l1_vy);
   fChain->SetBranchAddress("SkimBToKEE_l1_vz", SkimBToKEE_l1_vz, &b_SkimBToKEE_l1_vz);
   fChain->SetBranchAddress("SkimBToKEE_l1_charge", SkimBToKEE_l1_charge, &b_SkimBToKEE_l1_charge);
   fChain->SetBranchAddress("SkimBToKEE_l2_vx", SkimBToKEE_l2_vx, &b_SkimBToKEE_l2_vx);
   fChain->SetBranchAddress("SkimBToKEE_l2_vy", SkimBToKEE_l2_vy, &b_SkimBToKEE_l2_vy);
   fChain->SetBranchAddress("SkimBToKEE_l2_vz", SkimBToKEE_l2_vz, &b_SkimBToKEE_l2_vz);
   fChain->SetBranchAddress("SkimBToKEE_l2_charge", SkimBToKEE_l2_charge, &b_SkimBToKEE_l2_charge);
   fChain->SetBranchAddress("SkimBToKEE_k_vx", SkimBToKEE_k_vx, &b_SkimBToKEE_k_vx);
   fChain->SetBranchAddress("SkimBToKEE_k_vy", SkimBToKEE_k_vy, &b_SkimBToKEE_k_vy);
   fChain->SetBranchAddress("SkimBToKEE_k_vz", SkimBToKEE_k_vz, &b_SkimBToKEE_k_vz);
   fChain->SetBranchAddress("SkimBToKEE_k_dca_sig", SkimBToKEE_k_dca_sig, &b_SkimBToKEE_k_dca_sig);
   fChain->SetBranchAddress("SkimBToKEE_k_dz", SkimBToKEE_k_dz, &b_SkimBToKEE_k_dz);
   fChain->SetBranchAddress("SkimBToKEE_kMu_matched", SkimBToKEE_kMu_matched, &b_SkimBToKEE_kMu_matched);
   fChain->SetBranchAddress("SkimBToKEE_k_charge", SkimBToKEE_k_charge, &b_SkimBToKEE_k_charge);
   fChain->SetBranchAddress("SkimBToKEE_l_xy_sig", SkimBToKEE_l_xy_sig, &b_SkimBToKEE_l_xy_sig);
   fChain->SetBranchAddress("SkimBToKEE_l1l2Dz", SkimBToKEE_l1l2Dz, &b_SkimBToKEE_l1l2Dz);
   fChain->SetBranchAddress("SkimBToKEE_lKDz", SkimBToKEE_lKDz, &b_SkimBToKEE_lKDz);
   fChain->SetBranchAddress("SkimBToKEE_l1l2Dr", SkimBToKEE_l1l2Dr, &b_SkimBToKEE_l1l2Dr);
   fChain->SetBranchAddress("SkimBToKEE_lKDr", SkimBToKEE_lKDr, &b_SkimBToKEE_lKDr);
   fChain->SetBranchAddress("SkimBToKEE_kl_massKPi", SkimBToKEE_kl_massKPi, &b_SkimBToKEE_kl_massKPi);
   fChain->SetBranchAddress("SkimBToKEE_p_assymetry", SkimBToKEE_p_assymetry, &b_SkimBToKEE_p_assymetry);
   Notify();
}

Bool_t EventsBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EventsBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EventsBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Events_cxx
