//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul  7 12:41:07 2023 by ROOT version 6.26/11
// from TChain Events/
//////////////////////////////////////////////////////////

#ifndef EventsBaseMC_h
#define EventsBaseMC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
class EventsBaseMC {
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
   UInt_t          nGenPart;
   Float_t         genWeight;
   UInt_t          nPSWeight;
   Float_t         PSWeight[1]; 
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
   Float_t         genB_pdgId;
   Float_t         genE1_pdgId;
   Float_t         genE2_pdgId;
   Float_t         genK_pdgId;
   Float_t         genB_pt;
   Float_t         genE1_pt;
   Float_t         genE2_pt;
   Float_t         genK_pt;
   Float_t         genB_eta;
   Float_t         genE1_eta;
   Float_t         genE2_eta;
   Float_t         genK_eta;
   Float_t         genB_phi;
   Float_t         genE1_phi;
   Float_t         genE2_phi;
   Float_t         genK_phi;
   Float_t         genB_mass;
   Float_t         genE1_mass;
   Float_t         genE2_mass;
   Float_t         genK_mass;
   Float_t         genB_charge;
   Float_t         genE1_charge;
   Float_t         genE2_charge;
   Float_t         genK_charge;
   Float_t         recoE1_DR;
   Float_t         recoE1_Idx;
   Float_t         recoE1_pt;
   Float_t         recoE1_eta;
   Float_t         recoE1_phi;
   Float_t         recoE1_vx;
   Float_t         recoE1_vy;
   Float_t         recoE1_vz;
   Float_t         recoE1_isPF;
   Float_t         recoE1_isPFoverlap;
   Float_t         recoE1_charge;
   Float_t         recoE1_PFEleMvaID_Fall17NoIsoV2RawValue;
   Float_t         recoE1_PFEleMvaID_RetrainedRawValue;
   Float_t         recoE1_PFEleMvaID_Fall17NoIsoV2wpLoose;
   Float_t         recoE1_PFEleMvaID_Fall17NoIsoV2wp90;
   Float_t         recoE1_PFEleMvaID_Fall17NoIsoV2wp80;
   Float_t         recoE1_convVeto;
   Float_t         recoE2_DR;
   Float_t         recoE2_Idx;
   Float_t         recoE2_pt;
   Float_t         recoE2_eta;
   Float_t         recoE2_phi;
   Float_t         recoE2_vx;
   Float_t         recoE2_vy;
   Float_t         recoE2_vz;
   Float_t         recoE2_isPF;
   Float_t         recoE2_isPFoverlap;
   Float_t         recoE2_charge;
   Float_t         recoE2_PFEleMvaID_Fall17NoIsoV2RawValue;
   Float_t         recoE2_PFEleMvaID_RetrainedRawValue;
   Float_t         recoE2_PFEleMvaID_Fall17NoIsoV2wpLoose;
   Float_t         recoE2_PFEleMvaID_Fall17NoIsoV2wp90;
   Float_t         recoE2_PFEleMvaID_Fall17NoIsoV2wp80;
   Float_t         recoE2_convVeto;
   Float_t         recoK_DR;
   Float_t         recoK_Idx;
   Float_t         recoK_pt;
   Float_t         recoK_eta;
   Float_t         recoK_phi;
   Float_t         recoK_vx;
   Float_t         recoK_vy;
   Float_t         recoK_vz;
   Float_t         recoK_DCASig;
   Float_t         recoK_dzTrg;
   Float_t         recoK_isMatchedToMuon;
   Float_t         recoK_charge;
   UInt_t          npreSkimBToKEE;
   Double_t        preSkimBToKEE_fit_pt[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_fit_eta[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_fit_phi[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_D0_mass_LepToK_KToPi[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_D0_mass_LepToPi_KToK[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_fit_mass[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l_xy[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l_xy_unc[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_fit_cos2D[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_svprob[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_fit_massErr[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_b_iso04[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_mll_fullfit[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_vtx_x[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_vtx_y[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_vtx_z[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l1Idx[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l2Idx[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_kIdx[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_fit_k_pt[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_fit_k_eta[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_fit_k_phi[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_fit_l1_pt[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_fit_l1_eta[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_fit_l1_phi[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_fit_l2_pt[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_fit_l2_eta[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_fit_l2_phi[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l1_iso04[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l2_iso04[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l1_isPF[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l2_isPF[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_k_iso04[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l1_isPFoverlap[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l2_isPFoverlap[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l1_LooseId[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l2_LooseId[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l1_MediumId[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l2_MediumId[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l1_TightId[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l2_TightId[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l1_ConvVeto[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l2_ConvVeto[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l1_PFMvaID_Fall17[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l2_PFMvaID_Fall17[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l1_PFMvaID_retrained[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l2_PFMvaID_retrained[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l1_iso04_dca[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l2_iso04_dca[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_b_iso04_dca[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_k_iso04_dca[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_k_svip3d[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_k_svip3d_err[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l1_n_isotrk_dca[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_l2_n_isotrk_dca[302];   //[npreSkimBToKEE]
   Double_t        preSkimBToKEE_k_n_isotrk_dca[302];   //[npreSkimBToKEE]
   Float_t         recoB_fit_pt;
   Float_t         recoB_fit_eta;
   Float_t         recoB_fit_phi;
   Float_t         recoB_D0_mass_LepToK_KToPi;
   Float_t         recoB_D0_mass_LepToPi_KToK;
   Float_t         recoB_fit_mass;
   Float_t         recoB_l_xy;
   Float_t         recoB_l_xy_unc;
   Float_t         recoB_fit_cos2D;
   Float_t         recoB_svprob;
   Float_t         recoB_fit_massErr;
   Float_t         recoB_b_iso04;
   Float_t         recoB_mll_fullfit;
   Float_t         recoB_vtx_x;
   Float_t         recoB_vtx_y;
   Float_t         recoB_vtx_z;
   Float_t         recoB_l1Idx;
   Float_t         recoB_l2Idx;
   Float_t         recoB_kIdx;
   Float_t         recoB_fit_k_pt;
   Float_t         recoB_fit_k_eta;
   Float_t         recoB_fit_k_phi;
   Float_t         recoB_fit_l1_pt;
   Float_t         recoB_fit_l1_eta;
   Float_t         recoB_fit_l1_phi;
   Float_t         recoB_fit_l2_pt;
   Float_t         recoB_fit_l2_eta;
   Float_t         recoB_fit_l2_phi;
   Float_t         recoB_l1_iso04;
   Float_t         recoB_l2_iso04;
   Float_t         recoB_l1_isPF;
   Float_t         recoB_l2_isPF;
   Float_t         recoB_k_iso04;
   Float_t         recoB_l1_isPFoverlap;
   Float_t         recoB_l2_isPFoverlap;
   Float_t         recoB_l1_LooseId;
   Float_t         recoB_l2_LooseId;
   Float_t         recoB_l1_MediumId;
   Float_t         recoB_l2_MediumId;
   Float_t         recoB_l1_TightId;
   Float_t         recoB_l2_TightId;
   Float_t         recoB_l1_ConvVeto;
   Float_t         recoB_l2_ConvVeto;
   Float_t         recoB_l1_PFMvaID_Fall17;
   Float_t         recoB_l2_PFMvaID_Fall17;
   Float_t         recoB_l1_PFMvaID_retrained;
   Float_t         recoB_l2_PFMvaID_retrained;
   Float_t         recoB_l1_iso04_dca;
   Float_t         recoB_l2_iso04_dca;
   Float_t         recoB_b_iso04_dca;
   Float_t         recoB_k_iso04_dca;
   Float_t         recoB_k_svip3d;
   Float_t         recoB_k_svip3d_err;
   Float_t         recoB_l1_n_isotrk_dca;
   Float_t         recoB_l2_n_isotrk_dca;
   Float_t         recoB_k_n_isotrk_dca;
   Float_t         recoB_Idx;
   Float_t         recoB_l_xy_sig;
   Float_t         recoB_l1l2Dz;
   Float_t         recoB_lKDz;
   Float_t         recoB_l1l2Dr;
   Float_t         recoB_lKDr;
   Float_t         recoB_k_opp_l_mass;
   Float_t         recoB_p_assymetry;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_luminosityBlock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_nBToKEE;   //!
   TBranch        *b_nElectron;   //!
   TBranch        *b_nGenPart;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_nPSWeight;   //!
   TBranch        *b_PSWeight;   //!
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
   TBranch        *b_genB_pdgId;   //!
   TBranch        *b_genE1_pdgId;   //!
   TBranch        *b_genE2_pdgId;   //!
   TBranch        *b_genK_pdgId;   //!
   TBranch        *b_genB_pt;   //!
   TBranch        *b_genE1_pt;   //!
   TBranch        *b_genE2_pt;   //!
   TBranch        *b_genK_pt;   //!
   TBranch        *b_genB_eta;   //!
   TBranch        *b_genE1_eta;   //!
   TBranch        *b_genE2_eta;   //!
   TBranch        *b_genK_eta;   //!
   TBranch        *b_genB_phi;   //!
   TBranch        *b_genE1_phi;   //!
   TBranch        *b_genE2_phi;   //!
   TBranch        *b_genK_phi;   //!
   TBranch        *b_genB_mass;   //!
   TBranch        *b_genE1_mass;   //!
   TBranch        *b_genE2_mass;   //!
   TBranch        *b_genK_mass;   //!
   TBranch        *b_genB_charge;   //!
   TBranch        *b_genE1_charge;   //!
   TBranch        *b_genE2_charge;   //!
   TBranch        *b_genK_charge;   //!
   TBranch        *b_recoE1_DR;   //!
   TBranch        *b_recoE1_Idx;   //!
   TBranch        *b_recoE1_pt;   //!
   TBranch        *b_recoE1_eta;   //!
   TBranch        *b_recoE1_phi;   //!
   TBranch        *b_recoE1_vx;   //!
   TBranch        *b_recoE1_vy;   //!
   TBranch        *b_recoE1_vz;   //!
   TBranch        *b_recoE1_isPF;   //!
   TBranch        *b_recoE1_isPFoverlap;   //!
   TBranch        *b_recoE1_charge;   //!
   TBranch        *b_recoE1_PFEleMvaID_Fall17NoIsoV2RawValue;   //!
   TBranch        *b_recoE1_PFEleMvaID_RetrainedRawValue;   //!
   TBranch        *b_recoE1_PFEleMvaID_Fall17NoIsoV2wpLoose;   //!
   TBranch        *b_recoE1_PFEleMvaID_Fall17NoIsoV2wp90;   //!
   TBranch        *b_recoE1_PFEleMvaID_Fall17NoIsoV2wp80;   //!
   TBranch        *b_recoE1_convVeto;   //!
   TBranch        *b_recoE2_DR;   //!
   TBranch        *b_recoE2_Idx;   //!
   TBranch        *b_recoE2_pt;   //!
   TBranch        *b_recoE2_eta;   //!
   TBranch        *b_recoE2_phi;   //!
   TBranch        *b_recoE2_vx;   //!
   TBranch        *b_recoE2_vy;   //!
   TBranch        *b_recoE2_vz;   //!
   TBranch        *b_recoE2_isPF;   //!
   TBranch        *b_recoE2_isPFoverlap;   //!
   TBranch        *b_recoE2_charge;   //!
   TBranch        *b_recoE2_PFEleMvaID_Fall17NoIsoV2RawValue;   //!
   TBranch        *b_recoE2_PFEleMvaID_RetrainedRawValue;   //!
   TBranch        *b_recoE2_PFEleMvaID_Fall17NoIsoV2wpLoose;   //!
   TBranch        *b_recoE2_PFEleMvaID_Fall17NoIsoV2wp90;   //!
   TBranch        *b_recoE2_PFEleMvaID_Fall17NoIsoV2wp80;   //!
   TBranch        *b_recoE2_convVeto;   //!
   TBranch        *b_recoK_DR;   //!
   TBranch        *b_recoK_Idx;   //!
   TBranch        *b_recoK_pt;   //!
   TBranch        *b_recoK_eta;   //!
   TBranch        *b_recoK_phi;   //!
   TBranch        *b_recoK_vx;   //!
   TBranch        *b_recoK_vy;   //!
   TBranch        *b_recoK_vz;   //!
   TBranch        *b_recoK_DCASig;   //!
   TBranch        *b_recoK_dzTrg;   //!
   TBranch        *b_recoK_isMatchedToMuon;   //!
   TBranch        *b_recoK_charge;   //!
   TBranch        *b_npreSkimBToKEE;   //!
   TBranch        *b_preSkimBToKEE_fit_pt;   //!
   TBranch        *b_preSkimBToKEE_fit_eta;   //!
   TBranch        *b_preSkimBToKEE_fit_phi;   //!
   TBranch        *b_preSkimBToKEE_D0_mass_LepToK_KToPi;   //!
   TBranch        *b_preSkimBToKEE_D0_mass_LepToPi_KToK;   //!
   TBranch        *b_preSkimBToKEE_fit_mass;   //!
   TBranch        *b_preSkimBToKEE_l_xy;   //!
   TBranch        *b_preSkimBToKEE_l_xy_unc;   //!
   TBranch        *b_preSkimBToKEE_fit_cos2D;   //!
   TBranch        *b_preSkimBToKEE_svprob;   //!
   TBranch        *b_preSkimBToKEE_fit_massErr;   //!
   TBranch        *b_preSkimBToKEE_b_iso04;   //!
   TBranch        *b_preSkimBToKEE_mll_fullfit;   //!
   TBranch        *b_preSkimBToKEE_vtx_x;   //!
   TBranch        *b_preSkimBToKEE_vtx_y;   //!
   TBranch        *b_preSkimBToKEE_vtx_z;   //!
   TBranch        *b_preSkimBToKEE_l1Idx;   //!
   TBranch        *b_preSkimBToKEE_l2Idx;   //!
   TBranch        *b_preSkimBToKEE_kIdx;   //!
   TBranch        *b_preSkimBToKEE_fit_k_pt;   //!
   TBranch        *b_preSkimBToKEE_fit_k_eta;   //!
   TBranch        *b_preSkimBToKEE_fit_k_phi;   //!
   TBranch        *b_preSkimBToKEE_fit_l1_pt;   //!
   TBranch        *b_preSkimBToKEE_fit_l1_eta;   //!
   TBranch        *b_preSkimBToKEE_fit_l1_phi;   //!
   TBranch        *b_preSkimBToKEE_fit_l2_pt;   //!
   TBranch        *b_preSkimBToKEE_fit_l2_eta;   //!
   TBranch        *b_preSkimBToKEE_fit_l2_phi;   //!
   TBranch        *b_preSkimBToKEE_l1_iso04;   //!
   TBranch        *b_preSkimBToKEE_l2_iso04;   //!
   TBranch        *b_preSkimBToKEE_l1_isPF;   //!
   TBranch        *b_preSkimBToKEE_l2_isPF;   //!
   TBranch        *b_preSkimBToKEE_k_iso04;   //!
   TBranch        *b_preSkimBToKEE_l1_isPFoverlap;   //!
   TBranch        *b_preSkimBToKEE_l2_isPFoverlap;   //!
   TBranch        *b_preSkimBToKEE_l1_LooseId;   //!
   TBranch        *b_preSkimBToKEE_l2_LooseId;   //!
   TBranch        *b_preSkimBToKEE_l1_MediumId;   //!
   TBranch        *b_preSkimBToKEE_l2_MediumId;   //!
   TBranch        *b_preSkimBToKEE_l1_TightId;   //!
   TBranch        *b_preSkimBToKEE_l2_TightId;   //!
   TBranch        *b_preSkimBToKEE_l1_ConvVeto;   //!
   TBranch        *b_preSkimBToKEE_l2_ConvVeto;   //!
   TBranch        *b_preSkimBToKEE_l1_PFMvaID_Fall17;   //!
   TBranch        *b_preSkimBToKEE_l2_PFMvaID_Fall17;   //!
   TBranch        *b_preSkimBToKEE_l1_PFMvaID_retrained;   //!
   TBranch        *b_preSkimBToKEE_l2_PFMvaID_retrained;   //!
   TBranch        *b_preSkimBToKEE_l1_iso04_dca;   //!
   TBranch        *b_preSkimBToKEE_l2_iso04_dca;   //!
   TBranch        *b_preSkimBToKEE_b_iso04_dca;   //!
   TBranch        *b_preSkimBToKEE_k_iso04_dca;   //!
   TBranch        *b_preSkimBToKEE_k_svip3d;   //!
   TBranch        *b_preSkimBToKEE_k_svip3d_err;   //!
   TBranch        *b_preSkimBToKEE_l1_n_isotrk_dca;   //!
   TBranch        *b_preSkimBToKEE_l2_n_isotrk_dca;   //!
   TBranch        *b_preSkimBToKEE_k_n_isotrk_dca;   //!
   TBranch        *b_recoB_fit_pt;   //!
   TBranch        *b_recoB_fit_eta;   //!
   TBranch        *b_recoB_fit_phi;   //!
   TBranch        *b_recoB_D0_mass_LepToK_KToPi;   //!
   TBranch        *b_recoB_D0_mass_LepToPi_KToK;   //!
   TBranch        *b_recoB_fit_mass;   //!
   TBranch        *b_recoB_l_xy;   //!
   TBranch        *b_recoB_l_xy_unc;   //!
   TBranch        *b_recoB_fit_cos2D;   //!
   TBranch        *b_recoB_svprob;   //!
   TBranch        *b_recoB_fit_massErr;   //!
   TBranch        *b_recoB_b_iso04;   //!
   TBranch        *b_recoB_mll_fullfit;   //!
   TBranch        *b_recoB_vtx_x;   //!
   TBranch        *b_recoB_vtx_y;   //!
   TBranch        *b_recoB_vtx_z;   //!
   TBranch        *b_recoB_l1Idx;   //!
   TBranch        *b_recoB_l2Idx;   //!
   TBranch        *b_recoB_kIdx;   //!
   TBranch        *b_recoB_fit_k_pt;   //!
   TBranch        *b_recoB_fit_k_eta;   //!
   TBranch        *b_recoB_fit_k_phi;   //!
   TBranch        *b_recoB_fit_l1_pt;   //!
   TBranch        *b_recoB_fit_l1_eta;   //!
   TBranch        *b_recoB_fit_l1_phi;   //!
   TBranch        *b_recoB_fit_l2_pt;   //!
   TBranch        *b_recoB_fit_l2_eta;   //!
   TBranch        *b_recoB_fit_l2_phi;   //!
   TBranch        *b_recoB_l1_iso04;   //!
   TBranch        *b_recoB_l2_iso04;   //!
   TBranch        *b_recoB_l1_isPF;   //!
   TBranch        *b_recoB_l2_isPF;   //!
   TBranch        *b_recoB_k_iso04;   //!
   TBranch        *b_recoB_l1_isPFoverlap;   //!
   TBranch        *b_recoB_l2_isPFoverlap;   //!
   TBranch        *b_recoB_l1_LooseId;   //!
   TBranch        *b_recoB_l2_LooseId;   //!
   TBranch        *b_recoB_l1_MediumId;   //!
   TBranch        *b_recoB_l2_MediumId;   //!
   TBranch        *b_recoB_l1_TightId;   //!
   TBranch        *b_recoB_l2_TightId;   //!
   TBranch        *b_recoB_l1_ConvVeto;   //!
   TBranch        *b_recoB_l2_ConvVeto;   //!
   TBranch        *b_recoB_l1_PFMvaID_Fall17;   //!
   TBranch        *b_recoB_l2_PFMvaID_Fall17;   //!
   TBranch        *b_recoB_l1_PFMvaID_retrained;   //!
   TBranch        *b_recoB_l2_PFMvaID_retrained;   //!
   TBranch        *b_recoB_l1_iso04_dca;   //!
   TBranch        *b_recoB_l2_iso04_dca;   //!
   TBranch        *b_recoB_b_iso04_dca;   //!
   TBranch        *b_recoB_k_iso04_dca;   //!
   TBranch        *b_recoB_k_svip3d;   //!
   TBranch        *b_recoB_k_svip3d_err;   //!
   TBranch        *b_recoB_l1_n_isotrk_dca;   //!
   TBranch        *b_recoB_l2_n_isotrk_dca;   //!
   TBranch        *b_recoB_k_n_isotrk_dca;   //!
   TBranch        *b_recoB_Idx;   //!
   TBranch        *b_recoB_l_xy_sig;   //!
   TBranch        *b_recoB_l1l2Dz;   //!
   TBranch        *b_recoB_lKDz;   //!
   TBranch        *b_recoB_l1l2Dr;   //!
   TBranch        *b_recoB_lKDr;   //!
   TBranch        *b_recoB_k_opp_l_mass;   //!
   TBranch        *b_recoB_p_assymetry;   //!

   EventsBaseMC(TTree *tree=0);
   virtual ~EventsBaseMC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::string outname, std::string outdir);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef EventsBaseMC_cxx
EventsBaseMC::EventsBaseMC(TTree *tree) : fChain(0) 
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
      chain->Add("/eos/user/j/jodedra/fastforesttestergithub/CMSSW_13_1_0/src/BDTscoreAdder/BuToKJPsi_postEE.root_BuToKJPsi_postEE_all.root/Events");
      chain->Add("/eos/user/j/jodedra/fastforesttestergithub/CMSSW_13_1_0/src/BDTscoreAdder/BuToKJPsi_preEE.root_BuToKJPsi_preEE_all.root/Events");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

EventsBaseMC::~EventsBaseMC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EventsBaseMC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EventsBaseMC::LoadTree(Long64_t entry)
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

void EventsBaseMC::Init(TTree *tree)
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
   fChain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("nPSWeight", &nPSWeight, &b_nPSWeight);
   fChain->SetBranchAddress("PSWeight", PSWeight, &b_PSWeight);
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
   fChain->SetBranchAddress("genB_pdgId", &genB_pdgId, &b_genB_pdgId);
   fChain->SetBranchAddress("genE1_pdgId", &genE1_pdgId, &b_genE1_pdgId);
   fChain->SetBranchAddress("genE2_pdgId", &genE2_pdgId, &b_genE2_pdgId);
   fChain->SetBranchAddress("genK_pdgId", &genK_pdgId, &b_genK_pdgId);
   fChain->SetBranchAddress("genB_pt", &genB_pt, &b_genB_pt);
   fChain->SetBranchAddress("genE1_pt", &genE1_pt, &b_genE1_pt);
   fChain->SetBranchAddress("genE2_pt", &genE2_pt, &b_genE2_pt);
   fChain->SetBranchAddress("genK_pt", &genK_pt, &b_genK_pt);
   fChain->SetBranchAddress("genB_eta", &genB_eta, &b_genB_eta);
   fChain->SetBranchAddress("genE1_eta", &genE1_eta, &b_genE1_eta);
   fChain->SetBranchAddress("genE2_eta", &genE2_eta, &b_genE2_eta);
   fChain->SetBranchAddress("genK_eta", &genK_eta, &b_genK_eta);
   fChain->SetBranchAddress("genB_phi", &genB_phi, &b_genB_phi);
   fChain->SetBranchAddress("genE1_phi", &genE1_phi, &b_genE1_phi);
   fChain->SetBranchAddress("genE2_phi", &genE2_phi, &b_genE2_phi);
   fChain->SetBranchAddress("genK_phi", &genK_phi, &b_genK_phi);
   fChain->SetBranchAddress("genB_mass", &genB_mass, &b_genB_mass);
   fChain->SetBranchAddress("genE1_mass", &genE1_mass, &b_genE1_mass);
   fChain->SetBranchAddress("genE2_mass", &genE2_mass, &b_genE2_mass);
   fChain->SetBranchAddress("genK_mass", &genK_mass, &b_genK_mass);
   fChain->SetBranchAddress("genB_charge", &genB_charge, &b_genB_charge);
   fChain->SetBranchAddress("genE1_charge", &genE1_charge, &b_genE1_charge);
   fChain->SetBranchAddress("genE2_charge", &genE2_charge, &b_genE2_charge);
   fChain->SetBranchAddress("genK_charge", &genK_charge, &b_genK_charge);
   fChain->SetBranchAddress("recoE1_DR", &recoE1_DR, &b_recoE1_DR);
   fChain->SetBranchAddress("recoE1_Idx", &recoE1_Idx, &b_recoE1_Idx);
   fChain->SetBranchAddress("recoE1_pt", &recoE1_pt, &b_recoE1_pt);
   fChain->SetBranchAddress("recoE1_eta", &recoE1_eta, &b_recoE1_eta);
   fChain->SetBranchAddress("recoE1_phi", &recoE1_phi, &b_recoE1_phi);
   fChain->SetBranchAddress("recoE1_vx", &recoE1_vx, &b_recoE1_vx);
   fChain->SetBranchAddress("recoE1_vy", &recoE1_vy, &b_recoE1_vy);
   fChain->SetBranchAddress("recoE1_vz", &recoE1_vz, &b_recoE1_vz);
   fChain->SetBranchAddress("recoE1_isPF", &recoE1_isPF, &b_recoE1_isPF);
   fChain->SetBranchAddress("recoE1_isPFoverlap", &recoE1_isPFoverlap, &b_recoE1_isPFoverlap);
   fChain->SetBranchAddress("recoE1_charge", &recoE1_charge, &b_recoE1_charge);
   fChain->SetBranchAddress("recoE1_PFEleMvaID_Fall17NoIsoV2RawValue", &recoE1_PFEleMvaID_Fall17NoIsoV2RawValue, &b_recoE1_PFEleMvaID_Fall17NoIsoV2RawValue);
   fChain->SetBranchAddress("recoE1_PFEleMvaID_RetrainedRawValue", &recoE1_PFEleMvaID_RetrainedRawValue, &b_recoE1_PFEleMvaID_RetrainedRawValue);
   fChain->SetBranchAddress("recoE1_PFEleMvaID_Fall17NoIsoV2wpLoose", &recoE1_PFEleMvaID_Fall17NoIsoV2wpLoose, &b_recoE1_PFEleMvaID_Fall17NoIsoV2wpLoose);
   fChain->SetBranchAddress("recoE1_PFEleMvaID_Fall17NoIsoV2wp90", &recoE1_PFEleMvaID_Fall17NoIsoV2wp90, &b_recoE1_PFEleMvaID_Fall17NoIsoV2wp90);
   fChain->SetBranchAddress("recoE1_PFEleMvaID_Fall17NoIsoV2wp80", &recoE1_PFEleMvaID_Fall17NoIsoV2wp80, &b_recoE1_PFEleMvaID_Fall17NoIsoV2wp80);
   fChain->SetBranchAddress("recoE1_convVeto", &recoE1_convVeto, &b_recoE1_convVeto);
   fChain->SetBranchAddress("recoE2_DR", &recoE2_DR, &b_recoE2_DR);
   fChain->SetBranchAddress("recoE2_Idx", &recoE2_Idx, &b_recoE2_Idx);
   fChain->SetBranchAddress("recoE2_pt", &recoE2_pt, &b_recoE2_pt);
   fChain->SetBranchAddress("recoE2_eta", &recoE2_eta, &b_recoE2_eta);
   fChain->SetBranchAddress("recoE2_phi", &recoE2_phi, &b_recoE2_phi);
   fChain->SetBranchAddress("recoE2_vx", &recoE2_vx, &b_recoE2_vx);
   fChain->SetBranchAddress("recoE2_vy", &recoE2_vy, &b_recoE2_vy);
   fChain->SetBranchAddress("recoE2_vz", &recoE2_vz, &b_recoE2_vz);
   fChain->SetBranchAddress("recoE2_isPF", &recoE2_isPF, &b_recoE2_isPF);
   fChain->SetBranchAddress("recoE2_isPFoverlap", &recoE2_isPFoverlap, &b_recoE2_isPFoverlap);
   fChain->SetBranchAddress("recoE2_charge", &recoE2_charge, &b_recoE2_charge);
   fChain->SetBranchAddress("recoE2_PFEleMvaID_Fall17NoIsoV2RawValue", &recoE2_PFEleMvaID_Fall17NoIsoV2RawValue, &b_recoE2_PFEleMvaID_Fall17NoIsoV2RawValue);
   fChain->SetBranchAddress("recoE2_PFEleMvaID_RetrainedRawValue", &recoE2_PFEleMvaID_RetrainedRawValue, &b_recoE2_PFEleMvaID_RetrainedRawValue);
   fChain->SetBranchAddress("recoE2_PFEleMvaID_Fall17NoIsoV2wpLoose", &recoE2_PFEleMvaID_Fall17NoIsoV2wpLoose, &b_recoE2_PFEleMvaID_Fall17NoIsoV2wpLoose);
   fChain->SetBranchAddress("recoE2_PFEleMvaID_Fall17NoIsoV2wp90", &recoE2_PFEleMvaID_Fall17NoIsoV2wp90, &b_recoE2_PFEleMvaID_Fall17NoIsoV2wp90);
   fChain->SetBranchAddress("recoE2_PFEleMvaID_Fall17NoIsoV2wp80", &recoE2_PFEleMvaID_Fall17NoIsoV2wp80, &b_recoE2_PFEleMvaID_Fall17NoIsoV2wp80);
   fChain->SetBranchAddress("recoE2_convVeto", &recoE2_convVeto, &b_recoE2_convVeto);
   fChain->SetBranchAddress("recoK_DR", &recoK_DR, &b_recoK_DR);
   fChain->SetBranchAddress("recoK_Idx", &recoK_Idx, &b_recoK_Idx);
   fChain->SetBranchAddress("recoK_pt", &recoK_pt, &b_recoK_pt);
   fChain->SetBranchAddress("recoK_eta", &recoK_eta, &b_recoK_eta);
   fChain->SetBranchAddress("recoK_phi", &recoK_phi, &b_recoK_phi);
   fChain->SetBranchAddress("recoK_vx", &recoK_vx, &b_recoK_vx);
   fChain->SetBranchAddress("recoK_vy", &recoK_vy, &b_recoK_vy);
   fChain->SetBranchAddress("recoK_vz", &recoK_vz, &b_recoK_vz);
   fChain->SetBranchAddress("recoK_DCASig", &recoK_DCASig, &b_recoK_DCASig);
   fChain->SetBranchAddress("recoK_dzTrg", &recoK_dzTrg, &b_recoK_dzTrg);
   fChain->SetBranchAddress("recoK_isMatchedToMuon", &recoK_isMatchedToMuon, &b_recoK_isMatchedToMuon);
   fChain->SetBranchAddress("recoK_charge", &recoK_charge, &b_recoK_charge);
   fChain->SetBranchAddress("npreSkimBToKEE", &npreSkimBToKEE, &b_npreSkimBToKEE);
   fChain->SetBranchAddress("preSkimBToKEE_fit_pt", preSkimBToKEE_fit_pt, &b_preSkimBToKEE_fit_pt);
   fChain->SetBranchAddress("preSkimBToKEE_fit_eta", preSkimBToKEE_fit_eta, &b_preSkimBToKEE_fit_eta);
   fChain->SetBranchAddress("preSkimBToKEE_fit_phi", preSkimBToKEE_fit_phi, &b_preSkimBToKEE_fit_phi);
   fChain->SetBranchAddress("preSkimBToKEE_D0_mass_LepToK_KToPi", preSkimBToKEE_D0_mass_LepToK_KToPi, &b_preSkimBToKEE_D0_mass_LepToK_KToPi);
   fChain->SetBranchAddress("preSkimBToKEE_D0_mass_LepToPi_KToK", preSkimBToKEE_D0_mass_LepToPi_KToK, &b_preSkimBToKEE_D0_mass_LepToPi_KToK);
   fChain->SetBranchAddress("preSkimBToKEE_fit_mass", preSkimBToKEE_fit_mass, &b_preSkimBToKEE_fit_mass);
   fChain->SetBranchAddress("preSkimBToKEE_l_xy", preSkimBToKEE_l_xy, &b_preSkimBToKEE_l_xy);
   fChain->SetBranchAddress("preSkimBToKEE_l_xy_unc", preSkimBToKEE_l_xy_unc, &b_preSkimBToKEE_l_xy_unc);
   fChain->SetBranchAddress("preSkimBToKEE_fit_cos2D", preSkimBToKEE_fit_cos2D, &b_preSkimBToKEE_fit_cos2D);
   fChain->SetBranchAddress("preSkimBToKEE_svprob", preSkimBToKEE_svprob, &b_preSkimBToKEE_svprob);
   fChain->SetBranchAddress("preSkimBToKEE_fit_massErr", preSkimBToKEE_fit_massErr, &b_preSkimBToKEE_fit_massErr);
   fChain->SetBranchAddress("preSkimBToKEE_b_iso04", preSkimBToKEE_b_iso04, &b_preSkimBToKEE_b_iso04);
   fChain->SetBranchAddress("preSkimBToKEE_mll_fullfit", preSkimBToKEE_mll_fullfit, &b_preSkimBToKEE_mll_fullfit);
   fChain->SetBranchAddress("preSkimBToKEE_vtx_x", preSkimBToKEE_vtx_x, &b_preSkimBToKEE_vtx_x);
   fChain->SetBranchAddress("preSkimBToKEE_vtx_y", preSkimBToKEE_vtx_y, &b_preSkimBToKEE_vtx_y);
   fChain->SetBranchAddress("preSkimBToKEE_vtx_z", preSkimBToKEE_vtx_z, &b_preSkimBToKEE_vtx_z);
   fChain->SetBranchAddress("preSkimBToKEE_l1Idx", preSkimBToKEE_l1Idx, &b_preSkimBToKEE_l1Idx);
   fChain->SetBranchAddress("preSkimBToKEE_l2Idx", preSkimBToKEE_l2Idx, &b_preSkimBToKEE_l2Idx);
   fChain->SetBranchAddress("preSkimBToKEE_kIdx", preSkimBToKEE_kIdx, &b_preSkimBToKEE_kIdx);
   fChain->SetBranchAddress("preSkimBToKEE_fit_k_pt", preSkimBToKEE_fit_k_pt, &b_preSkimBToKEE_fit_k_pt);
   fChain->SetBranchAddress("preSkimBToKEE_fit_k_eta", preSkimBToKEE_fit_k_eta, &b_preSkimBToKEE_fit_k_eta);
   fChain->SetBranchAddress("preSkimBToKEE_fit_k_phi", preSkimBToKEE_fit_k_phi, &b_preSkimBToKEE_fit_k_phi);
   fChain->SetBranchAddress("preSkimBToKEE_fit_l1_pt", preSkimBToKEE_fit_l1_pt, &b_preSkimBToKEE_fit_l1_pt);
   fChain->SetBranchAddress("preSkimBToKEE_fit_l1_eta", preSkimBToKEE_fit_l1_eta, &b_preSkimBToKEE_fit_l1_eta);
   fChain->SetBranchAddress("preSkimBToKEE_fit_l1_phi", preSkimBToKEE_fit_l1_phi, &b_preSkimBToKEE_fit_l1_phi);
   fChain->SetBranchAddress("preSkimBToKEE_fit_l2_pt", preSkimBToKEE_fit_l2_pt, &b_preSkimBToKEE_fit_l2_pt);
   fChain->SetBranchAddress("preSkimBToKEE_fit_l2_eta", preSkimBToKEE_fit_l2_eta, &b_preSkimBToKEE_fit_l2_eta);
   fChain->SetBranchAddress("preSkimBToKEE_fit_l2_phi", preSkimBToKEE_fit_l2_phi, &b_preSkimBToKEE_fit_l2_phi);
   fChain->SetBranchAddress("preSkimBToKEE_l1_iso04", preSkimBToKEE_l1_iso04, &b_preSkimBToKEE_l1_iso04);
   fChain->SetBranchAddress("preSkimBToKEE_l2_iso04", preSkimBToKEE_l2_iso04, &b_preSkimBToKEE_l2_iso04);
   fChain->SetBranchAddress("preSkimBToKEE_l1_isPF", preSkimBToKEE_l1_isPF, &b_preSkimBToKEE_l1_isPF);
   fChain->SetBranchAddress("preSkimBToKEE_l2_isPF", preSkimBToKEE_l2_isPF, &b_preSkimBToKEE_l2_isPF);
   fChain->SetBranchAddress("preSkimBToKEE_k_iso04", preSkimBToKEE_k_iso04, &b_preSkimBToKEE_k_iso04);
   fChain->SetBranchAddress("preSkimBToKEE_l1_isPFoverlap", preSkimBToKEE_l1_isPFoverlap, &b_preSkimBToKEE_l1_isPFoverlap);
   fChain->SetBranchAddress("preSkimBToKEE_l2_isPFoverlap", preSkimBToKEE_l2_isPFoverlap, &b_preSkimBToKEE_l2_isPFoverlap);
   fChain->SetBranchAddress("preSkimBToKEE_l1_LooseId", preSkimBToKEE_l1_LooseId, &b_preSkimBToKEE_l1_LooseId);
   fChain->SetBranchAddress("preSkimBToKEE_l2_LooseId", preSkimBToKEE_l2_LooseId, &b_preSkimBToKEE_l2_LooseId);
   fChain->SetBranchAddress("preSkimBToKEE_l1_MediumId", preSkimBToKEE_l1_MediumId, &b_preSkimBToKEE_l1_MediumId);
   fChain->SetBranchAddress("preSkimBToKEE_l2_MediumId", preSkimBToKEE_l2_MediumId, &b_preSkimBToKEE_l2_MediumId);
   fChain->SetBranchAddress("preSkimBToKEE_l1_TightId", preSkimBToKEE_l1_TightId, &b_preSkimBToKEE_l1_TightId);
   fChain->SetBranchAddress("preSkimBToKEE_l2_TightId", preSkimBToKEE_l2_TightId, &b_preSkimBToKEE_l2_TightId);
   fChain->SetBranchAddress("preSkimBToKEE_l1_ConvVeto", preSkimBToKEE_l1_ConvVeto, &b_preSkimBToKEE_l1_ConvVeto);
   fChain->SetBranchAddress("preSkimBToKEE_l2_ConvVeto", preSkimBToKEE_l2_ConvVeto, &b_preSkimBToKEE_l2_ConvVeto);
   fChain->SetBranchAddress("preSkimBToKEE_l1_PFMvaID_Fall17", preSkimBToKEE_l1_PFMvaID_Fall17, &b_preSkimBToKEE_l1_PFMvaID_Fall17);
   fChain->SetBranchAddress("preSkimBToKEE_l2_PFMvaID_Fall17", preSkimBToKEE_l2_PFMvaID_Fall17, &b_preSkimBToKEE_l2_PFMvaID_Fall17);
   fChain->SetBranchAddress("preSkimBToKEE_l1_PFMvaID_retrained", preSkimBToKEE_l1_PFMvaID_retrained, &b_preSkimBToKEE_l1_PFMvaID_retrained);
   fChain->SetBranchAddress("preSkimBToKEE_l2_PFMvaID_retrained", preSkimBToKEE_l2_PFMvaID_retrained, &b_preSkimBToKEE_l2_PFMvaID_retrained);
   fChain->SetBranchAddress("preSkimBToKEE_l1_iso04_dca", preSkimBToKEE_l1_iso04_dca, &b_preSkimBToKEE_l1_iso04_dca);
   fChain->SetBranchAddress("preSkimBToKEE_l2_iso04_dca", preSkimBToKEE_l2_iso04_dca, &b_preSkimBToKEE_l2_iso04_dca);
   fChain->SetBranchAddress("preSkimBToKEE_b_iso04_dca", preSkimBToKEE_b_iso04_dca, &b_preSkimBToKEE_b_iso04_dca);
   fChain->SetBranchAddress("preSkimBToKEE_k_iso04_dca", preSkimBToKEE_k_iso04_dca, &b_preSkimBToKEE_k_iso04_dca);
   fChain->SetBranchAddress("preSkimBToKEE_k_svip3d", preSkimBToKEE_k_svip3d, &b_preSkimBToKEE_k_svip3d);
   fChain->SetBranchAddress("preSkimBToKEE_k_svip3d_err", preSkimBToKEE_k_svip3d_err, &b_preSkimBToKEE_k_svip3d_err);
   fChain->SetBranchAddress("preSkimBToKEE_l1_n_isotrk_dca", preSkimBToKEE_l1_n_isotrk_dca, &b_preSkimBToKEE_l1_n_isotrk_dca);
   fChain->SetBranchAddress("preSkimBToKEE_l2_n_isotrk_dca", preSkimBToKEE_l2_n_isotrk_dca, &b_preSkimBToKEE_l2_n_isotrk_dca);
   fChain->SetBranchAddress("preSkimBToKEE_k_n_isotrk_dca", preSkimBToKEE_k_n_isotrk_dca, &b_preSkimBToKEE_k_n_isotrk_dca);
   fChain->SetBranchAddress("recoB_fit_pt", &recoB_fit_pt, &b_recoB_fit_pt);
   fChain->SetBranchAddress("recoB_fit_eta", &recoB_fit_eta, &b_recoB_fit_eta);
   fChain->SetBranchAddress("recoB_fit_phi", &recoB_fit_phi, &b_recoB_fit_phi);
   fChain->SetBranchAddress("recoB_D0_mass_LepToK_KToPi", &recoB_D0_mass_LepToK_KToPi, &b_recoB_D0_mass_LepToK_KToPi);
   fChain->SetBranchAddress("recoB_D0_mass_LepToPi_KToK", &recoB_D0_mass_LepToPi_KToK, &b_recoB_D0_mass_LepToPi_KToK);
   fChain->SetBranchAddress("recoB_fit_mass", &recoB_fit_mass, &b_recoB_fit_mass);
   fChain->SetBranchAddress("recoB_l_xy", &recoB_l_xy, &b_recoB_l_xy);
   fChain->SetBranchAddress("recoB_l_xy_unc", &recoB_l_xy_unc, &b_recoB_l_xy_unc);
   fChain->SetBranchAddress("recoB_fit_cos2D", &recoB_fit_cos2D, &b_recoB_fit_cos2D);
   fChain->SetBranchAddress("recoB_svprob", &recoB_svprob, &b_recoB_svprob);
   fChain->SetBranchAddress("recoB_fit_massErr", &recoB_fit_massErr, &b_recoB_fit_massErr);
   fChain->SetBranchAddress("recoB_b_iso04", &recoB_b_iso04, &b_recoB_b_iso04);
   fChain->SetBranchAddress("recoB_mll_fullfit", &recoB_mll_fullfit, &b_recoB_mll_fullfit);
   fChain->SetBranchAddress("recoB_vtx_x", &recoB_vtx_x, &b_recoB_vtx_x);
   fChain->SetBranchAddress("recoB_vtx_y", &recoB_vtx_y, &b_recoB_vtx_y);
   fChain->SetBranchAddress("recoB_vtx_z", &recoB_vtx_z, &b_recoB_vtx_z);
   fChain->SetBranchAddress("recoB_l1Idx", &recoB_l1Idx, &b_recoB_l1Idx);
   fChain->SetBranchAddress("recoB_l2Idx", &recoB_l2Idx, &b_recoB_l2Idx);
   fChain->SetBranchAddress("recoB_kIdx", &recoB_kIdx, &b_recoB_kIdx);
   fChain->SetBranchAddress("recoB_fit_k_pt", &recoB_fit_k_pt, &b_recoB_fit_k_pt);
   fChain->SetBranchAddress("recoB_fit_k_eta", &recoB_fit_k_eta, &b_recoB_fit_k_eta);
   fChain->SetBranchAddress("recoB_fit_k_phi", &recoB_fit_k_phi, &b_recoB_fit_k_phi);
   fChain->SetBranchAddress("recoB_fit_l1_pt", &recoB_fit_l1_pt, &b_recoB_fit_l1_pt);
   fChain->SetBranchAddress("recoB_fit_l1_eta", &recoB_fit_l1_eta, &b_recoB_fit_l1_eta);
   fChain->SetBranchAddress("recoB_fit_l1_phi", &recoB_fit_l1_phi, &b_recoB_fit_l1_phi);
   fChain->SetBranchAddress("recoB_fit_l2_pt", &recoB_fit_l2_pt, &b_recoB_fit_l2_pt);
   fChain->SetBranchAddress("recoB_fit_l2_eta", &recoB_fit_l2_eta, &b_recoB_fit_l2_eta);
   fChain->SetBranchAddress("recoB_fit_l2_phi", &recoB_fit_l2_phi, &b_recoB_fit_l2_phi);
   fChain->SetBranchAddress("recoB_l1_iso04", &recoB_l1_iso04, &b_recoB_l1_iso04);
   fChain->SetBranchAddress("recoB_l2_iso04", &recoB_l2_iso04, &b_recoB_l2_iso04);
   fChain->SetBranchAddress("recoB_l1_isPF", &recoB_l1_isPF, &b_recoB_l1_isPF);
   fChain->SetBranchAddress("recoB_l2_isPF", &recoB_l2_isPF, &b_recoB_l2_isPF);
   fChain->SetBranchAddress("recoB_k_iso04", &recoB_k_iso04, &b_recoB_k_iso04);
   fChain->SetBranchAddress("recoB_l1_isPFoverlap", &recoB_l1_isPFoverlap, &b_recoB_l1_isPFoverlap);
   fChain->SetBranchAddress("recoB_l2_isPFoverlap", &recoB_l2_isPFoverlap, &b_recoB_l2_isPFoverlap);
   fChain->SetBranchAddress("recoB_l1_LooseId", &recoB_l1_LooseId, &b_recoB_l1_LooseId);
   fChain->SetBranchAddress("recoB_l2_LooseId", &recoB_l2_LooseId, &b_recoB_l2_LooseId);
   fChain->SetBranchAddress("recoB_l1_MediumId", &recoB_l1_MediumId, &b_recoB_l1_MediumId);
   fChain->SetBranchAddress("recoB_l2_MediumId", &recoB_l2_MediumId, &b_recoB_l2_MediumId);
   fChain->SetBranchAddress("recoB_l1_TightId", &recoB_l1_TightId, &b_recoB_l1_TightId);
   fChain->SetBranchAddress("recoB_l2_TightId", &recoB_l2_TightId, &b_recoB_l2_TightId);
   fChain->SetBranchAddress("recoB_l1_ConvVeto", &recoB_l1_ConvVeto, &b_recoB_l1_ConvVeto);
   fChain->SetBranchAddress("recoB_l2_ConvVeto", &recoB_l2_ConvVeto, &b_recoB_l2_ConvVeto);
   fChain->SetBranchAddress("recoB_l1_PFMvaID_Fall17", &recoB_l1_PFMvaID_Fall17, &b_recoB_l1_PFMvaID_Fall17);
   fChain->SetBranchAddress("recoB_l2_PFMvaID_Fall17", &recoB_l2_PFMvaID_Fall17, &b_recoB_l2_PFMvaID_Fall17);
   fChain->SetBranchAddress("recoB_l1_PFMvaID_retrained", &recoB_l1_PFMvaID_retrained, &b_recoB_l1_PFMvaID_retrained);
   fChain->SetBranchAddress("recoB_l2_PFMvaID_retrained", &recoB_l2_PFMvaID_retrained, &b_recoB_l2_PFMvaID_retrained);
   fChain->SetBranchAddress("recoB_l1_iso04_dca", &recoB_l1_iso04_dca, &b_recoB_l1_iso04_dca);
   fChain->SetBranchAddress("recoB_l2_iso04_dca", &recoB_l2_iso04_dca, &b_recoB_l2_iso04_dca);
   fChain->SetBranchAddress("recoB_b_iso04_dca", &recoB_b_iso04_dca, &b_recoB_b_iso04_dca);
   fChain->SetBranchAddress("recoB_k_iso04_dca", &recoB_k_iso04_dca, &b_recoB_k_iso04_dca);
   fChain->SetBranchAddress("recoB_k_svip3d", &recoB_k_svip3d, &b_recoB_k_svip3d);
   fChain->SetBranchAddress("recoB_k_svip3d_err", &recoB_k_svip3d_err, &b_recoB_k_svip3d_err);
   fChain->SetBranchAddress("recoB_l1_n_isotrk_dca", &recoB_l1_n_isotrk_dca, &b_recoB_l1_n_isotrk_dca);
   fChain->SetBranchAddress("recoB_l2_n_isotrk_dca", &recoB_l2_n_isotrk_dca, &b_recoB_l2_n_isotrk_dca);
   fChain->SetBranchAddress("recoB_k_n_isotrk_dca", &recoB_k_n_isotrk_dca, &b_recoB_k_n_isotrk_dca);
   fChain->SetBranchAddress("recoB_Idx", &recoB_Idx, &b_recoB_Idx);
   fChain->SetBranchAddress("recoB_l_xy_sig", &recoB_l_xy_sig, &b_recoB_l_xy_sig);
   fChain->SetBranchAddress("recoB_l1l2Dz", &recoB_l1l2Dz, &b_recoB_l1l2Dz);
   fChain->SetBranchAddress("recoB_lKDz", &recoB_lKDz, &b_recoB_lKDz);
   fChain->SetBranchAddress("recoB_l1l2Dr", &recoB_l1l2Dr, &b_recoB_l1l2Dr);
   fChain->SetBranchAddress("recoB_lKDr", &recoB_lKDr, &b_recoB_lKDr);
   fChain->SetBranchAddress("recoB_k_opp_l_mass", &recoB_k_opp_l_mass, &b_recoB_k_opp_l_mass);
   fChain->SetBranchAddress("recoB_p_assymetry", &recoB_p_assymetry, &b_recoB_p_assymetry);
   Notify();
}

Bool_t EventsBaseMC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EventsBaseMC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EventsBaseMC::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Events_cxx
