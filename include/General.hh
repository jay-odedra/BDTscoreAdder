// from TChain Events/
//////////////////////////////////////////////////////////

#ifndef General_hh
#define General_hh

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string>
#include "../include/GeneralBase.h"


class General : public GeneralBase {
public:
    General(TTree *tree=0);
    void Loop(std::string outname, std::string outdir);
    void Output(std::string outname, std::string outdir);
    void initVars();
    
    int           run_;
    int           luminosityBlock_;
    long long int         event_;
    int           bunchCrossing_;
    int           nBToKEE_;
    int           nElectron_;
    int           nGenPart_;
    float         genWeight_;
    int           nPSWeight_;
    int           nProbeTracks_;
    bool          HLT_DoubleEle10_eta1p22_mMax6_;
    bool          HLT_DoubleEle9p5_eta1p22_mMax6_;
    bool          HLT_DoubleEle9_eta1p22_mMax6_;
    bool          HLT_DoubleEle8p5_eta1p22_mMax6_;
    bool          HLT_DoubleEle8_eta1p22_mMax6_;
    bool          HLT_DoubleEle7p5_eta1p22_mMax6_;
    bool          HLT_DoubleEle7_eta1p22_mMax6_;
    bool          HLT_DoubleEle6p5_eta1p22_mMax6_;
    bool          HLT_DoubleEle6_eta1p22_mMax6_;
    bool          HLT_DoubleEle5p5_eta1p22_mMax6_;
    bool          HLT_DoubleEle5_eta1p22_mMax6_;
    bool          HLT_DoubleEle4p5_eta1p22_mMax6_;
    bool          HLT_DoubleEle4_eta1p22_mMax6_;
    bool          L1_DoubleEG11_er1p2_dR_Max0p6_;
    bool          L1_DoubleEG10p5_er1p2_dR_Max0p6_;
    bool          L1_DoubleEG10_er1p2_dR_Max0p6_;
    bool          L1_DoubleEG9p5_er1p2_dR_Max0p6_;
    bool          L1_DoubleEG9_er1p2_dR_Max0p7_;
    bool          L1_DoubleEG8p5_er1p2_dR_Max0p7_;
    bool          L1_DoubleEG8_er1p2_dR_Max0p7_;
    bool          L1_DoubleEG7p5_er1p2_dR_Max0p7_;
    bool          L1_DoubleEG7_er1p2_dR_Max0p8_;
    bool          L1_DoubleEG6p5_er1p2_dR_Max0p8_;
    bool          L1_DoubleEG6_er1p2_dR_Max0p8_;
    bool          L1_DoubleEG5p5_er1p2_dR_Max0p8_;
    bool          L1_DoubleEG5_er1p2_dR_Max0p9_;
    bool          L1_DoubleEG4p5_er1p2_dR_Max0p9_;
    bool          L1_DoubleEG4_er1p2_dR_Max0p9_;
    int           nTrigObj_;
    int           nOtherPV_;
    float           PV_ndof_;
    float           PV_x_;
    float           PV_y_;
    float           PV_z_;
    float           PV_chi2_;
    float           PV_score_;
    int            PV_npvs_;
    int            PV_npvsGood_;
    int           nSV_;
    float         genB_pdgId_;
    float         genE1_pdgId_;
    float         genE2_pdgId_;
    float         genK_pdgId_;
    float         genB_pt_;
    float         genE1_pt_;
    float         genE2_pt_;
    float         genK_pt_;
    float         genB_eta_;
    float         genE1_eta_;
    float         genE2_eta_;
    float         genK_eta_;
    float         genB_phi_;
    float         genE1_phi_;
    float         genE2_phi_;
    float         genK_phi_;
    float         genB_mass_;
    float         genE1_mass_;
    float         genE2_mass_;
    float         genK_mass_;
    float         genB_charge_;
    float         genE1_charge_;
    float         genE2_charge_;
    float         genK_charge_;
    int           nSkimBToKEE_;
    float         BToKEE_fit_pt_;
    float         BToKEE_fit_eta_;
    float         BToKEE_fit_phi_;
    float         BToKEE_D0_mass_LepToK_KToPi_;
    float         BToKEE_D0_mass_LepToPi_KToK_;
    float         BToKEE_fit_mass_;
    float         BToKEE_l_xy_;
    float         BToKEE_l_xy_unc_;
    float         BToKEE_fit_cos2D_;
    float         BToKEE_svprob_;
    float         BToKEE_fit_massErr_;
    float         BToKEE_b_iso04_;
    float         BToKEE_mll_fullfit_;
    float         BToKEE_vtx_x_;
    float         BToKEE_vtx_y_;
    float         BToKEE_vtx_z_;
    float         BToKEE_l1Idx_;
    float         BToKEE_l2Idx_;
    float         BToKEE_kIdx_;
    float         BToKEE_fit_k_pt_;
    float         BToKEE_fit_k_eta_;
    float         BToKEE_fit_k_phi_;
    float         BToKEE_fit_l1_pt_;
    float         BToKEE_fit_l1_eta_;
    float         BToKEE_fit_l1_phi_;
    float         BToKEE_fit_l2_pt_;
    float         BToKEE_fit_l2_eta_;
    float         BToKEE_fit_l2_phi_;
    float         BToKEE_l1_iso04_;
    float         BToKEE_l2_iso04_;
    float         BToKEE_l1_isPF_;
    float         BToKEE_l2_isPF_;
    float         BToKEE_k_iso04_;
    float         BToKEE_l1_isPFoverlap_;
    float         BToKEE_l2_isPFoverlap_;
    float         BToKEE_l1_LooseId_;
    float         BToKEE_l2_LooseId_;
    float         BToKEE_l1_MediumId_;
    float         BToKEE_l2_MediumId_;
    float         BToKEE_l1_TightId_;
    float         BToKEE_l2_TightId_;
    float         BToKEE_l1_ConvVeto_;
    float         BToKEE_l2_ConvVeto_;
    float         BToKEE_l1_PFMvaID_Fall17_;
    float         BToKEE_l2_PFMvaID_Fall17_;
    float         BToKEE_l1_PFMvaID_retrained_;
    float         BToKEE_l2_PFMvaID_retrained_;
    float         BToKEE_l1_iso04_dca_;
    float         BToKEE_l2_iso04_dca_;
    float         BToKEE_b_iso04_dca_;
    float         BToKEE_k_iso04_dca_;
    float         BToKEE_k_svip3d_;
    float         BToKEE_k_svip3d_err_;
    float         BToKEE_l1_n_isotrk_dca_;
    float         BToKEE_l2_n_isotrk_dca_;
    float         BToKEE_k_n_isotrk_dca_;
    float         BToKEE_l1_vx_;
    float         BToKEE_l1_vy_;
    float         BToKEE_l1_vz_;
    float         BToKEE_l1_charge_;
    float         BToKEE_l2_vx_;
    float         BToKEE_l2_vy_;
    float         BToKEE_l2_vz_;
    float         BToKEE_l2_charge_;
    float         BToKEE_k_vx_;
    float         BToKEE_k_vy_;
    float         BToKEE_k_vz_;
    float         BToKEE_k_dca_sig_;
    float         BToKEE_k_dz_;
    float         BToKEE_kMu_matched_;
    float         BToKEE_l_xy_sig_;
    float         BToKEE_l1l2Dz_;
    float         BToKEE_lKDz_;
    float         BToKEE_l1l2Dr_;
    float         BToKEE_lKDr_;
    float         BToKEE_kl_massKPi_;
    float         BToKEE_p_assymetry_;
    float         BDTSCORE_1_;

    TFile* outFile_=nullptr;
    TTree* outTree_=nullptr;



};

#endif
