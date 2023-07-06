#define Events_cxx
#include "Events.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TMath.h>
#include <TTree.h>
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector2.h>               // 2D vector class
#include <TMath.h>                  // ROOT math library
#include <TLeaf.h>
#include <list>
#include <TCanvas.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <typeinfo>
#include <TGraph.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TNtuple.h> 
#include <TH1.h>
#include <ROOT/RDataFrame.hxx>
#include "ROOT/RVec.hxx"
#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include "../FastForest/include/fastforest.h"
#include <TTreeReaderArray.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>

void LoadBDT(fastforest::FastForest& bdt,std::string modellocation,std::vector<std::string> features){
    bdt = fastforest::load_txt(modellocation.c_str(), features);
}

float evaluatemodel(fastforest::FastForest& bdt,std::vector<float> input){
    float score = bdt(input.data());
    return score;
}
void timing( int nentries, int jentry, auto start ) {
  auto now = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = now-start;
  std::chrono::duration<double> predicted_duration = elapsed_seconds * (nentries*1.)/(jentry*1.);
  std::chrono::system_clock::time_point end = start + std::chrono::seconds((int)predicted_duration.count());
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
  std::string tmp = std::ctime(&end_time); 
  tmp.resize(tmp.length()-1); // remove trailing \n
  std::cout << "Processed event " << jentry << " out of " << nentries << " (" << int(100.*jentry/nentries) << "%). "
	    << "Time elapsed: " << elapsed_seconds.count() << " s. "
	    << "Time remaining: " << (predicted_duration.count() - elapsed_seconds.count() ) << " s. "
	    << "Predicted finish time: " << tmp
	    << std::endl;
}

void Events::Loop(std::string outname, std::string outdir)
{
//   In a ROOT session, you can do:
//      root> .L Events.C
//      root> Events t
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

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   TString outputfilename(std::string(outdir)+"/"+std::string(outname)+"_.root");

   TFile* outFile_ = new TFile(outputfilename,"RECREATE");
   TTree* outTree_ = new TTree("tree","tree");
   int          run_;
   int          luminosityBlock_;
   long long int       event_;
   int          bunchCrossing_;
   int          nBToKEE_;
   int          nElectron_;
   int          nProbeTracks_;
   bool         HLT_DoubleEle10_eta1p22_mMax6_;
   bool         HLT_DoubleEle9p5_eta1p22_mMax6_;
   bool         HLT_DoubleEle9_eta1p22_mMax6_;
   bool         HLT_DoubleEle8p5_eta1p22_mMax6_;
   bool         HLT_DoubleEle8_eta1p22_mMax6_;
   bool         HLT_DoubleEle7p5_eta1p22_mMax6_;
   bool         HLT_DoubleEle7_eta1p22_mMax6_;
   bool         HLT_DoubleEle6p5_eta1p22_mMax6_;
   bool         HLT_DoubleEle6_eta1p22_mMax6_;
   bool         HLT_DoubleEle5p5_eta1p22_mMax6_;
   bool         HLT_DoubleEle5_eta1p22_mMax6_;
   bool         HLT_DoubleEle4p5_eta1p22_mMax6_;
   bool         HLT_DoubleEle4_eta1p22_mMax6_;
   bool         L1_DoubleEG11_er1p2_dR_Max0p6_;
   bool         L1_DoubleEG10p5_er1p2_dR_Max0p6_;
   bool         L1_DoubleEG10_er1p2_dR_Max0p6_;
   bool         L1_DoubleEG9p5_er1p2_dR_Max0p6_;
   bool         L1_DoubleEG9_er1p2_dR_Max0p7_;
   bool         L1_DoubleEG8p5_er1p2_dR_Max0p7_;
   bool         L1_DoubleEG8_er1p2_dR_Max0p7_;
   bool         L1_DoubleEG7p5_er1p2_dR_Max0p7_;
   bool         L1_DoubleEG7_er1p2_dR_Max0p8_;
   bool         L1_DoubleEG6p5_er1p2_dR_Max0p8_;
   bool         L1_DoubleEG6_er1p2_dR_Max0p8_;
   bool         L1_DoubleEG5p5_er1p2_dR_Max0p8_;
   bool         L1_DoubleEG5_er1p2_dR_Max0p9_;
   bool         L1_DoubleEG4p5_er1p2_dR_Max0p9_;
   bool         L1_DoubleEG4_er1p2_dR_Max0p9_;
   int          nTrigObj_;
   int          nOtherPV_;
   float         PV_ndof_;
   float         PV_x_;
   float         PV_y_;
   float         PV_z_;
   float         PV_chi2_;
   float         PV_score_;
   int          PV_npvs_;
   int          PV_npvsGood_;
   int          nSV_;
   int          nSkimBToKEE_;


   
   float        SkimBToKEE_fit_pt_;
   float        SkimBToKEE_fit_eta_;
   float        SkimBToKEE_fit_phi_;
   float        SkimBToKEE_D0_mass_LepToK_KToPi_;
   float        SkimBToKEE_D0_mass_LepToPi_KToK_;
   float        SkimBToKEE_fit_mass_;
   float        SkimBToKEE_l_xy_;
   float        SkimBToKEE_l_xy_unc_;
   float        SkimBToKEE_fit_cos2D_;
   float        SkimBToKEE_svprob_;
   float        SkimBToKEE_fit_massErr_;
   float        SkimBToKEE_b_iso04_;
   float        SkimBToKEE_mll_fullfit_;
   float        SkimBToKEE_vtx_x_;
   float        SkimBToKEE_vtx_y_;
   float        SkimBToKEE_vtx_z_;
   float        SkimBToKEE_l1Idx_;
   float        SkimBToKEE_l2Idx_;
   float        SkimBToKEE_kIdx_;
   float        SkimBToKEE_fit_k_pt_;
   float        SkimBToKEE_fit_k_eta_;
   float        SkimBToKEE_fit_k_phi_;
   float        SkimBToKEE_fit_l1_pt_;
   float        SkimBToKEE_fit_l1_eta_;
   float        SkimBToKEE_fit_l1_phi_;
   float        SkimBToKEE_fit_l2_pt_;
   float        SkimBToKEE_fit_l2_eta_;
   float        SkimBToKEE_fit_l2_phi_;
   float        SkimBToKEE_l1_iso04_;
   float        SkimBToKEE_l2_iso04_;
   float        SkimBToKEE_l1_isPF_;
   float        SkimBToKEE_l2_isPF_;
   float        SkimBToKEE_k_iso04_;
   float        SkimBToKEE_l1_isPFoverlap_;
   float        SkimBToKEE_l2_isPFoverlap_;
   float        SkimBToKEE_l1_LooseId_;
   float        SkimBToKEE_l2_LooseId_;
   float        SkimBToKEE_l1_MediumId_;
   float        SkimBToKEE_l2_MediumId_;
   float        SkimBToKEE_l1_TightId_;
   float        SkimBToKEE_l2_TightId_;
   float        SkimBToKEE_l1_ConvVeto_;
   float        SkimBToKEE_l2_ConvVeto_;
   float        SkimBToKEE_l1_PFMvaID_Fall17_;
   float        SkimBToKEE_l2_PFMvaID_Fall17_;
   float        SkimBToKEE_l1_PFMvaID_retrained_;
   float        SkimBToKEE_l2_PFMvaID_retrained_;
   float        SkimBToKEE_l1_iso04_dca_;
   float        SkimBToKEE_l2_iso04_dca_;
   float        SkimBToKEE_b_iso04_dca_;
   float        SkimBToKEE_k_iso04_dca_;
   float        SkimBToKEE_k_svip3d_;
   float        SkimBToKEE_k_svip3d_err_;
   float        SkimBToKEE_l1_n_isotrk_dca_;
   float        SkimBToKEE_l2_n_isotrk_dca_;
   float        SkimBToKEE_k_n_isotrk_dca_;
   float        SkimBToKEE_l1_vx_;
   float        SkimBToKEE_l1_vy_;
   float        SkimBToKEE_l1_vz_;
   float        SkimBToKEE_l1_charge_;
   float        SkimBToKEE_l2_vx_;
   float        SkimBToKEE_l2_vy_;
   float        SkimBToKEE_l2_vz_;
   float        SkimBToKEE_l2_charge_;
   float        SkimBToKEE_k_vx_;
   float        SkimBToKEE_k_vy_;
   float        SkimBToKEE_k_vz_;
   float        SkimBToKEE_k_dca_sig_;
   float        SkimBToKEE_k_dz_;
   float        SkimBToKEE_kMu_matched_;
   float        SkimBToKEE_k_charge_;
   float        SkimBToKEE_l_xy_sig_;
   float        SkimBToKEE_l1l2Dz_;
   float        SkimBToKEE_lKDz_;
   float        SkimBToKEE_l1l2Dr_;
   float        SkimBToKEE_lKDr_;
   float        SkimBToKEE_kl_massKPi_;
   float        SkimBToKEE_p_assymetry_;
   float         BDTSCORE_1;        

   outTree_->Branch("run", &run_);;
   outTree_->Branch("luminosityBlock", &luminosityBlock_);;
   outTree_->Branch(" event", &       event_);;
   outTree_->Branch("bunchCrossing", &bunchCrossing_);;
   outTree_->Branch("nBToKEE", &nBToKEE_);;
   outTree_->Branch("nElectron", &nElectron_);;
   outTree_->Branch("nProbeTracks", &nProbeTracks_);;
   outTree_->Branch("HLT_DoubleEle10_eta1p22_mMax6", &HLT_DoubleEle10_eta1p22_mMax6_);
   outTree_->Branch("HLT_DoubleEle9p5_eta1p22_mMax6", &HLT_DoubleEle9p5_eta1p22_mMax6_);
   outTree_->Branch("HLT_DoubleEle9_eta1p22_mMax6", &HLT_DoubleEle9_eta1p22_mMax6_);
   outTree_->Branch("HLT_DoubleEle8p5_eta1p22_mMax6", &HLT_DoubleEle8p5_eta1p22_mMax6_);
   outTree_->Branch("HLT_DoubleEle8_eta1p22_mMax6", &HLT_DoubleEle8_eta1p22_mMax6_);
   outTree_->Branch("HLT_DoubleEle7p5_eta1p22_mMax6", &HLT_DoubleEle7p5_eta1p22_mMax6_);
   outTree_->Branch("HLT_DoubleEle7_eta1p22_mMax6", &HLT_DoubleEle7_eta1p22_mMax6_);
   outTree_->Branch("HLT_DoubleEle6p5_eta1p22_mMax6", &HLT_DoubleEle6p5_eta1p22_mMax6_);
   outTree_->Branch("HLT_DoubleEle6_eta1p22_mMax6", &HLT_DoubleEle6_eta1p22_mMax6_);
   outTree_->Branch("HLT_DoubleEle5p5_eta1p22_mMax6", &HLT_DoubleEle5p5_eta1p22_mMax6_);
   outTree_->Branch("HLT_DoubleEle5_eta1p22_mMax6", &HLT_DoubleEle5_eta1p22_mMax6_);
   outTree_->Branch("HLT_DoubleEle4p5_eta1p22_mMax6", &HLT_DoubleEle4p5_eta1p22_mMax6_);
   outTree_->Branch("HLT_DoubleEle4_eta1p22_mMax6", &HLT_DoubleEle4_eta1p22_mMax6_);
   outTree_->Branch("L1_DoubleEG11_er1p2_dR_Max0p6", &L1_DoubleEG11_er1p2_dR_Max0p6_);
   outTree_->Branch("L1_DoubleEG10p5_er1p2_dR_Max0p6", &L1_DoubleEG10p5_er1p2_dR_Max0p6);
   outTree_->Branch("L1_DoubleEG10_er1p2_dR_Max0p6", &L1_DoubleEG10_er1p2_dR_Max0p6_);
   outTree_->Branch("L1_DoubleEG9p5_er1p2_dR_Max0p6", &L1_DoubleEG9p5_er1p2_dR_Max0p6_);
   outTree_->Branch("L1_DoubleEG9_er1p2_dR_Max0p7", &L1_DoubleEG9_er1p2_dR_Max0p7_);
   outTree_->Branch("L1_DoubleEG8p5_er1p2_dR_Max0p7", &L1_DoubleEG8p5_er1p2_dR_Max0p7_);
   outTree_->Branch("L1_DoubleEG8_er1p2_dR_Max0p7", &L1_DoubleEG8_er1p2_dR_Max0p7_);
   outTree_->Branch("L1_DoubleEG7p5_er1p2_dR_Max0p7", &L1_DoubleEG7p5_er1p2_dR_Max0p7_);
   outTree_->Branch("L1_DoubleEG7_er1p2_dR_Max0p8", &L1_DoubleEG7_er1p2_dR_Max0p8_);
   outTree_->Branch("L1_DoubleEG6p5_er1p2_dR_Max0p8", &L1_DoubleEG6p5_er1p2_dR_Max0p8_);
   outTree_->Branch("L1_DoubleEG6_er1p2_dR_Max0p8", &L1_DoubleEG6_er1p2_dR_Max0p8_);
   outTree_->Branch("L1_DoubleEG5p5_er1p2_dR_Max0p8", &L1_DoubleEG5p5_er1p2_dR_Max0p8_);
   outTree_->Branch("L1_DoubleEG5_er1p2_dR_Max0p9", &L1_DoubleEG5_er1p2_dR_Max0p9_);
   outTree_->Branch("L1_DoubleEG4p5_er1p2_dR_Max0p9", &L1_DoubleEG4p5_er1p2_dR_Max0p9_);
   outTree_->Branch("L1_DoubleEG4_er1p2_dR_Max0p9", &L1_DoubleEG4_er1p2_dR_Max0p9_);
   outTree_->Branch("nTrigObj", &nTrigObj_);
   outTree_->Branch("nOtherPV", &nOtherPV_);
   outTree_->Branch("  PV_ndof", &  PV_ndof_);
   outTree_->Branch("  PV_x", &  PV_x_);
   outTree_->Branch("  PV_y", &  PV_y_);
   outTree_->Branch("  PV_z", &  PV_z_);
   outTree_->Branch("  PV_chi2", &  PV_chi2_);
   outTree_->Branch("  PV_score", &  PV_score_);
   outTree_->Branch(" PV_npvs", & PV_npvs_);
   outTree_->Branch(" PV_npvsGood", & PV_npvsGood_);
   outTree_->Branch("nSV", &nSV_);
   outTree_->Branch("nSkimBToKEE", &nSkimBToKEE_);



   outTree_->Branch("SkimBToKEE_fit_pt", &SkimBToKEE_fit_pt_);
   outTree_->Branch("SkimBToKEE_fit_eta", &SkimBToKEE_fit_eta_);
   outTree_->Branch("SkimBToKEE_fit_phi", &SkimBToKEE_fit_phi_);
   outTree_->Branch("SkimBToKEE_D0_mass_LepToK_KToPi", &SkimBToKEE_D0_mass_LepToK_KToPi_);
   outTree_->Branch("SkimBToKEE_D0_mass_LepToPi_KToK", &SkimBToKEE_D0_mass_LepToPi_KToK_);
   outTree_->Branch("SkimBToKEE_fit_mass", &SkimBToKEE_fit_mass_);
   outTree_->Branch("SkimBToKEE_l_xy", &SkimBToKEE_l_xy_);
   outTree_->Branch("SkimBToKEE_l_xy_unc", &SkimBToKEE_l_xy_unc_);
   outTree_->Branch("SkimBToKEE_fit_cos2D", &SkimBToKEE_fit_cos2D_);
   outTree_->Branch("SkimBToKEE_svprob", &SkimBToKEE_svprob_);
   outTree_->Branch("SkimBToKEE_fit_massErr", &SkimBToKEE_fit_massErr_);
   outTree_->Branch("SkimBToKEE_b_iso04", &SkimBToKEE_b_iso04_);
   outTree_->Branch("SkimBToKEE_mll_fullfit", &SkimBToKEE_mll_fullfit_);
   outTree_->Branch("SkimBToKEE_vtx_x", &SkimBToKEE_vtx_x_);
   outTree_->Branch("SkimBToKEE_vtx_y", &SkimBToKEE_vtx_y_);
   outTree_->Branch("SkimBToKEE_vtx_z", &SkimBToKEE_vtx_z_);
   outTree_->Branch("SkimBToKEE_l1Idx", &SkimBToKEE_l1Idx_);
   outTree_->Branch("SkimBToKEE_l2Idx", &SkimBToKEE_l2Idx_);
   outTree_->Branch("SkimBToKEE_kIdx", &SkimBToKEE_kIdx_);
   outTree_->Branch("SkimBToKEE_fit_k_pt", &SkimBToKEE_fit_k_pt_);
   outTree_->Branch("SkimBToKEE_fit_k_eta", &SkimBToKEE_fit_k_eta_);
   outTree_->Branch("SkimBToKEE_fit_k_phi", &SkimBToKEE_fit_k_phi_);
   outTree_->Branch("SkimBToKEE_fit_l1_pt", &SkimBToKEE_fit_l1_pt_);
   outTree_->Branch("SkimBToKEE_fit_l1_eta", &SkimBToKEE_fit_l1_eta_);
   outTree_->Branch("SkimBToKEE_fit_l1_phi", &SkimBToKEE_fit_l1_phi_);
   outTree_->Branch("SkimBToKEE_fit_l2_pt", &SkimBToKEE_fit_l2_pt_);
   outTree_->Branch("SkimBToKEE_fit_l2_eta", &SkimBToKEE_fit_l2_eta_);
   outTree_->Branch("SkimBToKEE_fit_l2_phi", &SkimBToKEE_fit_l2_phi_);
   outTree_->Branch("SkimBToKEE_l1_iso04", &SkimBToKEE_l1_iso04_);
   outTree_->Branch("SkimBToKEE_l2_iso04", &SkimBToKEE_l2_iso04_);
   outTree_->Branch("SkimBToKEE_l1_isPF", &SkimBToKEE_l1_isPF_);
   outTree_->Branch("SkimBToKEE_l2_isPF", &SkimBToKEE_l2_isPF_);
   outTree_->Branch("SkimBToKEE_k_iso04", &SkimBToKEE_k_iso04_);
   outTree_->Branch("SkimBToKEE_l1_isPFoverlap", &SkimBToKEE_l1_isPFoverlap_);
   outTree_->Branch("SkimBToKEE_l2_isPFoverlap", &SkimBToKEE_l2_isPFoverlap_);
   outTree_->Branch("SkimBToKEE_l1_LooseId", &SkimBToKEE_l1_LooseId_);
   outTree_->Branch("SkimBToKEE_l2_LooseId", &SkimBToKEE_l2_LooseId_);
   outTree_->Branch("SkimBToKEE_l1_MediumId", &SkimBToKEE_l1_MediumId_);
   outTree_->Branch("SkimBToKEE_l2_MediumId", &SkimBToKEE_l2_MediumId_);
   outTree_->Branch("SkimBToKEE_l1_TightId", &SkimBToKEE_l1_TightId_);
   outTree_->Branch("SkimBToKEE_l2_TightId", &SkimBToKEE_l2_TightId_);
   outTree_->Branch("SkimBToKEE_l1_ConvVeto", &SkimBToKEE_l1_ConvVeto_);
   outTree_->Branch("SkimBToKEE_l2_ConvVeto", &SkimBToKEE_l2_ConvVeto_);
   outTree_->Branch("SkimBToKEE_l1_PFMvaID_Fall17", &SkimBToKEE_l1_PFMvaID_Fall17_);
   outTree_->Branch("SkimBToKEE_l2_PFMvaID_Fall17", &SkimBToKEE_l2_PFMvaID_Fall17_);
   outTree_->Branch("SkimBToKEE_l1_PFMvaID_retrained", &SkimBToKEE_l1_PFMvaID_retrained_);
   outTree_->Branch("SkimBToKEE_l2_PFMvaID_retrained", &SkimBToKEE_l2_PFMvaID_retrained_);
   outTree_->Branch("SkimBToKEE_l1_iso04_dca", &SkimBToKEE_l1_iso04_dca_);
   outTree_->Branch("SkimBToKEE_l2_iso04_dca", &SkimBToKEE_l2_iso04_dca_);
   outTree_->Branch("SkimBToKEE_b_iso04_dca", &SkimBToKEE_b_iso04_dca_);
   outTree_->Branch("SkimBToKEE_k_iso04_dca", &SkimBToKEE_k_iso04_dca_);
   outTree_->Branch("SkimBToKEE_k_svip3d", &SkimBToKEE_k_svip3d_);
   outTree_->Branch("SkimBToKEE_k_svip3d_err", &SkimBToKEE_k_svip3d_err_);
   outTree_->Branch("SkimBToKEE_l1_n_isotrk_dca", &SkimBToKEE_l1_n_isotrk_dca_);
   outTree_->Branch("SkimBToKEE_l2_n_isotrk_dca", &SkimBToKEE_l2_n_isotrk_dca_);
   outTree_->Branch("SkimBToKEE_k_n_isotrk_dca", &SkimBToKEE_k_n_isotrk_dca_);
   outTree_->Branch("SkimBToKEE_l1_vx", &SkimBToKEE_l1_vx_);
   outTree_->Branch("SkimBToKEE_l1_vy", &SkimBToKEE_l1_vy_);
   outTree_->Branch("SkimBToKEE_l1_vz", &SkimBToKEE_l1_vz_);
   outTree_->Branch("SkimBToKEE_l1_charge", &SkimBToKEE_l1_charge_);
   outTree_->Branch("SkimBToKEE_l2_vx", &SkimBToKEE_l2_vx_);
   outTree_->Branch("SkimBToKEE_l2_vy", &SkimBToKEE_l2_vy_);
   outTree_->Branch("SkimBToKEE_l2_vz", &SkimBToKEE_l2_vz_);
   outTree_->Branch("SkimBToKEE_l2_charge", &SkimBToKEE_l2_charge_);
   outTree_->Branch("SkimBToKEE_k_vx", &SkimBToKEE_k_vx_);
   outTree_->Branch("SkimBToKEE_k_vy", &SkimBToKEE_k_vy_);
   outTree_->Branch("SkimBToKEE_k_vz", &SkimBToKEE_k_vz_);
   outTree_->Branch("SkimBToKEE_k_dca_sig", &SkimBToKEE_k_dca_sig_);
   outTree_->Branch("SkimBToKEE_k_dz", &SkimBToKEE_k_dz_);
   outTree_->Branch("SkimBToKEE_kMu_matched", &SkimBToKEE_kMu_matched_);
   outTree_->Branch("SkimBToKEE_k_charge", &SkimBToKEE_k_charge_);
   outTree_->Branch("SkimBToKEE_l_xy_sig", &SkimBToKEE_l_xy_sig_);
   outTree_->Branch("SkimBToKEE_l1l2Dz", &SkimBToKEE_l1l2Dz_);
   outTree_->Branch("SkimBToKEE_lKDz", &SkimBToKEE_lKDz_);
   outTree_->Branch("SkimBToKEE_l1l2Dr", &SkimBToKEE_l1l2Dr_);
   outTree_->Branch("SkimBToKEE_lKDr", &SkimBToKEE_lKDr_);
   outTree_->Branch("SkimBToKEE_kl_massKPi", &SkimBToKEE_kl_massKPi_);
   outTree_->Branch("SkimBToKEE_p_assymetry", &SkimBToKEE_p_assymetry_);

   outTree_->Branch("BDTSCORE_1", &BDTSCORE_1);




   //initialize bdt
   fastforest::FastForest bdt1;
   std::string modellocation = "/eos/user/j/jodedra/fastforesttester/CMSSW_13_1_0/src/XGB_89_1.txt";
   std::vector<std::string> features{"SkimBToKEE_fit_pt","SkimBToKEE_fit_eta","SkimBToKEE_fit_phi","SkimBToKEE_D0_mass_LepToK_KToPi","SkimBToKEE_D0_mass_LepToPi_KToK","SkimBToKEE_fit_cos2D","SkimBToKEE_svprob","SkimBToKEE_b_iso04","SkimBToKEE_fit_k_pt","SkimBToKEE_fit_k_eta","SkimBToKEE_fit_k_phi","SkimBToKEE_fit_l1_pt","SkimBToKEE_fit_l1_eta","SkimBToKEE_fit_l1_phi","SkimBToKEE_fit_l2_pt","SkimBToKEE_fit_l2_eta","SkimBToKEE_fit_l2_phi","SkimBToKEE_l1_iso04","SkimBToKEE_l2_iso04","SkimBToKEE_k_iso04","SkimBToKEE_l1_PFMvaID_Fall17","SkimBToKEE_l2_PFMvaID_Fall17","SkimBToKEE_l1_PFMvaID_retrained","SkimBToKEE_l2_PFMvaID_retrained","SkimBToKEE_l1_iso04_dca","SkimBToKEE_l2_iso04_dca","SkimBToKEE_b_iso04_dca","SkimBToKEE_k_iso04_dca","SkimBToKEE_l1_n_isotrk_dca","SkimBToKEE_l2_n_isotrk_dca","SkimBToKEE_k_n_isotrk_dca","SkimBToKEE_l_xy_sig","SkimBToKEE_l1l2Dz","SkimBToKEE_lKDz","SkimBToKEE_l1l2Dr","SkimBToKEE_p_assymetry","SkimBToKEE_k_svip3d_sig","SkimBToKEE_fit_pt_over_mass","SkimBToKEE_fit_l1_pt_over_mass","SkimBToKEE_fit_l2_pt_over_mass","SkimBToKEE_fit_k_pt_over_mass"};
   LoadBDT(bdt1,modellocation,features); 
   int interval = 1000;
   auto start = std::chrono::system_clock::now();

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      if (jentry%interval==0) { timing(nentries,jentry,start); }

      run_ = run ;
      luminosityBlock_ = luminosityBlock ;
      event_ =        event ;
      bunchCrossing_ = bunchCrossing ;
      nBToKEE_ = nBToKEE ;
      nElectron_ = nElectron ;
      nProbeTracks_ = nProbeTracks ;
      HLT_DoubleEle10_eta1p22_mMax6_ = HLT_DoubleEle10_eta1p22_mMax6 ;
      HLT_DoubleEle9p5_eta1p22_mMax6_ = HLT_DoubleEle9p5_eta1p22_mMax6 ;
      HLT_DoubleEle9_eta1p22_mMax6_ = HLT_DoubleEle9_eta1p22_mMax6 ;
      HLT_DoubleEle8p5_eta1p22_mMax6_ = HLT_DoubleEle8p5_eta1p22_mMax6 ;
      HLT_DoubleEle8_eta1p22_mMax6_ = HLT_DoubleEle8_eta1p22_mMax6 ;
      HLT_DoubleEle7p5_eta1p22_mMax6_ = HLT_DoubleEle7p5_eta1p22_mMax6 ;
      HLT_DoubleEle7_eta1p22_mMax6_ = HLT_DoubleEle7_eta1p22_mMax6 ;
      HLT_DoubleEle6p5_eta1p22_mMax6_ = HLT_DoubleEle6p5_eta1p22_mMax6 ;
      HLT_DoubleEle6_eta1p22_mMax6_ = HLT_DoubleEle6_eta1p22_mMax6 ;
      HLT_DoubleEle5p5_eta1p22_mMax6_ = HLT_DoubleEle5p5_eta1p22_mMax6 ;
      HLT_DoubleEle5_eta1p22_mMax6_ = HLT_DoubleEle5_eta1p22_mMax6 ;
      HLT_DoubleEle4p5_eta1p22_mMax6_ = HLT_DoubleEle4p5_eta1p22_mMax6 ;
      HLT_DoubleEle4_eta1p22_mMax6_ = HLT_DoubleEle4_eta1p22_mMax6 ;
      L1_DoubleEG11_er1p2_dR_Max0p6_ = L1_DoubleEG11_er1p2_dR_Max0p6 ;
      L1_DoubleEG10p5_er1p2_dR_Max0p6 = L1_DoubleEG10p5_er1p2_dR_Max0p6 ;
      L1_DoubleEG10_er1p2_dR_Max0p6_ = L1_DoubleEG10_er1p2_dR_Max0p6 ;
      L1_DoubleEG9p5_er1p2_dR_Max0p6_ = L1_DoubleEG9p5_er1p2_dR_Max0p6 ;
      L1_DoubleEG9_er1p2_dR_Max0p7_ = L1_DoubleEG9_er1p2_dR_Max0p7 ;
      L1_DoubleEG8p5_er1p2_dR_Max0p7_ = L1_DoubleEG8p5_er1p2_dR_Max0p7 ;
      L1_DoubleEG8_er1p2_dR_Max0p7_ = L1_DoubleEG8_er1p2_dR_Max0p7 ;
      L1_DoubleEG7p5_er1p2_dR_Max0p7_ = L1_DoubleEG7p5_er1p2_dR_Max0p7 ;
      L1_DoubleEG7_er1p2_dR_Max0p8_ = L1_DoubleEG7_er1p2_dR_Max0p8 ;
      L1_DoubleEG6p5_er1p2_dR_Max0p8_ = L1_DoubleEG6p5_er1p2_dR_Max0p8 ;
      L1_DoubleEG6_er1p2_dR_Max0p8_ = L1_DoubleEG6_er1p2_dR_Max0p8 ;
      L1_DoubleEG5p5_er1p2_dR_Max0p8_ = L1_DoubleEG5p5_er1p2_dR_Max0p8 ;
      L1_DoubleEG5_er1p2_dR_Max0p9_ = L1_DoubleEG5_er1p2_dR_Max0p9 ;
      L1_DoubleEG4p5_er1p2_dR_Max0p9_ = L1_DoubleEG4p5_er1p2_dR_Max0p9 ;
      L1_DoubleEG4_er1p2_dR_Max0p9_ = L1_DoubleEG4_er1p2_dR_Max0p9 ;
      nTrigObj_ = nTrigObj ;
      nOtherPV_ = nOtherPV ;
      PV_ndof_ =   PV_ndof ;
      PV_x_ =   PV_x ;
      PV_y_ =   PV_y ;
      PV_z_ =   PV_z ;
      PV_chi2_ =   PV_chi2 ;
      PV_score_ =   PV_score ;
      PV_npvs_ =  PV_npvs ;
      PV_npvsGood_ =  PV_npvsGood ;
      nSV_ = nSV ;
      nSkimBToKEE_ = nSkimBToKEE ;
      for(int i=0;i<nSkimBToKEE_;i++){
         SkimBToKEE_fit_pt_ = SkimBToKEE_fit_pt[i];
         SkimBToKEE_fit_eta_ = SkimBToKEE_fit_eta[i];
         SkimBToKEE_fit_phi_ = SkimBToKEE_fit_phi[i];
         SkimBToKEE_D0_mass_LepToK_KToPi_ = SkimBToKEE_D0_mass_LepToK_KToPi[i];
         SkimBToKEE_D0_mass_LepToPi_KToK_ = SkimBToKEE_D0_mass_LepToPi_KToK[i];
         SkimBToKEE_fit_mass_ = SkimBToKEE_fit_mass[i];
         SkimBToKEE_l_xy_ = SkimBToKEE_l_xy[i];
         SkimBToKEE_l_xy_unc_ = SkimBToKEE_l_xy_unc[i];
         SkimBToKEE_fit_cos2D_ = SkimBToKEE_fit_cos2D[i];
         SkimBToKEE_svprob_ = SkimBToKEE_svprob[i];
         SkimBToKEE_fit_massErr_ = SkimBToKEE_fit_massErr[i];
         SkimBToKEE_b_iso04_ = SkimBToKEE_b_iso04[i];
         SkimBToKEE_mll_fullfit_ = SkimBToKEE_mll_fullfit[i];
         SkimBToKEE_vtx_x_ = SkimBToKEE_vtx_x[i];
         SkimBToKEE_vtx_y_ = SkimBToKEE_vtx_y[i];
         SkimBToKEE_vtx_z_ = SkimBToKEE_vtx_z[i];
         SkimBToKEE_l1Idx_ = SkimBToKEE_l1Idx[i];
         SkimBToKEE_l2Idx_ = SkimBToKEE_l2Idx[i];
         SkimBToKEE_kIdx_ = SkimBToKEE_kIdx[i];
         SkimBToKEE_fit_k_pt_ = SkimBToKEE_fit_k_pt[i];
         SkimBToKEE_fit_k_eta_ = SkimBToKEE_fit_k_eta[i];
         SkimBToKEE_fit_k_phi_ = SkimBToKEE_fit_k_phi[i];
         SkimBToKEE_fit_l1_pt_ = SkimBToKEE_fit_l1_pt[i];
         SkimBToKEE_fit_l1_eta_ = SkimBToKEE_fit_l1_eta[i];
         SkimBToKEE_fit_l1_phi_ = SkimBToKEE_fit_l1_phi[i];
         SkimBToKEE_fit_l2_pt_ = SkimBToKEE_fit_l2_pt[i];
         SkimBToKEE_fit_l2_eta_ = SkimBToKEE_fit_l2_eta[i];
         SkimBToKEE_fit_l2_phi_ = SkimBToKEE_fit_l2_phi[i];
         SkimBToKEE_l1_iso04_ = SkimBToKEE_l1_iso04[i];
         SkimBToKEE_l2_iso04_ = SkimBToKEE_l2_iso04[i];
         SkimBToKEE_l1_isPF_ = SkimBToKEE_l1_isPF[i];
         SkimBToKEE_l2_isPF_ = SkimBToKEE_l2_isPF[i];
         SkimBToKEE_k_iso04_ = SkimBToKEE_k_iso04[i];
         SkimBToKEE_l1_isPFoverlap_ = SkimBToKEE_l1_isPFoverlap[i];
         SkimBToKEE_l2_isPFoverlap_ = SkimBToKEE_l2_isPFoverlap[i];
         SkimBToKEE_l1_LooseId_ = SkimBToKEE_l1_LooseId[i];
         SkimBToKEE_l2_LooseId_ = SkimBToKEE_l2_LooseId[i];
         SkimBToKEE_l1_MediumId_ = SkimBToKEE_l1_MediumId[i];
         SkimBToKEE_l2_MediumId_ = SkimBToKEE_l2_MediumId[i];
         SkimBToKEE_l1_TightId_ = SkimBToKEE_l1_TightId[i];
         SkimBToKEE_l2_TightId_ = SkimBToKEE_l2_TightId[i];
         SkimBToKEE_l1_ConvVeto_ = SkimBToKEE_l1_ConvVeto[i];
         SkimBToKEE_l2_ConvVeto_ = SkimBToKEE_l2_ConvVeto[i];
         SkimBToKEE_l1_PFMvaID_Fall17_ = SkimBToKEE_l1_PFMvaID_Fall17[i];
         SkimBToKEE_l2_PFMvaID_Fall17_ = SkimBToKEE_l2_PFMvaID_Fall17[i];
         SkimBToKEE_l1_PFMvaID_retrained_ = SkimBToKEE_l1_PFMvaID_retrained[i];
         SkimBToKEE_l2_PFMvaID_retrained_ = SkimBToKEE_l2_PFMvaID_retrained[i];
         SkimBToKEE_l1_iso04_dca_ = SkimBToKEE_l1_iso04_dca[i];
         SkimBToKEE_l2_iso04_dca_ = SkimBToKEE_l2_iso04_dca[i];
         SkimBToKEE_b_iso04_dca_ = SkimBToKEE_b_iso04_dca[i];
         SkimBToKEE_k_iso04_dca_ = SkimBToKEE_k_iso04_dca[i];
         SkimBToKEE_k_svip3d_ = SkimBToKEE_k_svip3d[i];
         SkimBToKEE_k_svip3d_err_ = SkimBToKEE_k_svip3d_err[i];
         SkimBToKEE_l1_n_isotrk_dca_ = SkimBToKEE_l1_n_isotrk_dca[i];
         SkimBToKEE_l2_n_isotrk_dca_ = SkimBToKEE_l2_n_isotrk_dca[i];
         SkimBToKEE_k_n_isotrk_dca_ = SkimBToKEE_k_n_isotrk_dca[i];
         SkimBToKEE_l1_vx_ = SkimBToKEE_l1_vx[i];
         SkimBToKEE_l1_vy_ = SkimBToKEE_l1_vy[i];
         SkimBToKEE_l1_vz_ = SkimBToKEE_l1_vz[i];
         SkimBToKEE_l1_charge_ = SkimBToKEE_l1_charge[i];
         SkimBToKEE_l2_vx_ = SkimBToKEE_l2_vx[i];
         SkimBToKEE_l2_vy_ = SkimBToKEE_l2_vy[i];
         SkimBToKEE_l2_vz_ = SkimBToKEE_l2_vz[i];
         SkimBToKEE_l2_charge_ = SkimBToKEE_l2_charge[i];
         SkimBToKEE_k_vx_ = SkimBToKEE_k_vx[i];
         SkimBToKEE_k_vy_ = SkimBToKEE_k_vy[i];
         SkimBToKEE_k_vz_ = SkimBToKEE_k_vz[i];
         SkimBToKEE_k_dca_sig_ = SkimBToKEE_k_dca_sig[i];
         SkimBToKEE_k_dz_ = SkimBToKEE_k_dz[i];
         SkimBToKEE_kMu_matched_ = SkimBToKEE_kMu_matched[i];
         SkimBToKEE_k_charge_ = SkimBToKEE_k_charge[i];
         SkimBToKEE_l_xy_sig_ = SkimBToKEE_l_xy_sig[i];
         SkimBToKEE_l1l2Dz_ = SkimBToKEE_l1l2Dz[i];
         SkimBToKEE_lKDz_ = SkimBToKEE_lKDz[i];
         SkimBToKEE_l1l2Dr_ = SkimBToKEE_l1l2Dr[i];
         SkimBToKEE_lKDr_ = SkimBToKEE_lKDr[i];
         SkimBToKEE_kl_massKPi_ = SkimBToKEE_kl_massKPi[i];
         SkimBToKEE_p_assymetry_ = SkimBToKEE_p_assymetry[i];
         float SkimBToKEE_k_svip3d_sig_ = SkimBToKEE_k_svip3d_/SkimBToKEE_k_svip3d_err_;
         float SkimBToKEE_fit_pt_over_mass_= SkimBToKEE_fit_pt_/SkimBToKEE_fit_mass_;
         float SkimBToKEE_fit_l1_pt_over_mass_ = SkimBToKEE_fit_l1_pt_/SkimBToKEE_fit_mass_;
         float SkimBToKEE_fit_l2_pt_over_mass_ = SkimBToKEE_fit_l2_pt_/SkimBToKEE_fit_mass_;
         float SkimBToKEE_fit_k_pt_over_mass_ = SkimBToKEE_fit_k_pt_/SkimBToKEE_fit_mass_;
         
         //bdt1 inputvector
         std::vector<float> input{SkimBToKEE_fit_pt_,SkimBToKEE_fit_eta_,SkimBToKEE_fit_phi_,
            SkimBToKEE_D0_mass_LepToK_KToPi_,SkimBToKEE_D0_mass_LepToPi_KToK_,SkimBToKEE_fit_cos2D_,
            SkimBToKEE_svprob_,SkimBToKEE_b_iso04_,SkimBToKEE_fit_k_pt_,SkimBToKEE_fit_k_eta_,
            SkimBToKEE_fit_k_phi_,SkimBToKEE_fit_l1_pt_,SkimBToKEE_fit_l1_eta_,SkimBToKEE_fit_l1_phi_,
            SkimBToKEE_fit_l2_pt_,SkimBToKEE_fit_l2_eta_,SkimBToKEE_fit_l2_phi_,SkimBToKEE_l1_iso04_,
            SkimBToKEE_l2_iso04_,SkimBToKEE_k_iso04_,SkimBToKEE_l1_PFMvaID_Fall17_,SkimBToKEE_l2_PFMvaID_Fall17_,
            SkimBToKEE_l1_PFMvaID_retrained_,SkimBToKEE_l2_PFMvaID_retrained_,SkimBToKEE_l1_iso04_dca_,
            SkimBToKEE_l2_iso04_dca_,SkimBToKEE_b_iso04_dca_,SkimBToKEE_k_iso04_dca_,SkimBToKEE_l1_n_isotrk_dca_,
            SkimBToKEE_l2_n_isotrk_dca_,SkimBToKEE_k_n_isotrk_dca_,SkimBToKEE_l_xy_sig_,SkimBToKEE_l1l2Dz_,
            SkimBToKEE_lKDz_,SkimBToKEE_l1l2Dr_,SkimBToKEE_p_assymetry_,SkimBToKEE_k_svip3d_sig_,
            SkimBToKEE_fit_pt_over_mass_,SkimBToKEE_fit_l1_pt_over_mass_,SkimBToKEE_fit_l2_pt_over_mass_,SkimBToKEE_fit_k_pt_over_mass_};
         BDTSCORE_1 = bdt1(input.data());



         outTree_->Fill();
         SkimBToKEE_fit_pt_ = -1000.;
         SkimBToKEE_fit_eta_ = -1000.;
         SkimBToKEE_fit_phi_ = -1000.;
         SkimBToKEE_D0_mass_LepToK_KToPi_ = -1000.;
         SkimBToKEE_D0_mass_LepToPi_KToK_ = -1000.;
         SkimBToKEE_fit_mass_ = -1000.;
         SkimBToKEE_l_xy_ = -1000.;
         SkimBToKEE_l_xy_unc_ = -1000.;
         SkimBToKEE_fit_cos2D_ = -1000.;
         SkimBToKEE_svprob_ = -1000.;
         SkimBToKEE_fit_massErr_ = -1000.;
         SkimBToKEE_b_iso04_ = -1000.;
         SkimBToKEE_mll_fullfit_ = -1000.;
         SkimBToKEE_vtx_x_ = -1000.;
         SkimBToKEE_vtx_y_ = -1000.;
         SkimBToKEE_vtx_z_ = -1000.;
         SkimBToKEE_l1Idx_ = -1000.;
         SkimBToKEE_l2Idx_ = -1000.;
         SkimBToKEE_kIdx_ = -1000.;
         SkimBToKEE_fit_k_pt_ = -1000.;
         SkimBToKEE_fit_k_eta_ = -1000.;
         SkimBToKEE_fit_k_phi_ = -1000.;
         SkimBToKEE_fit_l1_pt_ = -1000.;
         SkimBToKEE_fit_l1_eta_ = -1000.;
         SkimBToKEE_fit_l1_phi_ = -1000.;
         SkimBToKEE_fit_l2_pt_ = -1000.;
         SkimBToKEE_fit_l2_eta_ = -1000.;
         SkimBToKEE_fit_l2_phi_ = -1000.;
         SkimBToKEE_l1_iso04_ = -1000.;
         SkimBToKEE_l2_iso04_ = -1000.;
         SkimBToKEE_l1_isPF_ = -1000.;
         SkimBToKEE_l2_isPF_ = -1000.;
         SkimBToKEE_k_iso04_ = -1000.;
         SkimBToKEE_l1_isPFoverlap_ = -1000.;
         SkimBToKEE_l2_isPFoverlap_ = -1000.;
         SkimBToKEE_l1_LooseId_ = -1000.;
         SkimBToKEE_l2_LooseId_ = -1000.;
         SkimBToKEE_l1_MediumId_ = -1000.;
         SkimBToKEE_l2_MediumId_ = -1000.;
         SkimBToKEE_l1_TightId_ = -1000.;
         SkimBToKEE_l2_TightId_ = -1000.;
         SkimBToKEE_l1_ConvVeto_ = -1000.;
         SkimBToKEE_l2_ConvVeto_ = -1000.;
         SkimBToKEE_l1_PFMvaID_Fall17_ = -1000.;
         SkimBToKEE_l2_PFMvaID_Fall17_ = -1000.;
         SkimBToKEE_l1_PFMvaID_retrained_ = -1000.;
         SkimBToKEE_l2_PFMvaID_retrained_ = -1000.;
         SkimBToKEE_l1_iso04_dca_ = -1000.;
         SkimBToKEE_l2_iso04_dca_ = -1000.;
         SkimBToKEE_b_iso04_dca_ = -1000.;
         SkimBToKEE_k_iso04_dca_ = -1000.;
         SkimBToKEE_k_svip3d_ = -1000.;
         SkimBToKEE_k_svip3d_err_ = -1000.;
         SkimBToKEE_l1_n_isotrk_dca_ = -1000.;
         SkimBToKEE_l2_n_isotrk_dca_ = -1000.;
         SkimBToKEE_k_n_isotrk_dca_ = -1000.;
         SkimBToKEE_l1_vx_ = -1000.;
         SkimBToKEE_l1_vy_ = -1000.;
         SkimBToKEE_l1_vz_ = -1000.;
         SkimBToKEE_l1_charge_ = -1000.;
         SkimBToKEE_l2_vx_ = -1000.;
         SkimBToKEE_l2_vy_ = -1000.;
         SkimBToKEE_l2_vz_ = -1000.;
         SkimBToKEE_l2_charge_ = -1000.;
         SkimBToKEE_k_vx_ = -1000.;
         SkimBToKEE_k_vy_ = -1000.;
         SkimBToKEE_k_vz_ = -1000.;
         SkimBToKEE_k_dca_sig_ = -1000.;
         SkimBToKEE_k_dz_ = -1000.;
         SkimBToKEE_kMu_matched_ = -1000.;
         SkimBToKEE_k_charge_ = -1000.;
         SkimBToKEE_l_xy_sig_ = -1000.;
         SkimBToKEE_l1l2Dz_ = -1000.;
         SkimBToKEE_lKDz_ = -1000.;
         SkimBToKEE_l1l2Dr_ = -1000.;
         SkimBToKEE_lKDr_ = -1000.;
         SkimBToKEE_kl_massKPi_ = -1000.;
         SkimBToKEE_p_assymetry_ = -1000.;             
         BDTSCORE_1 = -1000.;        
      }
      run_ = -1000;
      luminosityBlock_ = -1000;
      event_ = -1000;
      bunchCrossing_ = -1000;
      nBToKEE_ = -1000;
      nElectron_ = -1000;
      nProbeTracks_ = -1000;
      HLT_DoubleEle10_eta1p22_mMax6_ = false;
      HLT_DoubleEle9p5_eta1p22_mMax6_ = false;
      HLT_DoubleEle9_eta1p22_mMax6_ = false;
      HLT_DoubleEle8p5_eta1p22_mMax6_ = false;
      HLT_DoubleEle8_eta1p22_mMax6_ = false;
      HLT_DoubleEle7p5_eta1p22_mMax6_ = false;
      HLT_DoubleEle7_eta1p22_mMax6_ = false;
      HLT_DoubleEle6p5_eta1p22_mMax6_ = false;
      HLT_DoubleEle6_eta1p22_mMax6_ = false;
      HLT_DoubleEle5p5_eta1p22_mMax6_ = false;
      HLT_DoubleEle5_eta1p22_mMax6_ = false;
      HLT_DoubleEle4p5_eta1p22_mMax6_ = false;
      HLT_DoubleEle4_eta1p22_mMax6_ = false;
      L1_DoubleEG11_er1p2_dR_Max0p6_ = false;
      L1_DoubleEG10p5_er1p2_dR_Max0p6_ = false;
      L1_DoubleEG10_er1p2_dR_Max0p6_ = false;
      L1_DoubleEG9p5_er1p2_dR_Max0p6_ = false;
      L1_DoubleEG9_er1p2_dR_Max0p7_ = false;
      L1_DoubleEG8p5_er1p2_dR_Max0p7_ = false;
      L1_DoubleEG8_er1p2_dR_Max0p7_ = false;
      L1_DoubleEG7p5_er1p2_dR_Max0p7_ = false;
      L1_DoubleEG7_er1p2_dR_Max0p8_ = false;
      L1_DoubleEG6p5_er1p2_dR_Max0p8_ = false;
      L1_DoubleEG6_er1p2_dR_Max0p8_ = false;
      L1_DoubleEG5p5_er1p2_dR_Max0p8_ = false;
      L1_DoubleEG5_er1p2_dR_Max0p9_ = false;
      L1_DoubleEG4p5_er1p2_dR_Max0p9_ = false;
      L1_DoubleEG4_er1p2_dR_Max0p9_ = false;
      nTrigObj_ = -1000;
      nOtherPV_ = -1000;
      PV_ndof_ = -1000;
      PV_x_ = -1000;
      PV_y_ = -1000;
      PV_z_ = -1000;
      PV_chi2_ = -1000;
      PV_score_ = -1000;
      PV_npvs_ = -1000;
      PV_npvsGood_ = -1000;
      nSV_ = -1000;
      nSkimBToKEE_ = -1000;
      
      //std::cout<<HLT_DoubleEle10_eta1p22_mMax6_<<std::endl;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry>5000){
         break;
      }
   }
   outFile_->Write();
   outFile_->Close();
}




int main(int argc, char* argv[]){
    if (argc != 4) {
        std::cout << "Usage: program_name input_file_list output_file output_dir\n";
        return 1;
    }
    std::string inputFileList = argv[1];
    std::string outputFile = argv[2];
    std::string outputDir = argv[3];
    std::ifstream inputFile(inputFileList);
    std::cout<<outputFile<<std::endl;
    std::cout<<outputDir<<std::endl;
    if (!inputFile.is_open()) {
        std::cout << "Failed to open the input file list.\n";
        return 1;
    }

    TChain * chain = new TChain("Events","");
    std::vector<std::string> filenamesvector;
    std::string line;
    while (std::getline(inputFile, line)) {
        filenamesvector.push_back(line);
    }


    for(auto f: filenamesvector){
        std::cout<<f<<"reeeeeee"<<std::endl;
        TString iname(f+"/Events");
        chain->Add(iname);
    }
    //chain->Add("/eos/user/j/jodedra/fastforesttester/CMSSW_13_1_0/src/smallsample_Chunk0.root/Events");
    std::cout<<chain->GetEntries()<<std::endl;
    Events t(chain);
    t.Loop(outputFile,outputDir);
    return 0;
}
