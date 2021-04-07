#include <iostream>
#include "HHKinFit/HHKinFit/interface/HHKinFitMaster.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <math.h> 
#include "../interface/HelperFunctions.h"
#include <boost/program_options.hpp>

using boost::starts_with;
namespace po = boost::program_options;

int main(int argc, char* argv[])
{
  std::string input = "output.root";
  std::string folder = "mt_nominal";
  std::string tree = "ntuple";
  unsigned int first_entry = 0;
  unsigned int last_entry = 9;
  po::variables_map vm;
  po::options_description config("configuration");
  config.add_options()("input",
                       po::value<std::string>(&input)->default_value(input))(
      "folder", po::value<std::string>(&folder)->default_value(folder))(
      "tree", po::value<std::string>(&tree)->default_value(tree))(
      "first_entry",
      po::value<unsigned int>(&first_entry)->default_value(first_entry))(
      "last_entry",
      po::value<unsigned int>(&last_entry)->default_value(last_entry));
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);

  // Access input file and tree
  auto in = TFile::Open(input.c_str(), "read");
  auto dir = (TDirectoryFile *)in->Get(folder.c_str());
  auto inputtree = (TTree *)dir->Get(tree.c_str());

  // Initialize output file
  auto outputname =
      outputname_from_settings(input, folder, first_entry, last_entry);
  boost::filesystem::create_directories(filename_from_inputpath(input));
  auto out = TFile::Open(outputname.c_str(), "recreate");
  out->mkdir(folder.c_str());
  out->cd(folder.c_str());

  std::cout << outputname << std::endl;

  // Create output tree
  auto mkinfitfriend = new TTree("ntuple", "Tau trigger friend tree"); 


  Float_t kinfit_mH,kinfit_mh2,kinfit_chi2,kinfit_prob,kinfit_pull1,kinfit_pull2,kinfit_pullB;
  Int_t kinfit_convergence;

  Float_t kinfit_mH_HsmTune,kinfit_mh2_HsmTune,kinfit_chi2_HsmTune,kinfit_prob_HsmTune,kinfit_pull1_HsmTune,kinfit_pull2_HsmTune,kinfit_pullB_HsmTune;
  Int_t kinfit_convergence_HsmTune;

  
  //Float_t kinfit_res1, kinfit_res2; 

  Float_t bpt1_bRegRes, bpt2_bRegRes;
  //Float_t bReg1;
  


  mkinfitfriend->Branch("kinfit_mH",&kinfit_mH,"kinfit_mH/F");
  mkinfitfriend->Branch("kinfit_mh2",&kinfit_mh2,"kinfit_mh2/F");
  mkinfitfriend->Branch("kinfit_chi2",&kinfit_chi2,"kinfit_chi2/F");
  mkinfitfriend->Branch("kinfit_prob",&kinfit_prob,"kinfit_prob/F");
  mkinfitfriend->Branch("kinfit_pull1",&kinfit_pull1,"kinfit_pull1/F");
  mkinfitfriend->Branch("kinfit_pull2",&kinfit_pull2,"kinfit_pull2/F");
  mkinfitfriend->Branch("kinfit_pullB",&kinfit_pullB,"kinfit_pullB/F");
  mkinfitfriend->Branch("kinfit_convergence",&kinfit_convergence,"kinfit_convergence/I");

  mkinfitfriend->Branch("kinfit_mH_HsmTune",&kinfit_mH_HsmTune,"kinfit_mH_HsmTune/F");  
  mkinfitfriend->Branch("kinfit_mh2_HsmTune",&kinfit_mh2_HsmTune,"kinfit_mh2_HsmTune/F");
  mkinfitfriend->Branch("kinfit_chi2_HsmTune",&kinfit_chi2_HsmTune,"kinfit_chi2_HsmTune/F");
  mkinfitfriend->Branch("kinfit_prob_HsmTune",&kinfit_prob_HsmTune,"kinfit_prob_HsmTune/F");
  mkinfitfriend->Branch("kinfit_pull1_HsmTune",&kinfit_pull1_HsmTune,"kinfit_pull1_HsmTune/F");
  mkinfitfriend->Branch("kinfit_pull2_HsmTune",&kinfit_pull2_HsmTune,"kinfit_pull2_HsmTune/F");
  mkinfitfriend->Branch("kinfit_pullB_HsmTune",&kinfit_pullB_HsmTune,"kinfit_pullB_HsmTune/F");
  mkinfitfriend->Branch("kinfit_convergence_HsmTune",&kinfit_convergence_HsmTune,"kinfit_convergence_HsmTune/I");

  /*

  mkinfitfriend->Branch("kinfit_res1", &kinfit_res1, "kinfit_res1/F");
  mkinfitfriend->Branch("kinfit_res2", &kinfit_res2, "kinfit_res2/F");
  mkinfitfriend->Branch("bpt1_bRegRes", &bpt1_bRegRes, "bpt1_bRegRes/F");
  mkinfitfriend->Branch("bpt2_bRegRes", &bpt2_bRegRes, "bpt2_bRegRes/F");

  mkinfitfriend->Branch("bpt_bReg_1",&bReg1,"bReg1/F");
  */

   //Leaf types  
   Int_t           njets;
   Int_t           nbtag;
  //  Double_t        b1_dR;
  //  Float_t        bcsv_1;
   Float_t        bpt_1;
   Float_t        beta_1;
   Float_t        bphi_1;
   Float_t        bm_1;
   Float_t        jpt;
   Float_t        jeta;
   Float_t        jphi;
   Float_t        jm;
  //  Float_t        b2_dR;
  //  Float_t        bcsv_2;
   Float_t        bpt_2;
   Float_t        beta_2;
   Float_t        bphi_2;
   Float_t        bm_2;
   Float_t       pt_1;
   Float_t       eta_1;
   Float_t       phi_1;
   Float_t       m_1;
   Float_t       pt_2;
   Float_t       eta_2;
   Float_t       phi_2;
   Float_t       m_2;
   Float_t        puppimet;
   Float_t        puppimetphi;

   Float_t        puppimetcov00;
   Float_t        puppimetcov10;
   Float_t        puppimetcov01;
   Float_t        puppimetcov11;
   Float_t        byTightDeepTau2017v2p1VSjet_2;

  Float_t        bJet_bRegRes;
  Float_t        bJet2_bRegRes;

  Float_t       HsmTune_DiJetSystemJet1Pt_bReg;
  Float_t       HsmTune_DiJetSystemJet1Eta_bReg;
  Float_t       HsmTune_DiJetSystemJet1Phi_bReg;
  Float_t       HsmTune_DiJetSystemJet1mass_bReg;
  Float_t       HsmTune_DiJetSystemJet2Pt_bReg;
  Float_t       HsmTune_DiJetSystemJet2Eta_bReg;
  Float_t       HsmTune_DiJetSystemJet2Phi_bReg;
  Float_t       HsmTune_DiJetSystemJet2mass_bReg;
  Float_t       HsmTune_JetPlusBJetSystemBJetPt_bReg;
  Float_t       HsmTune_JetPlusBJetSystemBJetEta_bReg;
  Float_t       HsmTune_JetPlusBJetSystemBJetPhi_bReg;
  Float_t       HsmTune_JetPlusBJetSystemBJetmass_bReg;
  Float_t       HsmTune_JetPlusBJetSystemJetPt_bReg;
  Float_t       HsmTune_JetPlusBJetSystemJetEta_bReg;
  Float_t       HsmTune_JetPlusBJetSystemJetPhi_bReg;
  Float_t       HsmTune_JetPlusBJetSystemJetmass_bReg;

    // List of branches
   TBranch        *b_njets;
   TBranch        *b_nbtag;
  //  TBranch        *b_b1_dR;
  //  TBranch        *b_bcsv_1;
   TBranch        *b_bpt_1;
   TBranch        *b_beta_1;
   TBranch        *b_bphi_1;
   TBranch        *b_bm_1;
   TBranch        *b_jpt;
   TBranch        *b_jeta;
   TBranch        *b_jphi;
   TBranch        *b_jm;
  //  TBranch        *b_b2_dR;
  //  TBranch        *b_bcsv_2;
   TBranch        *b_bpt_2;
   TBranch        *b_beta_2;
   TBranch        *b_bphi_2;
   TBranch        *b_bm_2;
  //  TBranch        *b_b1_dR;
   TBranch        *b_pt_1;
   TBranch        *b_eta_1;
   TBranch        *b_phi_1;
   TBranch        *b_m_1;
  //  TBranch        *b_2_dR;
   TBranch        *b_pt_2;
   TBranch        *b_eta_2;
   TBranch        *b_phi_2;
   TBranch        *b_m_2;
   TBranch        *b_puppimet;
   TBranch        *b_puppimetphi;

   TBranch        *b_puppimetcov00;
   TBranch        *b_puppimetcov10;
   TBranch        *b_puppimetcov01;
   TBranch        *b_puppimetcov11;
   TBranch        *b_byTightDeepTau2017v2p1VSjet_2;

  TBranch         *b_bJet_bRegRes;
  TBranch         *b_bJet2_bRegRes;


  TBranch       *b_HsmTune_DiJetSystemJet1Pt_bReg;
  TBranch       *b_HsmTune_DiJetSystemJet1Eta_bReg;
  TBranch       *b_HsmTune_DiJetSystemJet1Phi_bReg;
  TBranch       *b_HsmTune_DiJetSystemJet1mass_bReg;
  TBranch       *b_HsmTune_DiJetSystemJet2Pt_bReg;
  TBranch       *b_HsmTune_DiJetSystemJet2Eta_bReg;
  TBranch       *b_HsmTune_DiJetSystemJet2Phi_bReg;
  TBranch       *b_HsmTune_DiJetSystemJet2mass_bReg;
  TBranch       *b_HsmTune_JetPlusBJetSystemBJetPt_bReg;
  TBranch       *b_HsmTune_JetPlusBJetSystemBJetEta_bReg;
  TBranch       *b_HsmTune_JetPlusBJetSystemBJetPhi_bReg;
  TBranch       *b_HsmTune_JetPlusBJetSystemBJetmass_bReg;
  TBranch       *b_HsmTune_JetPlusBJetSystemJetPt_bReg;
  TBranch       *b_HsmTune_JetPlusBJetSystemJetEta_bReg;
  TBranch       *b_HsmTune_JetPlusBJetSystemJetPhi_bReg;
  TBranch       *b_HsmTune_JetPlusBJetSystemJetmass_bReg;

   inputtree->SetBranchAddress("njets", &njets, &b_njets);
   inputtree->SetBranchAddress("nbtag", &nbtag, &b_nbtag);
   
  //  inputtree->SetBranchAddress("bcsv_1", &bcsv_1, &b_bcsv_1);
   inputtree->SetBranchAddress("bpt_bReg_1", &bpt_1, &b_bpt_1);
   inputtree->SetBranchAddress("bphi_bReg_1", &bphi_1, &b_bphi_1);
   inputtree->SetBranchAddress("beta_bReg_1", &beta_1, &b_beta_1);
   inputtree->SetBranchAddress("bm_bReg_1", &bm_1, &b_bm_1);
  //  inputtree->SetBranchAddress("bcsv_2", &bcsv_2, &b_bcsv_2);
   inputtree->SetBranchAddress("bpt_bReg_2", &bpt_2, &b_bpt_2);
   inputtree->SetBranchAddress("bphi_bReg_2", &bphi_2, &b_bphi_2);
   inputtree->SetBranchAddress("beta_bReg_2", &beta_2, &b_beta_2);
   inputtree->SetBranchAddress("bm_bReg_2", &bm_2, &b_bm_2);
   inputtree->SetBranchAddress("highCSVjetUsedFordiBJetSystemPt_bReg", &jpt, &b_jpt);
   inputtree->SetBranchAddress("highCSVjetUsedFordiBJetSystemEta_bReg", &jeta, &b_jeta);
   inputtree->SetBranchAddress("highCSVjetUsedFordiBJetSystemPhi_bReg", &jphi, &b_jphi);
   inputtree->SetBranchAddress("highCSVjetUsedFordiBJetSystemMass_bReg", &jm, &b_jm);
   inputtree->SetBranchAddress("pt_1", &pt_1, &b_pt_1);
   inputtree->SetBranchAddress("phi_1", &phi_1, &b_phi_1);
   inputtree->SetBranchAddress("eta_1", &eta_1, &b_eta_1);
   inputtree->SetBranchAddress("m_1", &m_1, &b_m_1);
   inputtree->SetBranchAddress("pt_2", &pt_2, &b_pt_2);
   inputtree->SetBranchAddress("phi_2", &phi_2, &b_phi_2);
   inputtree->SetBranchAddress("eta_2", &eta_2, &b_eta_2);
   inputtree->SetBranchAddress("m_2", &m_2, &b_m_2);
   inputtree->SetBranchAddress("puppimet", &puppimet, &b_puppimet);
   inputtree->SetBranchAddress("puppimetphi", &puppimetphi, &b_puppimetphi);
   inputtree->SetBranchAddress("puppimetcov00", &puppimetcov00, &b_puppimetcov00);
   inputtree->SetBranchAddress("puppimetcov10", &puppimetcov10, &b_puppimetcov10);
   inputtree->SetBranchAddress("puppimetcov01", &puppimetcov01, &b_puppimetcov01);
   inputtree->SetBranchAddress("puppimetcov11", &puppimetcov11, &b_puppimetcov11);
   inputtree->SetBranchAddress("byTightDeepTau2017v2p1VSjet_2", &byTightDeepTau2017v2p1VSjet_2, &b_byTightDeepTau2017v2p1VSjet_2);
   inputtree->SetBranchAddress("bJet_bRegRes", &bJet_bRegRes, &b_bJet_bRegRes);
   inputtree->SetBranchAddress("bJet2_bRegRes", &bJet2_bRegRes, &b_bJet2_bRegRes);

   inputtree->SetBranchAddress("HsmTune_DiJetSystemJet1Pt_bReg", &HsmTune_DiJetSystemJet1Pt_bReg, &b_HsmTune_DiJetSystemJet1Pt_bReg);
   inputtree->SetBranchAddress("HsmTune_DiJetSystemJet1Eta_bReg", &HsmTune_DiJetSystemJet1Eta_bReg, &b_HsmTune_DiJetSystemJet1Eta_bReg);
   inputtree->SetBranchAddress("HsmTune_DiJetSystemJet1Phi_bReg", &HsmTune_DiJetSystemJet1Phi_bReg, &b_HsmTune_DiJetSystemJet1Phi_bReg);
   inputtree->SetBranchAddress("HsmTune_DiJetSystemJet1mass_bReg", &HsmTune_DiJetSystemJet1mass_bReg, &b_HsmTune_DiJetSystemJet1mass_bReg);
   inputtree->SetBranchAddress("HsmTune_DiJetSystemJet2Pt_bReg", &HsmTune_DiJetSystemJet2Pt_bReg, &b_HsmTune_DiJetSystemJet2Pt_bReg);
   inputtree->SetBranchAddress("HsmTune_DiJetSystemJet2Eta_bReg", &HsmTune_DiJetSystemJet2Eta_bReg, &b_HsmTune_DiJetSystemJet2Eta_bReg);
   inputtree->SetBranchAddress("HsmTune_DiJetSystemJet2Phi_bReg", &HsmTune_DiJetSystemJet2Phi_bReg, &b_HsmTune_DiJetSystemJet2Phi_bReg);
   inputtree->SetBranchAddress("HsmTune_DiJetSystemJet2mass_bReg", &HsmTune_DiJetSystemJet2mass_bReg, &b_HsmTune_DiJetSystemJet2mass_bReg);

   inputtree->SetBranchAddress("HsmTune_JetPlusBJetSystemBJetPt_bReg", &HsmTune_JetPlusBJetSystemBJetPt_bReg, &b_HsmTune_JetPlusBJetSystemBJetPt_bReg);
   inputtree->SetBranchAddress("HsmTune_JetPlusBJetSystemBJetEta_bReg", &HsmTune_JetPlusBJetSystemBJetEta_bReg, &b_HsmTune_JetPlusBJetSystemBJetEta_bReg);
   inputtree->SetBranchAddress("HsmTune_JetPlusBJetSystemBJetPhi_bReg", &HsmTune_JetPlusBJetSystemBJetPhi_bReg, &b_HsmTune_JetPlusBJetSystemBJetPhi_bReg);
   inputtree->SetBranchAddress("HsmTune_JetPlusBJetSystemBJetmass_bReg", &HsmTune_JetPlusBJetSystemBJetmass_bReg, &b_HsmTune_JetPlusBJetSystemBJetmass_bReg);
   inputtree->SetBranchAddress("HsmTune_JetPlusBJetSystemJetPt_bReg", &HsmTune_JetPlusBJetSystemJetPt_bReg, &b_HsmTune_JetPlusBJetSystemJetPt_bReg);
   inputtree->SetBranchAddress("HsmTune_JetPlusBJetSystemJetEta_bReg", &HsmTune_JetPlusBJetSystemJetEta_bReg, &b_HsmTune_JetPlusBJetSystemJetEta_bReg);
   inputtree->SetBranchAddress("HsmTune_JetPlusBJetSystemJetPhi_bReg", &HsmTune_JetPlusBJetSystemJetPhi_bReg, &b_HsmTune_JetPlusBJetSystemJetPhi_bReg);
   inputtree->SetBranchAddress("HsmTune_JetPlusBJetSystemJetmass_bReg", &HsmTune_JetPlusBJetSystemJetmass_bReg, &b_HsmTune_JetPlusBJetSystemJetmass_bReg);

   //define the testd hypotheses
   std::vector<Int_t> hypo_mh1;
   hypo_mh1.push_back(125);

  
   std::vector<Int_t> hypo_mh2 = {5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,180,210,240,270,300,330,360,390,420,450,480,510,540,570,600,630,660,690,720,750,780,810,840,870,900,950,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2500,3000};
   
   std::cout << "first: " << first_entry << " last: " << last_entry << std::endl;

   // event loop
  for (unsigned int i = first_entry; i < last_entry; i++) {
    inputtree->GetEntry(i);
    std::cout << "Event " << i-first_entry << " (" << round(100*float(i-first_entry)/float(last_entry-first_entry)) << "%)" << std::endl;
    if (nbtag < 1) {
      kinfit_mH = -10;
      kinfit_mh2 = -10;
      kinfit_chi2 = -10;
      kinfit_prob = -10;
      kinfit_pull1 = -10;
      kinfit_pull2 = -10;
      kinfit_pullB = -10;
      kinfit_convergence = -10;

      kinfit_mH_HsmTune = -10;
      kinfit_mh2_HsmTune = -10;
      kinfit_chi2_HsmTune = -10;
      kinfit_prob_HsmTune = -10;
      kinfit_pull1_HsmTune = -10;
      kinfit_pull2_HsmTune = -10;
      kinfit_pullB_HsmTune = -10;
      kinfit_convergence_HsmTune = -10;

      mkinfitfriend->Fill();      
      continue;
    }   
    if ((nbtag==1) && (jpt<0)) { // happens when no non-bjet exists or non-bjet is too close in dR to b-jet
      kinfit_mH = -10;
      kinfit_mh2 = -10;
      kinfit_chi2 = -10;
      kinfit_prob = -10;
      kinfit_pull1 = -10;
      kinfit_pull2 = -10;
      kinfit_pullB = -10;
      kinfit_convergence = -10;

      kinfit_mH_HsmTune = -10;
      kinfit_mh2_HsmTune = -10;
      kinfit_chi2_HsmTune = -10;
      kinfit_prob_HsmTune = -10;
      kinfit_pull1_HsmTune = -10;
      kinfit_pull2_HsmTune = -10;
      kinfit_pullB_HsmTune = -10;
      kinfit_convergence_HsmTune = -10;


      /*
      kinfit_res1 = -10;
      bpt1_bRegRes = -10;   
      kinfit_res2 = -10; 
      bpt2_bRegRes = -10;
      */
      mkinfitfriend->Fill();      
      continue;
    }
    
    //define input vectors
    TLorentzVector b1      = TLorentzVector(); b1.SetPtEtaPhiM(bpt_1,beta_1,bphi_1,bm_1);
    TLorentzVector b2      = TLorentzVector(); 
    TLorentzVector b1_HsmTune = TLorentzVector(); 
    TLorentzVector b2_HsmTune = TLorentzVector(); 
    if (nbtag==1)
    { 
      b2.SetPtEtaPhiM(jpt,jeta,jphi,jm);
      b1_HsmTune.SetPtEtaPhiM(HsmTune_JetPlusBJetSystemBJetPt_bReg, HsmTune_JetPlusBJetSystemBJetEta_bReg, HsmTune_JetPlusBJetSystemBJetPhi_bReg, HsmTune_JetPlusBJetSystemBJetmass_bReg);
      b2_HsmTune.SetPtEtaPhiM(HsmTune_JetPlusBJetSystemJetPt_bReg, HsmTune_JetPlusBJetSystemJetEta_bReg, HsmTune_JetPlusBJetSystemJetPhi_bReg, HsmTune_JetPlusBJetSystemJetmass_bReg);
    }
    else {
      b2.SetPtEtaPhiM(bpt_2,beta_2,bphi_2,bm_2);
      b1_HsmTune.SetPtEtaPhiM(HsmTune_DiJetSystemJet1Pt_bReg,HsmTune_DiJetSystemJet1Eta_bReg,HsmTune_DiJetSystemJet1Phi_bReg, HsmTune_DiJetSystemJet1mass_bReg);
      b2_HsmTune.SetPtEtaPhiM(HsmTune_DiJetSystemJet2Pt_bReg,HsmTune_DiJetSystemJet2Eta_bReg,HsmTune_DiJetSystemJet2Phi_bReg, HsmTune_DiJetSystemJet2mass_bReg);
    }
    TLorentzVector tau1vis = TLorentzVector(); tau1vis.SetPtEtaPhiM(pt_1,eta_1,phi_1,m_1);
    TLorentzVector tau2vis = TLorentzVector(); tau2vis.SetPtEtaPhiM(pt_2,eta_2,phi_2,m_2);

    TLorentzVector ptmiss  = TLorentzVector(); ptmiss.SetPtEtaPhiE(puppimet,0,puppimetphi,puppimet);
    TMatrixD puppimetcov(2,2);
    puppimetcov(0,0)=puppimetcov00;
    puppimetcov(1,0)=puppimetcov10;
    puppimetcov(0,1)=puppimetcov01;    
    puppimetcov(1,1)=puppimetcov11;
    

    bool ifbReg = false;

    
    bpt1_bRegRes = bpt_1*bJet_bRegRes;
    bpt2_bRegRes = bpt_2*bJet2_bRegRes;

  
    
    
    //intance of fitter master class

    ifbReg = false; // IMPORTANT!! 
    
    HHKinFitMaster kinFits = HHKinFitMaster(&b1,&b2,&tau1vis,&tau2vis,ifbReg,bpt1_bRegRes,bpt2_bRegRes);
    kinFits.setAdvancedBalance(&ptmiss,puppimetcov);
    kinFits.addMh1Hypothesis(hypo_mh2);
    kinFits.addMh2Hypothesis(hypo_mh1);
    kinFits.doFullFit(ifbReg, bpt1_bRegRes,bpt2_bRegRes);       
    std::pair<Int_t, Int_t> bestHypo = kinFits.getBestHypoFullFit();

    

    HHKinFitMaster kinFits_HsmTune = HHKinFitMaster(&b1_HsmTune,&b2_HsmTune,&tau1vis,&tau2vis,ifbReg,bpt1_bRegRes,bpt2_bRegRes);
    kinFits_HsmTune.setAdvancedBalance(&ptmiss,puppimetcov);
    kinFits_HsmTune.addMh1Hypothesis(hypo_mh2);
    kinFits_HsmTune.addMh2Hypothesis(hypo_mh1);
    kinFits_HsmTune.doFullFit(ifbReg, bpt1_bRegRes,bpt2_bRegRes);       
    std::pair<Int_t, Int_t> bestHypo_HsmTune = kinFits_HsmTune.getBestHypoFullFit();

    

    /* 

    kinfit_res1 = kinFits.GetBjetResolution(b1.Eta(), b1.Et(), false, 0);
    kinfit_res2 = kinFits.GetBjetResolution(b2.Eta(), b2.Et(), false, 0);

    
    if (bpt_1 == -10 && bJet_bRegRes == -10){
      bpt1_bRegRes = -10;
      kinfit_res1  = -10;
    }

    if (bpt_2 == -10 && bJet2_bRegRes == -10){
      bpt2_bRegRes = -10;
      kinfit_res2 = -10; 
    }
    */
  
    
    if(bestHypo.second>0) {
        std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_chi2 = kinFits.getChi2FullFit();
        std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_fitprob = kinFits.getFitProbFullFit();
        std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_mH = kinFits.getMHFullFit();
        std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_pull_b1 = kinFits.getPullB1FullFit();
        std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_pull_b2 = kinFits.getPullB2FullFit();
        std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_pull_balance = kinFits.getPullBalanceFullFit();
        std::map< std::pair<Int_t, Int_t>, Int_t> fit_convergence = kinFits.getConvergenceFullFit();

        kinfit_convergence = fit_convergence.at(bestHypo);
        kinfit_mH = fit_results_mH.at(bestHypo);
        kinfit_mh2 = bestHypo.first;
        kinfit_chi2 = fit_results_chi2.at(bestHypo);
        kinfit_prob = fit_results_fitprob.at(bestHypo);
        kinfit_pull1 = fit_results_pull_b1.at(bestHypo);
        kinfit_pull2 = fit_results_pull_b2.at(bestHypo);
        kinfit_pullB = fit_results_pull_balance.at(bestHypo);
    }
    else {
      kinfit_mH = -10;
      kinfit_mh2 = -10;
      kinfit_chi2 = 999.;
      kinfit_prob = 0.0;
      kinfit_pull1 = -10;
      kinfit_pull2 = -10;
      kinfit_pullB = -10;
      kinfit_convergence = -1;
    }
    
    if(bestHypo_HsmTune.second>0) {
        std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_chi2_HsmTune = kinFits_HsmTune.getChi2FullFit();
        std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_fitprob_HsmTune = kinFits_HsmTune.getFitProbFullFit();
        std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_mH_HsmTune = kinFits_HsmTune.getMHFullFit();
        std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_pull_b1_HsmTune = kinFits_HsmTune.getPullB1FullFit();
        std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_pull_b2_HsmTune = kinFits_HsmTune.getPullB2FullFit();
        std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_pull_balance_HsmTune = kinFits_HsmTune.getPullBalanceFullFit();
        std::map< std::pair<Int_t, Int_t>, Int_t> fit_convergence_HsmTune = kinFits_HsmTune.getConvergenceFullFit();

        kinfit_convergence_HsmTune = fit_convergence_HsmTune.at(bestHypo_HsmTune);
        kinfit_mH_HsmTune = fit_results_mH_HsmTune.at(bestHypo_HsmTune);
        kinfit_mh2_HsmTune = bestHypo_HsmTune.first;
        kinfit_chi2_HsmTune = fit_results_chi2_HsmTune.at(bestHypo_HsmTune);
        kinfit_prob_HsmTune = fit_results_fitprob_HsmTune.at(bestHypo_HsmTune);
        kinfit_pull1_HsmTune = fit_results_pull_b1_HsmTune.at(bestHypo_HsmTune);
        kinfit_pull2_HsmTune = fit_results_pull_b2_HsmTune.at(bestHypo_HsmTune);
        kinfit_pullB_HsmTune = fit_results_pull_balance_HsmTune.at(bestHypo_HsmTune);
    }
    else {
      kinfit_mH_HsmTune = -10;
      kinfit_mh2_HsmTune = -10;
      kinfit_chi2_HsmTune = 999.;
      kinfit_prob_HsmTune = 0.0;
      kinfit_pull1_HsmTune = -10;
      kinfit_pull2_HsmTune = -10;
      kinfit_pullB_HsmTune = -10;
      kinfit_convergence_HsmTune = -1;
    }    
    
    mkinfitfriend->Fill();
  }
  mkinfitfriend->Write();
  out->Close();
  in->Close();
  return (0);
}


// checklist for ttbb <-> bbtt : second <-> first, mh1<->mh2
