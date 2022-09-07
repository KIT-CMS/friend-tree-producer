#include "TauAnalysis/ClassicSVfit/interface/FastMTT.h"
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

#include "TH1F.h"
#include "TVector2.h"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

#include <stdexcept>

#include "HiggsAnalysis/friend-tree-producer/interface/HelperFunctions.h"

using namespace classic_svFit;
using boost::starts_with;
namespace po = boost::program_options;

std::pair<MeasuredTauLepton::kDecayType, MeasuredTauLepton::kDecayType>
folder_to_ditaudecay(std::string channel) {
  if (channel == "em")
    return std::make_pair(MeasuredTauLepton::kTauToElecDecay,
                          MeasuredTauLepton::kTauToMuDecay);
  else if (channel == "ee")
    return std::make_pair(MeasuredTauLepton::kTauToElecDecay,
                          MeasuredTauLepton::kTauToElecDecay);
  else if (channel == "mm")
    return std::make_pair(MeasuredTauLepton::kTauToMuDecay,
                          MeasuredTauLepton::kTauToMuDecay);
  else if (channel == "et")
    return std::make_pair(MeasuredTauLepton::kTauToElecDecay,
                          MeasuredTauLepton::kTauToHadDecay);
  else if (channel == "mt")
    return std::make_pair(MeasuredTauLepton::kTauToMuDecay,
                          MeasuredTauLepton::kTauToHadDecay);
  else if (channel == "tt")
    return std::make_pair(MeasuredTauLepton::kTauToHadDecay,
                          MeasuredTauLepton::kTauToHadDecay);
  else {
    std::cout << "Channel determined in a wrong way. Exiting";
    exit(1);
  }
  return std::make_pair(MeasuredTauLepton::kUndefinedDecayType,
                        MeasuredTauLepton::kUndefinedDecayType);
}

int main(int argc, char **argv) {
  std::string input = "output.root";
  std::string output_dir = "";
  std::string folder = "mt_nominal";
  std::string tree = "ntuple";
  std::string channel = "mt";
  std::string era = "2017";
  int first_entry = 0;
  int last_entry = -1;
  bool organize_outputs = true;
  po::variables_map vm;
  po::options_description config("configuration");
  config.add_options()("input",
                       po::value<std::string>(&input)->default_value(input))(
      "output_dir",
      po::value<std::string>(&output_dir)->default_value(output_dir))(
      "folder", po::value<std::string>(&folder)->default_value(folder))(
      "tree", po::value<std::string>(&tree)->default_value(tree))(
      "first_entry", po::value<int>(&first_entry)->default_value(first_entry))(
      "last_entry", po::value<int>(&last_entry)->default_value(last_entry))(
      "channel", po::value<std::string>(&channel)->default_value(channel))(
      "era", po::value<std::string>(&era)->default_value(era))(
      "organize_outputs",
      po::value<bool>(&organize_outputs)->default_value(organize_outputs));
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);
  std::cout << "Setting up FastMTT for channel " << channel << " and era "
            << era << " pipeline " << folder << std::endl;
  // Access input file and tree
  TFile *in = TFile::Open(input.c_str(), "read");
  // TDirectoryFile *dir = (TDirectoryFile *)in->Get(folder.c_str());
  TTree *inputtree = (TTree *)in->Get(tree.c_str());

  // Restrict input tree to needed branches
  inputtree->SetBranchStatus("*", 0);

  // Setting events processing ranges
  std::vector<std::string> required_quantities{
      "pt_1",     "eta_1",    "phi_1",    "mass_1",      "pt_2",
      "eta_2",    "phi_2",    "mass_2",   "decaymode_2", "pfmet",
      "metcov00", "metcov01", "metcov10", "metcov11",    "pfmetphi"};
  int include_last_ev = 1;
  bool run_required = false;
  if (last_entry < 0 || last_entry >= inputtree->GetEntries()) {
    last_entry = inputtree->GetEntries();
    include_last_ev = 0;
  }

  std::vector<std::string> quantities =
      find_quantities(inputtree, required_quantities);
  std::map<std::string, std::string> quantity_names = build_quantities_map(
      folder, quantities, required_quantities, run_required, false);

  // Initialize output file
  std::cout << "Initializing output file" << std::endl;
  std::cout << "Passing args: " << input << " " << folder << " " << era << " "
            << channel << " " << first_entry << " " << last_entry << " "
            << organize_outputs << std::endl;
  std::string outputname = outputname_from_settings_crown(
      input, folder, first_entry, last_entry, output_dir, era, channel,
      organize_outputs);
  if (organize_outputs) {
    std::cout << "Output file: " << outputname << std::endl;
    std::cout << "Creating folder: " << filename_from_inputpath(input)
              << std::endl;
    boost::filesystem::create_directories(era + "_" + channel + "_" +
                                          filename_from_inputpath(input));
  }
  TFile *out = TFile::Open(outputname.c_str(), "recreate");
  std::cout << "Output file: " << outputname << std::endl;
  std::cout << "input file: " << input << std::endl;
  std::cout << "filename_from_inputpath " << filename_from_inputpath(input)
            << std::endl;
  // out->mkdir(folder.c_str());
  // out->cd(folder.c_str());

  // Create output tree
  TTree *svfitfriend = new TTree("ntuple", "svfit friend tree");

  // do not rerun everything for jes and jer
  if (folder.find("jes") != std::string::npos ||
      folder.find("jer") != std::string::npos) {
    run_required = false;
  }

  if (run_required == false) {
    std::cout << "No calculation needed, will write an empty tree" << std::endl;
    svfitfriend->Write("", TObject::kOverwrite);
    out->Close();
    in->Close();
    std::cout << "done" << std::endl;
    return 0;
  } else {
    // Now we can set the branches to be used
    // Quantities of first lepton
    std::cout << "Setting branches for first lepton" << std::endl;

    inputtree->SetBranchStatus(quantity_names["pt_1"].c_str(), 1);
    inputtree->SetBranchStatus(quantity_names["eta_1"].c_str(), 1);
    inputtree->SetBranchStatus(quantity_names["phi_1"].c_str(), 1);
    inputtree->SetBranchStatus(quantity_names["mass_1"].c_str(), 1);
    Float_t pt_1, eta_1, phi_1, m_1;
    Int_t decayMode_1;
    inputtree->SetBranchAddress(quantity_names["pt_1"].c_str(), &pt_1);
    inputtree->SetBranchAddress(quantity_names["eta_1"].c_str(), &eta_1);
    inputtree->SetBranchAddress(quantity_names["phi_1"].c_str(), &phi_1);
    inputtree->SetBranchAddress(quantity_names["mass_1"].c_str(), &m_1);

    // Quantities of second lepton
    std::cout << "Setting branches for second lepton" << std::endl;
    inputtree->SetBranchStatus(quantity_names["pt_2"].c_str(), 1);
    inputtree->SetBranchStatus(quantity_names["eta_2"].c_str(), 1);
    inputtree->SetBranchStatus(quantity_names["phi_2"].c_str(), 1);
    inputtree->SetBranchStatus(quantity_names["mass_2"].c_str(), 1);
    inputtree->SetBranchStatus(quantity_names["decaymode_2"].c_str(), 1);
    Float_t pt_2, eta_2, phi_2, m_2;
    Int_t decayMode_2;
    inputtree->SetBranchAddress(quantity_names["pt_2"].c_str(), &pt_2);
    inputtree->SetBranchAddress(quantity_names["eta_2"].c_str(), &eta_2);
    inputtree->SetBranchAddress(quantity_names["phi_2"].c_str(), &phi_2);
    inputtree->SetBranchAddress(quantity_names["mass_2"].c_str(), &m_2);
    inputtree->SetBranchAddress(quantity_names["decaymode_2"].c_str(),
                                &decayMode_2);

    // Quantities of MET
    std::cout << "Setting branches for MET" << std::endl;
    inputtree->SetBranchStatus(quantity_names["pfmet"].c_str(), 1);
    inputtree->SetBranchStatus(quantity_names["metcov00"].c_str(), 1);
    inputtree->SetBranchStatus(quantity_names["metcov01"].c_str(), 1);
    inputtree->SetBranchStatus(quantity_names["metcov10"].c_str(), 1);
    inputtree->SetBranchStatus(quantity_names["metcov11"].c_str(), 1);
    inputtree->SetBranchStatus(quantity_names["pfmetphi"].c_str(), 1);
    Float_t met, metcov00, metcov01, metcov10, metcov11, metphi;
    inputtree->SetBranchAddress(quantity_names["pfmet"].c_str(), &met);
    inputtree->SetBranchAddress(quantity_names["metcov00"].c_str(), &metcov00);
    inputtree->SetBranchAddress(quantity_names["metcov01"].c_str(), &metcov01);
    inputtree->SetBranchAddress(quantity_names["metcov10"].c_str(), &metcov10);
    inputtree->SetBranchAddress(quantity_names["metcov11"].c_str(), &metcov11);
    inputtree->SetBranchAddress(quantity_names["pfmetphi"].c_str(), &metphi);

    if (channel == "tt") {
      inputtree->SetBranchStatus("decaymode_1", 1);
      inputtree->SetBranchAddress("decaymode_1", &decayMode_1);
    }

    // // Quantities of puppi MET
    // inputtree->SetBranchStatus(quantity_names["puppimet"].c_str(),1);
    // inputtree->SetBranchStatus(quantity_names["puppimetcov00"].c_str(),1);
    // inputtree->SetBranchStatus(quantity_names["puppimetcov01"].c_str(),1);
    // inputtree->SetBranchStatus(quantity_names["puppimetcov10"].c_str(),1);
    // inputtree->SetBranchStatus(quantity_names["puppimetcov11"].c_str(),1);
    // inputtree->SetBranchStatus(quantity_names["puppimetphi"].c_str(),1);
    // Float_t
    // puppimet,puppimetcov00,puppimetcov01,puppimetcov10,puppimetcov11,puppimetphi;
    // inputtree->SetBranchAddress(quantity_names["puppimet"],.c_str()&puppimet);
    // inputtree->SetBranchAddress(quantity_names["puppimetcov00"],.c_str()&puppimetcov00);
    // inputtree->SetBranchAddress(quantity_names["puppimetcov01"],.c_str()&puppimetcov01);
    // inputtree->SetBranchAddress(quantity_names["puppimetcov10"],.c_str()&puppimetcov10);
    // inputtree->SetBranchAddress(quantity_names["puppimetcov11"],.c_str()&puppimetcov11);
    // inputtree->SetBranchAddress(quantity_names["puppimetphi"],.c_str()&puppimetphi);

    // FastMTT outputs
    Float_t pt_fastmtt, eta_fastmtt, phi_fastmtt, m_fastmtt;
    std::string branch_name_pt, branch_name_eta, branch_name_phi, branch_name_m,
        branch_name_pt_datatype, branch_name_eta_datatype,
        branch_name_phi_datatype, branch_name_m_datatype;

    if (folder == "nominal") {
      std::string branch_name_pt = "pt_fastmtt";
      std::string branch_name_eta = "eta_fastmtt";
      std::string branch_name_phi = "phi_fastmtt";
      std::string branch_name_m = "m_fastmtt";
      std::string branch_name_pt_datatype = "pt_fastmtt/F";
      std::string branch_name_eta_datatype = "eta_fastmtt/F";
      std::string branch_name_phi_datatype = "phi_fastmtt/F";
      std::string branch_name_m_datatype = "m_fastmtt/F";
      svfitfriend->Branch(branch_name_pt.c_str(), &pt_fastmtt,
                          branch_name_pt_datatype.c_str());
      svfitfriend->Branch(branch_name_eta.c_str(), &eta_fastmtt,
                          branch_name_eta_datatype.c_str());
      svfitfriend->Branch(branch_name_phi.c_str(), &phi_fastmtt,
                          branch_name_phi_datatype.c_str());
      svfitfriend->Branch(branch_name_m.c_str(), &m_fastmtt,
                          branch_name_m_datatype.c_str());
    } else {
      std::string branch_name_pt =
          "pt_fastmtt" + std::string("__") + folder.c_str();
      std::string branch_name_eta =
          "eta_fastmtt" + std::string("__") + folder.c_str();
      std::string branch_name_phi =
          "phi_fastmtt" + std::string("__") + folder.c_str();
      std::string branch_name_m =
          "m_fastmtt" + std::string("__") + folder.c_str();
      std::string branch_name_pt_datatype =
          "pt_fastmtt" + std::string("__") + folder.c_str() + std::string("/F");
      std::string branch_name_eta_datatype = "eta_fastmtt" + std::string("__") +
                                             folder.c_str() + std::string("/F");
      std::string branch_name_phi_datatype = "phi_fastmtt" + std::string("__") +
                                             folder.c_str() + std::string("/F");
      std::string branch_name_m_datatype =
          "m_fastmtt" + std::string("__") + folder.c_str() + std::string("/F");
      svfitfriend->Branch(branch_name_pt.c_str(), &pt_fastmtt,
                          branch_name_pt_datatype.c_str());
      svfitfriend->Branch(branch_name_eta.c_str(), &eta_fastmtt,
                          branch_name_eta_datatype.c_str());
      svfitfriend->Branch(branch_name_phi.c_str(), &phi_fastmtt,
                          branch_name_phi_datatype.c_str());
      svfitfriend->Branch(branch_name_m.c_str(), &m_fastmtt,
                          branch_name_m_datatype.c_str());
    }

    // Float_t
    // pt_fastmtt_puppi,eta_fastmtt_puppi,phi_fastmtt_puppi,m_fastmtt_puppi;
    // svfitfriend->Branch("pt_fastmtt_puppi",&pt_fastmtt_puppi,"pt_fastmtt_puppi/F");
    // svfitfriend->Branch("eta_fastmtt_puppi",&eta_fastmtt_puppi,"eta_fastmtt_puppi/F");
    // svfitfriend->Branch("phi_fastmtt_puppi",&phi_fastmtt_puppi,"phi_fastmtt_puppi/F");
    // svfitfriend->Branch("m_fastmtt_puppi",&m_fastmtt_puppi,"m_fastmtt_puppi/F");

    // Initialize SVFit settings
    std::pair<MeasuredTauLepton::kDecayType, MeasuredTauLepton::kDecayType>
        ditaudecay = folder_to_ditaudecay(channel);

    // Initialize FastMTT
    FastMTT aFastMTTAlgo;

    // Loop over desired events of the input tree & compute outputs
    for (int i = first_entry; i < last_entry + include_last_ev; i++) {

      if ((i - first_entry) % 1000 == 0)
        std::cout << " processing event " << i
                  << ", last event: " << last_entry + include_last_ev << " ["
                  << float(i - first_entry) /
                         (last_entry + include_last_ev - first_entry)
                  << " %]" << std::endl;
      int bytes_read = inputtree->GetEntry(i);
      if (bytes_read == -1) {
        std::string message =
            "I/O error while reading event " + std::to_string(i) + ".";
        throw std::runtime_error(message);
      }

      // define MET;
      TVector2 metVec;
      metVec.SetMagPhi(met, metphi);

      // define MET covariance
      TMatrixD covMET(2, 2);
      covMET[0][0] = metcov00;
      covMET[1][0] = metcov10;
      covMET[0][1] = metcov01;
      covMET[1][1] = metcov11;

      // // define puppi MET;
      // TVector2 puppimetVec;
      // puppimetVec.SetMagPhi(puppimet,puppimetphi);

      // // define puppi MET covariance
      // TMatrixD puppicovMET(2, 2);
      // puppicovMET[0][0] = puppimetcov00;
      // puppicovMET[1][0] = puppimetcov10;
      // puppicovMET[0][1] = puppimetcov01;
      // puppicovMET[1][1] = puppimetcov11;

      // determine the right mass convention for the TauLepton decay products
      Float_t mass_1, mass_2;
      if (ditaudecay.first == MeasuredTauLepton::kTauToElecDecay)
        mass_1 = 0.51100e-3;
      else if (ditaudecay.first == MeasuredTauLepton::kTauToElecDecay)
        mass_1 = 105.658e-3;
      else
        mass_1 = m_1;

      if (ditaudecay.second == MeasuredTauLepton::kTauToElecDecay)
        mass_2 = 0.51100e-3;
      else if (ditaudecay.second == MeasuredTauLepton::kTauToElecDecay)
        mass_2 = 105.658e-3;
      else
        mass_2 = m_2;

      // define lepton four vectors
      std::vector<MeasuredTauLepton> measuredTauLeptons;
      if (ditaudecay.second == MeasuredTauLepton::kTauToHadDecay) {
        measuredTauLeptons.push_back(
            MeasuredTauLepton(ditaudecay.first, pt_1, eta_1, phi_1, mass_1,
                              decayMode_1 >= 0 ? decayMode_1 : -1));
      } else {
        measuredTauLeptons.push_back(MeasuredTauLepton(
            ditaudecay.first, pt_1, eta_1, phi_1, mass_1, -1));
      }

      measuredTauLeptons.push_back(
          MeasuredTauLepton(ditaudecay.second, pt_2, eta_2, phi_2, mass_2,
                            decayMode_2 >= 0 ? decayMode_2 : -1));

      // Run FastMTT with puppi
      aFastMTTAlgo.run(measuredTauLeptons, metVec.X(), metVec.Y(), covMET);
      LorentzVector P4 = aFastMTTAlgo.getBestP4();
      pt_fastmtt = P4.Pt();
      eta_fastmtt = P4.Eta();
      phi_fastmtt = P4.Phi();
      m_fastmtt = P4.M();

      // std::cout << "pt: " << pt_fastmtt << " eta: " << eta_fastmtt << " phi:
      // "
      //           << phi_fastmtt << " m: " << m_fastmtt << std::endl;

      // // Run FastMTT with puppi
      // aFastMTTAlgo.run(measuredTauLeptons,  puppimetVec.X(), puppimetVec.Y(),
      // puppicovMET); LorentzVector puppittP4 = aFastMTTAlgo.getBestP4();
      // pt_fastmtt_puppi = puppittP4.Pt();
      // eta_fastmtt_puppi = puppittP4.Eta();
      // phi_fastmtt_puppi = puppittP4.Phi();
      // m_fastmtt_puppi = puppittP4.M();
      // Fill output tree
      svfitfriend->Fill();
    }

    // Fill output file
    std::cout << "Writing output file ..." << std::endl;
    svfitfriend->Write("", TObject::kOverwrite);
    out->Close();
    in->Close();
    std::cout << "done" << std::endl;

    return 0;
  }
}
