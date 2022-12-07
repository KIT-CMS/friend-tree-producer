#include "correction.h"

#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

#include "HiggsAnalysis/friend-tree-producer/interface/HelperFunctions.h"

namespace po = boost::program_options;


int main(int argc, char **argv) {
  std::string input = "output.root";
  std::string ff_file = "fake_factors.json";
  std::string output_dir = "";
  std::string folder = "mt_nominal";
  std::string tree = "ntuple";
  std::string channel = "mt";
  std::string era = "2018";
  int first_entry = 0;
  int last_entry = -1;
  bool organize_outputs = true;
  po::variables_map vm;
  po::options_description config("configuration");
  config.add_options()("input",
                       po::value<std::string>(&input)->default_value(input))(
      "ff_file",
      po::value<std::string>(&ff_file)->default_value(ff_file))(
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
  std::cout << "Setting up FakeFactors for channel " << channel << " and era "
            << era << " pipeline " << folder << std::endl;
  // Access input file and tree
  TFile *in = TFile::Open(input.c_str(), "read");
  // TDirectoryFile *dir = (TDirectoryFile *)in->Get(folder.c_str());
  TTree *inputtree = (TTree *)in->Get(tree.c_str());

  // Restrict input tree to needed branches
  inputtree->SetBranchStatus("*", 0);

  // Setting events processing ranges
  std::vector<std::string> required_quantities{
      "pt_2",     "mt_1",    "njets",    "nbtag",      "deltaR_ditaupair", "m_vis",
    };
  int include_last_ev = 1;
  bool run_required = true;
  if (last_entry < 0 || last_entry >= inputtree->GetEntries()) {
    last_entry = inputtree->GetEntries();
    include_last_ev = 0;
  }

  std::vector<std::string> quantities =
      find_quantities(inputtree, required_quantities);
  std::map<std::string, std::string> quantity_names = build_quantities_map(
      folder, quantities, required_quantities, run_required, false);

  // Initialize output file
  std::string outputname = outputname_from_settings_crown(
      input, folder, first_entry, last_entry, output_dir, era, channel,
      organize_outputs);
  if (organize_outputs) {
    boost::filesystem::create_directories(era + "_" + channel + "_" +
                                          filename_from_inputpath(input));
  }
  TFile *out = TFile::Open(outputname.c_str(), "recreate");
  std::cout << "Initializing Output file: " << outputname << std::endl;

  // Create output tree
  TTree *FFfriend = new TTree("ntuple", "fake factor friend tree");

  if (run_required == false) {
    std::cout << "No calculation needed, will write an empty tree" << std::endl;
    FFfriend->Write("", TObject::kOverwrite);
    out->Close();
    in->Close();
    std::cout << "done" << std::endl;
    return 0;
  } else {
    // Now we can set the branches to be used
    // Quantities of first lepton
    std::cout << "Setting branches for first lepton" << std::endl;

    inputtree->SetBranchStatus(quantity_names["pt_2"].c_str(), 1);
    inputtree->SetBranchStatus(quantity_names["mt_1"].c_str(), 1);
    inputtree->SetBranchStatus(quantity_names["njets"].c_str(), 1);
    inputtree->SetBranchStatus(quantity_names["nbtag"].c_str(), 1);
    inputtree->SetBranchStatus(quantity_names["deltaR_ditaupair"].c_str(), 1);
    Float_t pt_2, mt_1, m_vis;
    Double_t deltaR_ditaupair;
    Int_t njets, nbtag;
    inputtree->SetBranchAddress(quantity_names["pt_2"].c_str(), &pt_2);
    inputtree->SetBranchAddress(quantity_names["mt_1"].c_str(), &mt_1);
    inputtree->SetBranchAddress(quantity_names["njets"].c_str(), &njets);
    inputtree->SetBranchAddress(quantity_names["nbtag"].c_str(), &nbtag);
    inputtree->SetBranchAddress(quantity_names["deltaR_ditaupair"].c_str(), &deltaR_ditaupair);

    if (channel == "tt") {
      inputtree->SetBranchStatus("m_vis", 1);
      inputtree->SetBranchAddress("m_vis", &m_vis);
    }

    Float_t fake_factor;
    std::string branch_name_FF, branch_name_FF_datatype;

    if (folder == "nominal") {
      std::string branch_name_FF = "fake_factor";
      std::string branch_name_FF_datatype = "fake_factor/F";
      FFfriend->Branch(branch_name_FF.c_str(), &fake_factor,
                          branch_name_FF_datatype.c_str());
    } else {
      std::string branch_name_FF =
          "fake_factor" + std::string("__") + folder.c_str();
      std::string branch_name_FF_datatype =
          "fake_factor" + std::string("__") + folder.c_str() + std::string("/F");
      FFfriend->Branch(branch_name_FF.c_str(), &fake_factor,
                          branch_name_FF_datatype.c_str());
    }

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

    auto qcd = correction::CorrectionSet::from_file(ff_file)->at("QCD_fake_factors");
    auto wjets = correction::CorrectionSet::from_file(ff_file)->at("Wjets_fake_factors");
    auto ttbar = correction::CorrectionSet::from_file(ff_file)->at("ttbar_fake_factors");
    auto frac = correction::CorrectionSet::from_file(ff_file)->at("process_fractions");
    Float_t qcd_ff = 1.;
    Float_t wjets_ff = 1.;
    Float_t ttbar_ff = 1.;
    qcd_ff = qcd->evaluate({pt_2, static_cast<Float_t>(njets)});
    wjets_ff = wjets->evaluate({pt_2, static_cast<Float_t>(njets)});
    ttbar_ff = ttbar->evaluate({pt_2, static_cast<Float_t>(njets)});
    Float_t qcd_frac = 1.;
    Float_t wjets_frac = 1.;
    Float_t ttbar_frac = 1.;
    qcd_frac = frac->evaluate({"QCD", mt_1, static_cast<Float_t>(nbtag)});
    wjets_frac = frac->evaluate({"Wjets", mt_1, static_cast<Float_t>(nbtag)});
    ttbar_frac = frac->evaluate({"ttbar", mt_1, static_cast<Float_t>(nbtag)});
    fake_factor = qcd_frac*qcd_ff+wjets_frac*wjets_ff+ttbar_frac*ttbar_ff;

    FFfriend->Fill();
    }

    // Fill output file
    std::cout << "Writing output file ..." << std::endl;
    FFfriend->Write("", TObject::kOverwrite);
    out->Close();
    in->Close();
    std::cout << "done" << std::endl;

    return 0;
  }
}