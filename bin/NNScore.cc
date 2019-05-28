#include "lwtnn/LightweightNeuralNetwork.hh"
#include "lwtnn/parse_json.hh"
#include <fstream>

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TVector2.h"
#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/PxPyPzM4D.h"
#include "Math/Vector4Dfwd.h"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>

#include <iostream>

#include "HiggsAnalysis/friend-tree-producer/interface/HelperFunctions.h"

using boost::starts_with;
namespace po = boost::program_options;

int main(int argc, char **argv) {
  std::string input = "output.root";
  std::vector<std::string> input_friends  = {};
  std::string folder = "mt_nominal";
  std::string tree = "ntuple";
  std::string lwtnn_config = "HiggsAnalysis/friend-tree-producer/data/inputs_lwtnn/";
  unsigned int first_entry = 0;
  unsigned int last_entry = 9;
  po::variables_map vm;
  po::options_description config("configuration");
  config.add_options()
     ("input",         po::value<std::string>(&input)->default_value(input))
     ("input_friends", po::value<std::vector<std::string>>(&input_friends)->multitoken())
     ("folder",        po::value<std::string>(&folder)->default_value(folder))
     ("tree",          po::value<std::string>(&tree)->default_value(tree))
     ("first_entry",   po::value<unsigned int>(&first_entry)->default_value(first_entry))
     ("last_entry",    po::value<unsigned int>(&last_entry)->default_value(last_entry))
     ("lwtnn_config",  po::value<std::string>(&lwtnn_config)->default_value(lwtnn_config));
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);
  // Add additional info inferred from options above
  std::vector<std::string> folder_split;
  boost::split(folder_split, folder, boost::is_any_of("_"));
  std::string channel = folder_split.at(0);

  // Access input file and tree
  auto in = TFile::Open(input.c_str(), "read");
  auto dir = (TDirectoryFile *)in->Get(folder.c_str());
  auto inputtree = (TTree *)dir->Get(tree.c_str());
  for(std::vector<std::string>::const_iterator s = input_friends.begin(); s != input_friends.end(); ++s)
  {
    inputtree->AddFriend((folder+"/"+tree).c_str(), (&(*s))->c_str());
  }

  // Set up lwtnn
  if (!boost::filesystem::exists(lwtnn_config)) {
      throw std::runtime_error("LWTNN config file does not exist.");
  }
  std::map<int,lwt::LightweightNeuralNetwork*> models;

  std::ifstream config_file0(lwtnn_config+"/fold0_lwtnn.json");
  auto nnconfig0 = lwt::parse_json(config_file0);
  models[0] = new lwt::LightweightNeuralNetwork(nnconfig0.inputs, nnconfig0.layers, nnconfig0.outputs);

  std::ifstream config_file1(lwtnn_config+"/fold1_lwtnn.json");
  auto nnconfig1 = lwt::parse_json(config_file1);
  models[1] = new lwt::LightweightNeuralNetwork(nnconfig1.inputs, nnconfig1.layers, nnconfig1.outputs);

  // Initialize inputs
  std::map<std::string, Float_t> inputs;
  for(size_t n=0; n < nnconfig0.inputs.size(); n++)
  {
    inputs[nnconfig0.inputs.at(n).name] = 0.0;
    inputtree->SetBranchAddress((nnconfig0.inputs.at(n).name).c_str(), &(inputs.find(nnconfig0.inputs.at(n).name)->second));
  }
  ULong64_t event;
  inputtree->SetBranchAddress("event", &event);

  // Initialize output file
  auto outputname =
      outputname_from_settings(input, folder, first_entry, last_entry);
  boost::filesystem::create_directories(filename_from_inputpath(input));
  auto out = TFile::Open(outputname.c_str(), "recreate");
  out->mkdir(folder.c_str());
  out->cd(folder.c_str());

  // Create output tree
  auto nnfriend = new TTree("ntuple", "NN score friend tree");

  // Initialize outputs for the tree
  std::map<std::string, Float_t> outputs;
  for(size_t n=0; n < nnconfig0.outputs.size(); n++)
  {
    outputs[nnconfig0.outputs.at(n)] = 0.0;
    nnfriend->Branch((channel+"_"+nnconfig0.outputs.at(n)).c_str(), &(outputs.find(nnconfig0.outputs.at(n))->second), (channel+"_"+nnconfig0.outputs.at(n)+"/F").c_str());
  }
  Float_t max_score = default_float;
  Float_t max_index = 0.0;
  std::string max_score_name = channel+"_max_score";
  std::string max_index_name = channel+"_max_index";
  nnfriend->Branch(max_score_name.c_str(), &max_score, (max_score_name+"/F").c_str());
  nnfriend->Branch(max_index_name.c_str(), &max_index, (max_index_name+"/F").c_str());

  // Loop over desired events of the input tree & compute outputs
  for (unsigned int i = first_entry; i <= last_entry; i++) {
    // Get entry
    inputtree->GetEntry(i);
    std::cout << "Event: " <<  event << std::endl;

    // Convert the inputs from Float_t to double
    std::map<std::string, double> model_inputs;
    for(auto &in : inputs)
    {
      model_inputs[in.first] = in.second;
    }
   
    // Apply model on inputs
    auto model_outputs = models[event % 2]->compute(model_inputs);

    // Fill output map
    for(size_t index=0; index < nnconfig0.outputs.size(); index++)
    {
      auto output_value = model_outputs[nnconfig0.outputs.at(index)];
      outputs[nnconfig0.outputs.at(index)] = output_value;
      if (output_value > max_score)
      {
        max_score = output_value;
        max_index = index;
      }
    }

    // Fill output tree
    nnfriend->Fill();

    // Reset max quantities
    max_score = default_float;
    max_index = 0.0;
  }

  // Fill output file
  out->cd(folder.c_str());
  nnfriend->Write("", TObject::kOverwrite);
  out->Close();
  in->Close();

  return 0;
}
