#include "lwtnn/LightweightNeuralNetwork.hh"
#include "lwtnn/parse_json.hh"
#include <fstream>

#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/PxPyPzM4D.h"
#include "Math/Vector4Dfwd.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TVector2.h"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/regex.hpp>

#include <iostream>

#include "HiggsAnalysis/friend-tree-producer/interface/HelperFunctions.h"

using boost::starts_with;
namespace po = boost::program_options;

int main(int argc, char **argv) {
  std::string input = "output.root";
  std::vector<std::string> input_friends = {};
  std::string folder = "mt_nominal";
  std::string tree = "ntuple";
  std::string channel = "mt";
  std::string era = "2018";
  std::string trainings_folder = "";
  unsigned int first_entry = 0;
  unsigned int last_entry = 9;
  bool run_required = false;
  po::variables_map vm;
  po::options_description config("configuration");
  config.add_options()("input",
                       po::value<std::string>(&input)->default_value(input))(
      "input_friends",
      po::value<std::vector<std::string>>(&input_friends)->multitoken())(
      "folder", po::value<std::string>(&folder)->default_value(folder))(
      "tree", po::value<std::string>(&tree)->default_value(tree))(
      "first_entry",
      po::value<unsigned int>(&first_entry)->default_value(first_entry))(
      "last_entry",
      po::value<unsigned int>(&last_entry)->default_value(last_entry))(
      "trainings_folder", po::value<std::string>(&trainings_folder)
                              ->default_value(trainings_folder))(
      "channel", po::value<std::string>(&channel)->default_value(channel))(
      "era", po::value<std::string>(&era)->default_value(era));
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);

  // Access input file and tree
  std::cout << "Access input file and tree" << std::endl;
  TFile *in = TFile::Open(input.c_str(), "read");
  TTree *inputtree = (TTree *)in->Get(tree.c_str());
  std::map<std::string, TTree *> input_friends_trees;
  // now loop over all input friends and add them to the map
  for (auto &input_friend : input_friends) {
    std::cout << "Access input friend file and tree" << std::endl;
    TFile *in_friend = TFile::Open(input_friend.c_str(), "read");
    TTree *inputtree_friend = (TTree *)in_friend->Get(tree.c_str());
    input_friends_trees[input_friend] = inputtree_friend;
    inputtree->AddFriend(input_friends_trees[input_friend]);
  }

  // for (std::vector<std::string>::const_iterator s = input_friends.begin();
  //      s != input_friends.end(); ++s) {
  //   auto f = std::unique_ptr<TFile>(TFile::Open((&(*s))->c_str()));
  //   input_friends_trees[&(*s)] = f->Get<TTree>(tree.c_str());
  //   std::cout << "Adding friend tree " << *s << " from tree " << tree
  //             << std::endl;

  //   inputtree->AddFriend(input_friends_trees[&(*s)]);
  // }
  // print all quantities in the tree
  inputtree->GetListOfFriends()->Print();
  std::cout << "Systematic: " << folder << std::endl;
  // Set up lwtnn
  std::cout << "Set up lwtnn from " << trainings_folder << std::endl;
  if (!boost::filesystem::exists(trainings_folder)) {
    throw std::runtime_error("LWTNN config file does not exist.");
  }
  std::map<int, lwt::LightweightNeuralNetwork *> models;

  std::ifstream config_file0(trainings_folder + "/fold0_lwtnn.json");
  auto nnconfig0 = lwt::parse_json(config_file0);
  models[1] = new lwt::LightweightNeuralNetwork(
      nnconfig0.inputs, nnconfig0.layers, nnconfig0.outputs);
  std::cout
      << "Loading fold0 model for application on ODD events (event % 2 == 1)"
      << std::endl;

  std::ifstream config_file1(trainings_folder + "/fold1_lwtnn.json");
  auto nnconfig1 = lwt::parse_json(config_file1);
  models[0] = new lwt::LightweightNeuralNetwork(
      nnconfig1.inputs, nnconfig1.layers, nnconfig1.outputs);
  std::cout
      << "Loading fold1 model for application on EVEN events (event % 2 == 0)"
      << std::endl;

  // Initialize inputs
  std::cout << "Initialize input" << std::endl;
  // get a vector of all the input variables from the nnconfig0 file
  std::vector<std::string> required_quantities;
  for (auto &input : nnconfig0.inputs) {
    required_quantities.push_back(input.name);
  }
  std::map<std::string, Float_t> float_inputs;
  std::map<std::string, Double_t> double_inputs;
  std::map<std::string, Int_t> int_inputs;
  size_t input_size = nnconfig0.inputs.size();
  std::vector<std::string> available_qunatities;
  TObjArray *leavescopy = inputtree->GetListOfLeaves();
  for (int i = 0; i < leavescopy->GetEntries(); ++i) {
    available_qunatities.push_back(leavescopy->At(i)->GetName());
  }
  // add all the quantites from the friends as well
  for (auto &input_friend : input_friends) {
    TObjArray *leavescopy_friend =
        input_friends_trees[input_friend]->GetListOfLeaves();
    for (int i = 0; i < leavescopy_friend->GetEntries(); ++i) {
      available_qunatities.push_back(leavescopy_friend->At(i)->GetName());
    }
  }
  std::vector<std::string> quantities =
      find_quantities(available_qunatities, required_quantities);
  std::map<std::string, std::string> quantity_names = build_quantities_map(
      folder, quantities, required_quantities, run_required, false);
  // print the quantities map
  std::cout << "Quantities map:" << std::endl;
  for (auto &kv : quantity_names) {
    std::cout << "    " << kv.first << " -> " << kv.second << std::endl;
  }
  // sanity check, the len of the quantities map should be the same as the input
  // size
  if (quantity_names.size() != input_size) {
    // print the quantities that are missing
    std::cout << "Missing quantities: " << std::endl;
    for (auto &input : nnconfig0.inputs) {
      if (quantity_names.find(input.name) == quantity_names.end()) {
        std::cout << input.name << std::endl;
      }
    }
    throw std::runtime_error("The number of input variables does not match the "
                             "number of inputs in the nnconfig0 file.");
  }
  if (!run_required) {
    std::cout << "No Execution needed, aborting ..";
    return 0;
  }

  inputtree->SetBranchStatus("*",0);
  for (size_t n = 0; n < input_size; n++) {
    std::string unshifted_name = nnconfig0.inputs.at(n).name;
    std::string quantity_name = quantity_names[unshifted_name];
    std::string input_type =
        inputtree->GetLeaf(quantity_name.c_str())->GetTypeName();
    inputtree->SetBranchStatus(quantity_name.c_str(), 1);
    std::cout << "Input " << n << ": " << quantity_name << " (" << input_type
              << ")" << std::endl;
    if (input_type == "Float_t") {
        float_inputs[unshifted_name] = -10.0;
        inputtree->SetBranchAddress(
            quantity_name.c_str(),
            &(float_inputs.find(unshifted_name)->second));
    } else if (input_type == "Int_t") {
      int_inputs[unshifted_name] = 0;
      inputtree->SetBranchAddress(quantity_name.c_str(),
                                  &(int_inputs.find(unshifted_name)->second));
    } else if (input_type == "Double_t") {
      double_inputs[unshifted_name] = 0.0;
      inputtree->SetBranchAddress(
          quantity_name.c_str(), &(double_inputs.find(unshifted_name)->second));
    } else {
      std::cout << "Type " << input_type
                << " not implemented! Exiting with error code '1'" << std::endl;
      return 1;
    }
  }
  ULong64_t event;
  inputtree->SetBranchAddress("event", &event);
  std::cout << "Inputs intialized" << std::endl;
  // Initialize output file
  std::cout << "Initializing output file" << std::endl;
  std::cout << "Passing args: " << input << " " << folder << " " << era << " "
            << channel << " " << first_entry << " " << last_entry << " "
            << std::endl;
  std::string outputname = outputname_from_settings_crown(
      input, folder, first_entry, last_entry, "", era, channel, true);
  std::cout << "Output file: " << outputname << std::endl;
  std::cout << "Creating folder: " << filename_from_inputpath(input)
            << std::endl;
  boost::filesystem::create_directories(era + "_" + channel + "_" +
                                        filename_from_inputpath(input));
  TFile *out = TFile::Open(outputname.c_str(), "recreate");

  // Create output tree
  auto nnfriend = new TTree("ntuple", "NN score friend tree");

  // Initialize outputs for the tree
  std::map<std::string, Float_t> outputs;
  for (size_t n = 0; n < nnconfig0.outputs.size(); n++) {
    outputs[nnconfig0.outputs.at(n)] = 0.0;
    if (folder == "nominal") {
      std::string output_name = channel + "_" + nnconfig0.outputs.at(n);
      nnfriend->Branch((output_name).c_str(),
                       &(outputs.find(nnconfig0.outputs.at(n))->second),
                       (output_name + "/F").c_str());
    } else {
      std::string output_name =
          channel + "_" + nnconfig0.outputs.at(n) + "__" + folder;
      nnfriend->Branch((output_name).c_str(),
                       &(outputs.find(nnconfig0.outputs.at(n))->second),
                       (output_name + "/F").c_str());
    }
  }
  Float_t max_score = default_float;
  Float_t max_index = 0.0;
  std::string max_score_name, max_index_name;
  if (folder == "nominal") {
    max_score_name = channel + "_max_score";
    max_index_name = channel + "_max_index";
    nnfriend->Branch(max_score_name.c_str(), &max_score,
                     (max_score_name + "/F").c_str());
    nnfriend->Branch(max_index_name.c_str(), &max_index,
                     (max_index_name + "/F").c_str());
  } else {
    max_score_name = channel + "_max_score__" + folder;
    max_index_name = channel + "_max_index__" + folder;
    nnfriend->Branch(max_score_name.c_str(), &max_score,
                     (max_score_name + "/F").c_str());
    nnfriend->Branch(max_index_name.c_str(), &max_index,
                     (max_index_name + "/F").c_str());
  }

  std::cout << "Setup done" << std::endl;
  // Loop over desired events of the input tree & compute outputs
  int total_events = last_entry - first_entry;
  // std::cout << "Starting processing of " << total_events << " events"
  //           << std::endl;
  // std::cout << "first entry: " << first_entry << " last entry: " <<
  // last_entry
  //           << std::endl;
  for (unsigned int i = first_entry; i <= last_entry; i++) {
    // periodically print progress every 5000 events
    if ((i - first_entry) % 5000 == 0) {
      std::cout << "Processing event " << i - first_entry << " of "
                << total_events << " ("
                << 100 * (i - first_entry) / (total_events) << "%)"
                << std::endl;
    }
    // Get entry
    inputtree->GetEntry(i);

    // Convert the inputs from Float_t to double
    std::map<std::string, double> model_inputs;
    for (auto &in : float_inputs) {
      model_inputs[in.first] = in.second;
    }
    for (auto &in : int_inputs) {
      model_inputs[in.first] = in.second;
    }
    for (auto &in : double_inputs) {
      model_inputs[in.first] = in.second;
    }
    // print the model_inputs
    // std::cout << "Model inputs: " << std::endl;
    // for (auto &input : model_inputs) {
    //   std::cout << "    " << input.first << " " << input.second << std::endl;
    // }

    // Apply model on inputs
    auto model_outputs = models[event % 2]->compute(model_inputs);

    // Fill output map
    for (size_t index = 0; index < nnconfig0.outputs.size(); index++) {
      auto output_value = model_outputs[nnconfig0.outputs.at(index)];
      outputs[nnconfig0.outputs.at(index)] = output_value;
      if (output_value > max_score) {
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
  std::cout << "Outputs Computed" << std::endl;
  // Fill output file
  nnfriend->Write("", TObject::kOverwrite);
  out->Close();
  in->Close();

  return 0;
}