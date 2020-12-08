#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TVector2.h"
#include "RooWorkspace.h"
#include "RooAbsReal.h"

#include <math.h>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>

#include <iostream>
#include <vector>
#include <map>

#include "HiggsAnalysis/friend-tree-producer/interface/HelperFunctions.h"

using boost::starts_with;
namespace po = boost::program_options;

class SingletauTriggerWeightsProducer
{
    public:
        SingletauTriggerWeightsProducer(std::string inputfile,
                                        int era,
                                        unsigned int first_entry,
                                        unsigned int last_entry,
                                        std::string workspace,
                                        std::string outputfile,
                                        std::string config,
                                        std::string treename,
                                        std::string channel,
                                        std::string pipelines,
                                        std::string friend_file_name);
        ~SingletauTriggerWeightsProducer();
        void run();
        Float_t calculate_mc_med_region(std::map<std::string, Float_t>&);
        Float_t calculate_mc_high_region(std::map<std::string, Float_t>&);
        Float_t calculate_data_med_region(std::map<std::string, Float_t>&, bool is_emb=false);
        Float_t calculate_data_high_region(std::map<std::string, Float_t>&, bool is_emb=false);
        Float_t calculate_emb_med_region(std::map<std::string, Float_t>&);
        Float_t calculate_emb_high_region(std::map<std::string, Float_t>&);

    private:
        Float_t get_quantity(const std::string&, const std::vector<std::string>&,
                              std::map<std::string, std::string>&);

        TFile* _inputfile;
        int _era;
        unsigned int _first_entry;
        unsigned int _last_entry;
        std::string _outputname;
        std::string _treename;
        std::string _channel;
        std::string _pipelines;
        std::string _friend_file;
        RooWorkspace* _workspace;
        boost::property_tree::ptree _config;

        std::map<std::string, Float_t> _float_inputs = {
            {"pt_1", 0.},
            {"pt_2", 0.},
            {"eta_1", 0.},
            {"eta_2", 0.},
            {"crossTriggerDataEfficiencyWeight_medium_DeepTau_1", 0.},
            {"crossTriggerDataEfficiencyWeight_medium_DeepTau_2", 0.},
        };
        std::map<std::string, Int_t> _int_inputs = {
            {"decayMode_1", 0},
            {"decayMode_2", 0}
        };
};

SingletauTriggerWeightsProducer::SingletauTriggerWeightsProducer(std::string inputfile,
                                                                 int era,
                                                                 unsigned int first_entry,
                                                                 unsigned int last_entry,
                                                                 std::string workspace,
                                                                 std::string outputfile,
                                                                 std::string config,
                                                                 std::string treename,
                                                                 std::string channel,
                                                                 std::string pipelines,
                                                                 std::string friend_file_name
                                                                 ) :
    _inputfile(TFile::Open(inputfile.c_str(), "read")),
    _era(era),
    _first_entry(first_entry),
    _last_entry(last_entry),
    _outputname(outputfile),
    _treename(treename),
    _channel(channel),
    _pipelines(pipelines),
    _friend_file(friend_file_name)
{
        auto workspacefile = TFile::Open(workspace.c_str(), "read");
        _workspace = (RooWorkspace*)workspacefile->Get("w");
        workspacefile->Close();
        boost::property_tree::read_json(config, _config);
        auto nickname = filename_from_inputpath(inputfile);
        bool is_embedded = (nickname.find("Embedding") != std::string::npos);
        if (is_embedded)
        {
            _float_inputs["crossTriggerCorrectedMCEfficiencyWeight_medium_DeepTau_1"] = 0.;
            _float_inputs["crossTriggerCorrectedMCEfficiencyWeight_medium_DeepTau_2"] = 0.;
        }
        else
        {
            _float_inputs["crossTriggerMCEfficiencyWeight_medium_DeepTau_1"] = 0.;
            _float_inputs["crossTriggerMCEfficiencyWeight_medium_DeepTau_2"] = 0.;
        }
}

Float_t SingletauTriggerWeightsProducer::get_quantity(const std::string& function, const std::vector<std::string>& arguments,
                                                       std::map<std::string, std::string>& variables_map)
{
    // Translate branch entry names to workspace argument names
    std::vector<std::string> wspace_args;
    for (auto arg: arguments)
    {
        wspace_args.push_back(variables_map[arg]);
    }
    // Build comma separated list of parameter names in workspace.
    std::ostringstream ss;
    std::copy(wspace_args.begin(), std::prev(wspace_args.end()), std::ostream_iterator<std::string>(ss, ","));
    ss << wspace_args.back();

    RooArgSet argset = _workspace->argSet(ss.str().c_str());
    RooAbsReal* roofunction = _workspace->function(function.c_str());

    for (auto par: arguments)
    {
        // std::cout << "Parameter is: " << par << std::endl;
        // std::cout << "Mapped parameter is: " << variables_map[par] << std::endl;
        if (par.find("decayMode") != std::string::npos)
        {
            argset.setRealValue(variables_map[par].c_str(), _int_inputs[par]);
        }
        else
        {
            argset.setRealValue(variables_map[par].c_str(), _float_inputs[par]);
        }
    }
    Float_t result = roofunction->getVal(argset);
    return result;
}

void SingletauTriggerWeightsProducer::run()
{
  auto nickname = filename_from_inputpath(_inputfile->GetName());
  // maybe add loop over pipelines
  auto dir = (TDirectoryFile *)_inputfile->Get(_pipelines.c_str());
  auto inputtree = (TTree *)dir->Get(_treename.c_str());

  if (_friend_file != "None")
  {
      inputtree->AddFriend((_pipelines + "/" + _treename).c_str(),
                           _friend_file.c_str());
  }

  auto out = TFile::Open(_outputname.c_str(), "recreate");
  out->mkdir(_pipelines.c_str());
  out->cd(_pipelines.c_str());
  //
  // Create output tree
  auto output_tree = new TTree(_treename.c_str(), _treename.c_str());

  bool is_embedded = (nickname.find("Embedding") != std::string::npos);
  std::vector<std::string> additional_variables;
  if (is_embedded)
  {
        additional_variables.push_back("singleTauTriggerEMBWeightHighPt");
        additional_variables.push_back("singleTauTriggerEMBWeightMedPt");
        additional_variables.push_back("singleTauTriggerDataWeightHighPt");
        additional_variables.push_back("singleTauTriggerDataWeightMedPt");
        additional_variables.push_back("singleTauTriggerSFWeightHighPt");
        additional_variables.push_back("singleTauTriggerSFWeightMedPt");
  }
  else
  {
        additional_variables.push_back("singleTauTriggerMCWeightHighPt");
        additional_variables.push_back("singleTauTriggerMCWeightMedPt");
        additional_variables.push_back("singleTauTriggerDataWeightHighPt");
        additional_variables.push_back("singleTauTriggerDataWeightMedPt");
        additional_variables.push_back("singleTauTriggerSFWeightHighPt");
        additional_variables.push_back("singleTauTriggerSFWeightMedPt");
  }

  // Create map from variables to workspace object names
  std::map<std::string, std::string> var_map;
  BOOST_FOREACH(const boost::property_tree::ptree::value_type& child,
                _config.get_child("map_arguments"))
  {
      var_map[child.first] = child.second.get_value<std::string>();
  }
  // for (auto trans: var_map)
  //     std::cout << "Association: " << trans.first << " -> " << trans.second << std::endl;
  // create output buffer objects
  std::vector<std::string> branches;
  BOOST_FOREACH(const boost::property_tree::ptree::value_type& child,
                _config.get_child("rooworkspace"))
  {
      branches.push_back(child.first);
  }
  std::map<std::string, Float_t> output_buffer;
  for (auto branch: branches)
  {
      output_buffer[branch] = 0.0;
      output_tree->Branch(branch.c_str(), &(output_buffer[branch]), (branch + "/F").c_str());
  }
  for (auto branch: additional_variables)
  {
      output_buffer[branch] = 0.0;
      output_tree->Branch(branch.c_str(), &(output_buffer[branch]), (branch + "/F").c_str());
  }
  // Different readout for embedded samples.
  // Set up read out of quantities from event.
  for (auto & br: _float_inputs)
  {
      inputtree->SetBranchAddress(br.first.c_str(), &(br.second));
  }
  for (auto & br: _int_inputs)
  {
      inputtree->SetBranchAddress(br.first.c_str(), &(br.second));
  }

  // Start event loop here.
  for (unsigned int i = _first_entry; i <= _last_entry; i++)
  {
    // Get entry
    inputtree->GetEntry(i);
    for (auto branch: branches)
    {
        std::vector<std::string> arguments;
        BOOST_FOREACH(const boost::property_tree::ptree::value_type& child,
                      _config.get_child("rooworkspace."+branch+".WorkspaceObjectArguments"))
        {
            arguments.push_back(child.second.get_value<std::string>());
        }
        // for (auto par: arguments)
        // {
        //     std::cout << par << std::endl;
        // }
        output_buffer[branch] = get_quantity(_config.get<std::string>("rooworkspace."+branch+".WorkspaceWeightNames"), arguments, var_map);
        // std::cout << "Returned get_quantity for branch " << branch << ": " << get_quantity(_config.get<std::string>("rooworkspace."+branch+".WorkspaceWeightNames"), arguments, var_map) << std::endl;
        // std::cout << output_buffer[branch] << std::endl;
        if (is_embedded)
        {
            output_buffer["singleTauTriggerEMBWeightHighPt"] = calculate_emb_high_region(output_buffer);
            output_buffer["singleTauTriggerEMBWeightMedPt"] = calculate_emb_high_region(output_buffer);
            output_buffer["singleTauTriggerDataWeightHighPt"] = calculate_data_high_region(output_buffer, true);
            output_buffer["singleTauTriggerDataWeightMedPt"] = calculate_data_med_region(output_buffer, true);
            output_buffer["singleTauTriggerSFWeightHighPt"] = output_buffer["singleTauTriggerDataWeightHighPt"]/output_buffer["singleTauTriggerEMBWeightHighPt"];
            output_buffer["singleTauTriggerSFWeightMedPt"] = output_buffer["singleTauTriggerDataWeightMedPt"]/output_buffer["singleTauTriggerEMBWeightMedPt"];
        }
        else
        {
            output_buffer["singleTauTriggerMCWeightHighPt"] = calculate_mc_high_region(output_buffer);
            output_buffer["singleTauTriggerMCWeightMedPt"] = calculate_mc_med_region(output_buffer);
            output_buffer["singleTauTriggerDataWeightHighPt"] = calculate_data_high_region(output_buffer);
            output_buffer["singleTauTriggerDataWeightMedPt"] = calculate_data_med_region(output_buffer);
            output_buffer["singleTauTriggerSFWeightHighPt"] = output_buffer["singleTauTriggerDataWeightHighPt"]/output_buffer["singleTauTriggerMCWeightHighPt"];
            output_buffer["singleTauTriggerSFWeightMedPt"] = output_buffer["singleTauTriggerDataWeightMedPt"]/output_buffer["singleTauTriggerMCWeightMedPt"];
        }
    }
    output_tree->Fill();
  }
  output_tree->Write();
  std::cout << "Tree successfully written." << std::endl;
  out->Close();
}

SingletauTriggerWeightsProducer::~SingletauTriggerWeightsProducer()
{
}

Float_t SingletauTriggerWeightsProducer::calculate_emb_med_region(std::map<std::string, Float_t>& outbuf)
{
    Float_t eff_ditau_1 = outbuf["t_trg_double_t_emb_1"];
    Float_t eff_ditau_2 = outbuf["t_trg_double_t_emb_2"];
    Float_t eff_stau_1 = outbuf["t_trg_single_t_emb_1"];
    Float_t result = (eff_ditau_1 * eff_ditau_2) + eff_stau_1 - (eff_ditau_2 * eff_stau_1);
    return result;
}


Float_t SingletauTriggerWeightsProducer::calculate_emb_high_region(std::map<std::string, Float_t>& outbuf)
{
    Float_t eff_ditau_1 = outbuf["t_trg_double_t_emb_1"];
    Float_t eff_ditau_2 = outbuf["t_trg_double_t_emb_2"];
    Float_t eff_stau_1 = outbuf["t_trg_single_t_emb_1"];
    Float_t eff_stau_2 = outbuf["t_trg_single_t_emb_2"];
    Float_t result = eff_ditau_1*eff_ditau_2 + eff_stau_1 + eff_stau_2 \
                     - eff_ditau_1*eff_stau_2 - eff_ditau_2*eff_stau_1;
    return result;
}


Float_t SingletauTriggerWeightsProducer::calculate_mc_med_region(std::map<std::string, Float_t>& outbuf)
{
    Float_t eff_ditau_1 = outbuf["t_trg_double_t_mc_1"];
    Float_t eff_ditau_2 = outbuf["t_trg_double_t_mc_2"];
    Float_t eff_stau_1 = outbuf["t_trg_single_t_mc_1"];
    Float_t eff_stau_and_ditau_1 = outbuf["t_trg_singleANDdouble_t_mc_1"];
    Float_t result = eff_ditau_1*eff_ditau_2 + eff_stau_1 - eff_ditau_2*eff_stau_and_ditau_1;
    return result;
}


Float_t SingletauTriggerWeightsProducer::calculate_mc_high_region(std::map<std::string, Float_t>& outbuf)
{
    Float_t eff_ditau_1 = outbuf["t_trg_double_t_mc_1"];
    Float_t eff_ditau_2 = outbuf["t_trg_double_t_mc_2"];
    Float_t eff_stau_1 = outbuf["t_trg_single_t_mc_1"];
    Float_t eff_stau_2 = outbuf["t_trg_single_t_mc_2"];
    Float_t eff_stau_and_ditau_1 = outbuf["t_trg_singleANDdouble_t_mc_1"];
    Float_t eff_stau_and_ditau_2 = outbuf["t_trg_singleANDdouble_t_mc_2"];
    Float_t result = eff_ditau_1*eff_ditau_2 + eff_stau_1 + eff_stau_2 \
             - eff_ditau_1*eff_stau_and_ditau_2 - eff_ditau_2*eff_stau_and_ditau_1 - eff_stau_1*eff_stau_2 \
             + eff_stau_and_ditau_1*eff_stau_and_ditau_2;
    return result;
}


Float_t SingletauTriggerWeightsProducer::calculate_data_med_region(std::map<std::string, Float_t>& outbuf, bool is_emb)
{
    Float_t eff_ditau_1 = _float_inputs["crossTriggerDataEfficiencyWeight_medium_DeepTau_1"];
    Float_t eff_ditau_2 = _float_inputs["crossTriggerDataEfficiencyWeight_medium_DeepTau_2"];
    Float_t sf_ditau_1, sf_ditau_2;
    if (is_emb)
    {
        sf_ditau_1 = eff_ditau_1/_float_inputs["crossTriggerCorrectedMCEfficiencyWeight_medium_DeepTau_1"];
        sf_ditau_2 = eff_ditau_2/_float_inputs["crossTriggerCorrectedMCEfficiencyWeight_medium_DeepTau_2"];
    }
    else
    {
        sf_ditau_1 = eff_ditau_1/_float_inputs["crossTriggerMCEfficiencyWeight_medium_DeepTau_1"];
        sf_ditau_2 = eff_ditau_2/_float_inputs["crossTriggerMCEfficiencyWeight_medium_DeepTau_2"];
    }
    Float_t sf_ditau = sf_ditau_1*sf_ditau_2;
    Float_t eff_ditau_1_mc = outbuf["t_trg_double_t_mc_1"];
    Float_t eff_ditau_2_mc = outbuf["t_trg_double_t_mc_2"];
    Float_t eff_stau_1 = outbuf["t_trg_single_t_wstar_sf_1"]*outbuf["t_trg_single_t_mc_1"];
    Float_t eff_stau_and_ditau_1 = outbuf["t_trg_single_t_wstar_sf_1"]*outbuf["t_trg_singleANDdouble_t_mc_1"];
    Float_t result = sf_ditau*eff_ditau_1_mc*eff_ditau_2_mc + eff_stau_1 - sf_ditau*eff_ditau_2_mc*eff_stau_and_ditau_1;
    return result;
}


Float_t SingletauTriggerWeightsProducer::calculate_data_high_region(std::map<std::string, Float_t>&  outbuf, bool is_emb)
{
    Float_t eff_ditau_1 = _float_inputs["crossTriggerDataEfficiencyWeight_medium_DeepTau_1"];
    Float_t eff_ditau_2 = _float_inputs["crossTriggerDataEfficiencyWeight_medium_DeepTau_2"];
    Float_t sf_ditau_1, sf_ditau_2;
    if (is_emb)
    {
        sf_ditau_1 = eff_ditau_1/_float_inputs["crossTriggerCorrectedMCEfficiencyWeight_medium_DeepTau_1"];
        sf_ditau_2 = eff_ditau_2/_float_inputs["crossTriggerCorrectedMCEfficiencyWeight_medium_DeepTau_2"];
    }
    else
    {
        sf_ditau_1 = eff_ditau_1/_float_inputs["crossTriggerMCEfficiencyWeight_medium_DeepTau_1"];
        sf_ditau_2 = eff_ditau_2/_float_inputs["crossTriggerMCEfficiencyWeight_medium_DeepTau_2"];
    }
    Float_t sf_ditau = sf_ditau_1*sf_ditau_2;
    Float_t eff_ditau_1_mc = outbuf["t_trg_double_t_mc_1"];
    Float_t eff_ditau_2_mc = outbuf["t_trg_double_t_mc_2"];
    Float_t eff_stau_1 = outbuf["t_trg_single_t_wstar_sf_1"]*outbuf["t_trg_single_t_mc_1"];
    Float_t eff_stau_2 = outbuf["t_trg_single_t_wstar_sf_2"]*outbuf["t_trg_single_t_mc_2"];
    Float_t eff_stau_and_ditau_1 = outbuf["t_trg_single_t_wstar_sf_1"]*outbuf["t_trg_singleANDdouble_t_mc_1"];
    Float_t eff_stau_and_ditau_2 = outbuf["t_trg_single_t_wstar_sf_2"]*outbuf["t_trg_singleANDdouble_t_mc_2"];
    Float_t result = sf_ditau*eff_ditau_1_mc*eff_ditau_2_mc + eff_stau_1 + eff_stau_2 \
             - sf_ditau*eff_ditau_1_mc*eff_stau_and_ditau_2 - sf_ditau*eff_ditau_2_mc*eff_stau_and_ditau_1 - eff_stau_1*eff_stau_2 \
             + sf_ditau*eff_ditau_1_mc*eff_ditau_2_mc*eff_stau_and_ditau_1*eff_stau_and_ditau_2;
    return result;
}

int main(int argc, char **argv) {
  std::string input = "output.root";
  std::string folder = "mt_nominal";
  std::string tree = "ntuple";
  std::string workspacepath = "workspace.root";
  unsigned int first_entry = 0;
  unsigned int last_entry = 9;
  std::string friend_file_name = "None";
  po::variables_map vm;
  po::options_description config("configuration");
  config.add_options()("input",
                       po::value<std::string>(&input)->default_value(input))(
      "folder", po::value<std::string>(&folder)->default_value(folder))(
      "tree", po::value<std::string>(&tree)->default_value(tree))(
      "first_entry",
      po::value<unsigned int>(&first_entry)->default_value(first_entry))(
      "last_entry",
      po::value<unsigned int>(&last_entry)->default_value(last_entry))(
      "rooworkspace-file",
      po::value<std::string>(&workspacepath)->default_value(workspacepath))(
      "input_friends",
      po::value<std::string>(&friend_file_name)->default_value(friend_file_name));
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);

  std::string cmsswbase = std::getenv("CMSSW_BASE");
   // Get channel from pipeline name
  std::vector<std::string> fields;
  boost::split(fields, folder, boost::is_any_of("_"));
  std::string channel =  fields[0];

  auto nickname = filename_from_inputpath(input);
  boost::property_tree::ptree pt;
  boost::property_tree::read_json(cmsswbase+"/src/HiggsAnalysis/friend-tree-producer/data/input_params/datasets.json", pt);
  int year = pt.get<int>(nickname + ".year");

  std::string configpath = cmsswbase + \
            "/src/HiggsAnalysis/friend-tree-producer/data/config_singletauTriggerWeightsProducer.json";

  auto outputname = \
      outputname_from_settings(input, folder, first_entry, last_entry);
  boost::filesystem::create_directories(filename_from_inputpath(input));

  SingletauTriggerWeightsProducer producer = SingletauTriggerWeightsProducer(input,
                                                                             year,
                                                                             first_entry,
                                                                             last_entry,
                                                                             workspacepath,
                                                                             outputname,
                                                                             configpath,
                                                                             tree,
                                                                             channel,
                                                                             folder,
                                                                             friend_file_name);
  producer.run();

  return 0;
}
