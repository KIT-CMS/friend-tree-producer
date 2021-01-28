#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooFunctor.h"

#include <math.h>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <iostream>
#include <vector>
#include <map>
#include <regex>

#include "HiggsAnalysis/friend-tree-producer/interface/HelperFunctions.h"

using boost::starts_with;
namespace po = boost::program_options;

class NLOReweightingWeightsProducer
{
public:
    NLOReweightingWeightsProducer(std::string inputfile,
                                  int era,
                                  unsigned int first_entry,
                                  unsigned int last_entry,
                                  std::string workspace,
                                  std::string outputfile,
                                  std::string treename,
                                  std::string channel,
                                  std::string pipelines);
    ~NLOReweightingWeightsProducer();
    void run();

private:
    Float_t get_quantity(const std::string&);

    TFile *_inputfile;
    int _era;
    unsigned int _first_entry;
    unsigned int _last_entry;
    std::string _outputname;
    std::string _treename;
    std::string _channel;
    std::string _pipelines;
    RooWorkspace *_workspace;
};

NLOReweightingWeightsProducer::NLOReweightingWeightsProducer(std::string inputfile,
                                                             int era,
                                                             unsigned int first_entry,
                                                             unsigned int last_entry,
                                                             std::string workspace,
                                                             std::string outputfile,
                                                             std::string treename,
                                                             std::string channel,
                                                             std::string pipelines) : _inputfile(TFile::Open(inputfile.c_str(), "read")),
                                                                                  _era(era),
                                                                                  _first_entry(first_entry),
                                                                                  _last_entry(last_entry),
                                                                                  _outputname(outputfile),
                                                                                  _treename(treename),
                                                                                  _channel(channel),
                                                                                  _pipelines(pipelines)
{
    auto workspacefile = TFile::Open(workspace.c_str(), "read");
    _workspace = (RooWorkspace *)workspacefile->Get("w");
    workspacefile->Close();
    auto nickname = filename_from_inputpath(inputfile);
}

Float_t NLOReweightingWeightsProducer::get_quantity(const std::string& func_name)
{
    return _workspace->function(func_name.c_str())->getVal();
}

void NLOReweightingWeightsProducer::run()
{
    auto nickname = filename_from_inputpath(_inputfile->GetName());
    std::regex mass_regex("SUSYGluGluToHToTauTauM([0-9]+)_");
    std::smatch mass_match;
    if (std::regex_search(nickname, mass_match, mass_regex)) {
        for (size_t i = 0; i < mass_match.size(); i++) {
            std::cout << i << ": " << mass_match[i] << "\n";
        }
    }
    else {
        std::cout << "No match of expression in nickname " << nickname << std::endl;
    }
    std::string mass = mass_match[1];
    // maybe add loop over pipelines
    auto dir = (TDirectoryFile *)_inputfile->Get(_pipelines.c_str());
    auto inputtree = (TTree *)dir->Get(_treename.c_str());

    auto out = TFile::Open(_outputname.c_str(), "recreate");
    out->mkdir(_pipelines.c_str());
    out->cd(_pipelines.c_str());
    //
    // Create output tree
    auto output_tree = new TTree(_treename.c_str(), _treename.c_str());
    
    // Set up branches to be read out.
    std::vector<std::string> higgs_bosons = {"h", "A"};
    std::vector<std::string> contributions = {"t", "b", "i"};
    std::vector<std::string> variations = {"", "_scale_up", "_scale_down", "_hdamp_up", "_hdamp_down"};
    std::map<std::string, std::string> branches;
    for (auto var: variations) {
        for (auto contrib: contributions) {
            for (auto h_boson: higgs_bosons) {
                branches["gg" + h_boson + "_" + contrib + "_weight" + var] = h_boson + "_" + contrib + "_ratio" + var;
            }
        }
    }

    std::map<std::string, Float_t> output_buffer;
    for (auto branch : branches)
    {
        output_buffer[branch.first] = 0.0;
        output_tree->Branch(branch.first.c_str(), &(output_buffer[branch.first]), (branch.first + "/F").c_str());
    }
    // Different readout for embedded samples.
    // Set up read out of quantities from event.
    Float_t higgs_pt;
    inputtree->SetBranchAddress("genbosonpt", &(higgs_pt));

    // Start event loop here.
    for (unsigned int i = _first_entry; i <= _last_entry; i++)
    {
        if ((i - _first_entry) % 1000 == 0)
            std::cout << " processing event " << i << ", last event: " << _last_entry << " [" << 100 * (float(i - _first_entry) / (_last_entry - _first_entry)) << " %]" << std::endl;
        // Get entry
        inputtree->GetEntry(i);
        _workspace->var("h_pt")->setVal(higgs_pt);
        _workspace->var("h_mass")->setVal(std::stof(mass));
        for (auto branch : branches)
        {
            output_buffer[branch.first] = get_quantity(branch.second);
        }
        output_tree->Fill();
    }
    output_tree->Write();
    std::cout << "Tree successfully written." << std::endl;
    out->Close();
}

NLOReweightingWeightsProducer::~NLOReweightingWeightsProducer()
{
}

int main(int argc, char **argv)
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

    std::string cmsswbase = std::getenv("CMSSW_BASE");
    // Get channel from pipeline name
    std::vector<std::string> fields;
    boost::split(fields, folder, boost::is_any_of("_"));
    std::string channel = fields[0];

    auto nickname = filename_from_inputpath(input);
    boost::property_tree::ptree pt;
    boost::property_tree::read_json(cmsswbase + "/src/HiggsAnalysis/friend-tree-producer/data/input_params/datasets.json", pt);
    int year = pt.get<int>(nickname + ".year");
    
    std::string workspacepath;
    // 2017 and 2018 share the same input workspace
    if (year == 2016) {
        workspacepath = cmsswbase + "/src/HiggsAnalysis/friend-tree-producer/data/input_params/higgs_pt_" + std::to_string(year) + "_v0.root";
    }
    else if (year == 2017 || year == 2018) {
        workspacepath = cmsswbase + "/src/HiggsAnalysis/friend-tree-producer/data/input_params/higgs_pt_v0.root";
    }
    else {
        std::cout << "Era " << year << " of the sample is not supported. Aborting..." << std::endl;
        assert(0);
    }

    std::cout << "Using " << workspacepath << std::endl;
    auto outputname =
        outputname_from_settings(input, folder, first_entry, last_entry);
    boost::filesystem::create_directories(filename_from_inputpath(input));

    NLOReweightingWeightsProducer producer = NLOReweightingWeightsProducer(input,
                                                                           year,
                                                                           first_entry,
                                                                           last_entry,
                                                                           workspacepath,
                                                                           outputname,
                                                                           tree,
                                                                           channel,
                                                                           folder);
    producer.run();

    return 0;
}
