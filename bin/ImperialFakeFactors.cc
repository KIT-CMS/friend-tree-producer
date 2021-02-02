#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TVector2.h"
#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "Math/Vector4D.h"
#include "Math/Vector4Dfwd.h"
#include "Math/VectorUtil.h"
#include "RooFunctor.h"

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

class ImperialFakeFactorsProducer
{
public:
    ImperialFakeFactorsProducer(std::string inputfile,
                                int era,
                                unsigned int first_entry,
                                unsigned int last_entry,
                                std::string workspace,
                                std::string outputfile,
                                std::string config,
                                std::string treename,
                                std::string channel,
                                std::string pipelines);
    ~ImperialFakeFactorsProducer();
    void run();
    Float_t calculate_met_var_qcd(std::map<std::string, Float_t> &);
    Float_t calculate_met_var_w(std::map<std::string, Float_t> &);
    Float_t calculate_os(std::map<std::string, Float_t> &);

private:
    Float_t get_quantity(const std::shared_ptr<RooFunctor> &, const std::vector<std::string> &,
                         std::map<std::string, std::string> &);

    TFile *_inputfile;
    int _era;
    unsigned int _first_entry;
    unsigned int _last_entry;
    std::string _outputname;
    std::string _treename;
    std::string _channel;
    std::string _pipelines;
    RooWorkspace *_workspace;
    boost::property_tree::ptree _config;

    std::map<std::string, Float_t> _float_inputs = {
        {"pt_1", 0.},
        {"pt_2", 0.},
        {"taujet_pt_1", 0.},
        {"taujet_pt_2", 0.},
        {"iso_1", 0.},
        {"eta_1", 0.},
        {"eta_2", 0.},
        {"phi_1", 0.},
        {"phi_2", 0.},
        {"puppimet", 0.},
        {"puppimetphi", 0.},
        {"mt_1_puppi", 0.},
        {"mt_2_puppi", 0.},
        {"DiTauDeltaR", 0.},
        {"q_1", 0.},
        {"q_2", 0.},
        {"mt_tot_puppi", 0.},};

    std::map<std::string, Double_t> _calc_inputs = {
        {"met_var_qcd", 0.},
        {"met_var_w", 0.},
        {"os", 0.}};
    std::map<std::string, Int_t> _int_inputs = {
        {"nbtag", 0},
        {"njetspt20eta2p4", 0},
        {"njetspt20eta2p5", 0}};
};

ImperialFakeFactorsProducer::ImperialFakeFactorsProducer(std::string inputfile,
                                                         int era,
                                                         unsigned int first_entry,
                                                         unsigned int last_entry,
                                                         std::string workspace,
                                                         std::string outputfile,
                                                         std::string config,
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
    boost::property_tree::read_json(config, _config);
    auto nickname = filename_from_inputpath(inputfile);
}

Float_t ImperialFakeFactorsProducer::get_quantity(const std::shared_ptr<RooFunctor> &function, const std::vector<std::string> &arguments, std::map<std::string, std::string> &variables_map)
{
    auto argvalues = std::vector<double>{};
    for (auto par : arguments)
    {
        // std::cout << "Setting  " << variables_map[par] << " for " << par << std::endl;
        if (par == "met_var_qcd")
        {
            argvalues.push_back(ImperialFakeFactorsProducer::calculate_met_var_qcd(_float_inputs));
        }
        else if (par == "met_var_w")
        {
            argvalues.push_back(ImperialFakeFactorsProducer::calculate_met_var_w(_float_inputs));
        }
        else if (par == "os")
        {
            argvalues.push_back(ImperialFakeFactorsProducer::calculate_os(_float_inputs));
        }
        else if (par == "nbjets" || par == "njets")
        {
            argvalues.push_back(_int_inputs[variables_map[par]]);
        }
        else
        {
            argvalues.push_back(_float_inputs[variables_map[par]]);
        }
    }
    // for (size_t index = 0; index < arguments.size(); ++index)
    //     std::cout << arguments[index] << " --> " << argvalues[index] << std::endl;
    Float_t result = function->eval(argvalues.data());
    return result;
}

void ImperialFakeFactorsProducer::run()
{
    auto nickname = filename_from_inputpath(_inputfile->GetName());
    // maybe add loop over pipelines
    auto dir = (TDirectoryFile *)_inputfile->Get(_pipelines.c_str());
    auto inputtree = (TTree *)dir->Get(_treename.c_str());

    auto out = TFile::Open(_outputname.c_str(), "recreate");
    out->mkdir(_pipelines.c_str());
    out->cd(_pipelines.c_str());
    //
    // Create output tree
    auto output_tree = new TTree(_treename.c_str(), _treename.c_str());

    // Create map from variables to workspace object names
    std::map<std::string, std::string> var_map;
    BOOST_FOREACH (const boost::property_tree::ptree::value_type &child,
                   _config.get_child("map_arguments." + std::to_string(_era) + "." + _channel))
    {
        var_map[child.first] = child.second.get_value<std::string>();
    }

    // create output buffer objects
    // get list of needed variables
    std::vector<std::string> arguments;
    BOOST_FOREACH (const boost::property_tree::ptree::value_type &child,
                   _config.get_child("rooworkspace." + _channel + ".arguments"))
    {
        arguments.push_back(child.second.get_value<std::string>());
    }
    // build varstring needed for RooFunctor
    std::ostringstream ss;
    std::copy(arguments.begin(), std::prev(arguments.end()), std::ostream_iterator<std::string>(ss, ","));
    ss << arguments.back();
    std::vector<std::string> branches;
    std::map<std::string, std::shared_ptr<RooFunctor>> fns_;
    BOOST_FOREACH (const boost::property_tree::ptree::value_type &child,
                   _config.get_child("rooworkspace." + _channel + ".functions.main"))
    {
        branches.push_back(child.first);
        if (child.first != "ff_total"){
            std::string toErase = ",mt_tot";
            std::string argumentlist = ss.str();
            size_t pos = argumentlist.find(toErase);
            argumentlist.erase(pos, toErase.length());
            fns_[child.first.c_str()] = std::shared_ptr<RooFunctor>(_workspace->function(child.first.c_str())->functor(_workspace->argSet(argumentlist.c_str())));
        }
        else{
            fns_[child.first.c_str()] = std::shared_ptr<RooFunctor>(_workspace->function(child.first.c_str())->functor(_workspace->argSet(ss.str().c_str())));
        }


    }
    BOOST_FOREACH (const boost::property_tree::ptree::value_type &child,
                   _config.get_child("rooworkspace." + _channel + ".functions.uncertainties"))
    {
        std::string upvariation = child.first + "_up";
        std::string downvariation = child.first + "_down";
        branches.push_back(upvariation);
        branches.push_back(downvariation);
        fns_[upvariation.c_str()] = std::shared_ptr<RooFunctor>(_workspace->function(upvariation.c_str())->functor(_workspace->argSet(ss.str().c_str())));
        fns_[downvariation.c_str()] = std::shared_ptr<RooFunctor>(_workspace->function(downvariation.c_str())->functor(_workspace->argSet(ss.str().c_str())));
    }
    std::map<std::string, Float_t> output_buffer;
    for (auto branch : branches)
    {
        output_buffer[branch] = 0.0;
        output_tree->Branch(branch.c_str(), &(output_buffer[branch]), (branch + "/F").c_str());
    }
    for (auto &branch : _calc_inputs)
    {
        output_buffer[branch.first] = 0.0;
        output_tree->Branch(branch.first.c_str(), &(output_buffer[branch.first]), (branch.first + "/F").c_str());
    }
    // Different readout for embedded samples.
    // Set up read out of quantities from event.
    for (auto &br : _float_inputs)
    {
        inputtree->SetBranchAddress(br.first.c_str(), &(br.second));
    }
    for (auto &br : _int_inputs)
    {
        inputtree->SetBranchAddress(br.first.c_str(), &(br.second));
    }

    // Start event loop here.
    for (unsigned int i = _first_entry; i <= _last_entry; i++)
    {
        if ((i - _first_entry) % 1000 == 0)
            std::cout << " processing event " << i << ", last event: " << _last_entry << " [" << 100 * (float(i - _first_entry) / (_last_entry - _first_entry)) << " %]" << std::endl;
        // Get entry
        inputtree->GetEntry(i);
        for (auto branch : branches)
        {
            output_buffer[branch] = get_quantity(fns_[branch], arguments, var_map);
            output_buffer["met_var_qcd"] = calculate_met_var_qcd(_float_inputs);
            output_buffer["met_var_w"] = calculate_met_var_w(_float_inputs);
            output_buffer["os"] = calculate_os(_float_inputs);
        }
        output_tree->Fill();
    }
    output_tree->Write();
    std::cout << "Tree successfully written." << std::endl;
    out->Close();
}

ImperialFakeFactorsProducer::~ImperialFakeFactorsProducer()
{
}

Float_t ImperialFakeFactorsProducer::calculate_met_var_qcd(std::map<std::string, Float_t> &outbuf)
{
    auto hadron_vec = ROOT::Math::PtEtaPhiMVector(outbuf["pt_2"], outbuf["eta_2"], outbuf["phi_2"], outbuf["mt_2_puppi"]);
    auto met_vec = ROOT::Math::PtEtaPhiEVector(outbuf["puppimet"], 0, outbuf["puppimetphi"], outbuf["puppimet"]);
    auto delta_phi = ROOT::Math::VectorUtil::DeltaPhi(met_vec, hadron_vec);
    auto met_pt = sqrt(met_vec.px() * met_vec.px() + met_vec.py() * met_vec.py());
    return met_pt / outbuf["pt_2"] * cos(delta_phi);
}

Float_t ImperialFakeFactorsProducer::calculate_met_var_w(std::map<std::string, Float_t> &outbuf)
{
    // construct met and lepton, tau 4 vectors
    auto lepton_vec = ROOT::Math::PtEtaPhiMVector(outbuf["pt_1"], outbuf["eta_1"], outbuf["phi_1"], outbuf["mt_1_puppi"]);
    auto hadron_vec = ROOT::Math::PtEtaPhiMVector(outbuf["pt_2"], outbuf["eta_2"], outbuf["phi_2"], outbuf["mt_2_puppi"]);
    auto met_vec = ROOT::Math::PtEtaPhiEVector(outbuf["puppimet"], 0, outbuf["puppimetphi"], outbuf["puppimet"]);
    // make the vectorial sum of the two, and use the resulting vector instead
    // the pure met vector
    auto comb = lepton_vec + met_vec;
    auto met_pt = sqrt(comb.px() * comb.px() + comb.py() * comb.py());
    auto delta_phi = ROOT::Math::VectorUtil::DeltaPhi(comb, hadron_vec);
    return met_pt / outbuf["pt_2"] * cos(delta_phi);
}

Float_t ImperialFakeFactorsProducer::calculate_os(std::map<std::string, Float_t> &outbuf)
{
    // same sign = 0 / opposite sign = 1
    double sign = 0.0;
    if (outbuf["q_1"] != outbuf["q_2"])
        sign = 1.0;
    return sign;
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

    std::string configpath = cmsswbase +
                             "/src/HiggsAnalysis/friend-tree-producer/data/config_imperialFakeFactorProducer.json";

    std::string workspacepath = cmsswbase + "/src/HiggsAnalysis/friend-tree-producer/data/imperial_ff/fakefactors_ws_" + channel + "_mssm_" +  std::to_string(year) + "_v2.root";
    std::cout << "Using " << workspacepath << std::endl;
    auto outputname =
        outputname_from_settings(input, folder, first_entry, last_entry);
    boost::filesystem::create_directories(filename_from_inputpath(input));

    ImperialFakeFactorsProducer producer = ImperialFakeFactorsProducer(input,
                                                                       year,
                                                                       first_entry,
                                                                       last_entry,
                                                                       workspacepath,
                                                                       outputname,
                                                                       configpath,
                                                                       tree,
                                                                       channel,
                                                                       folder);
    producer.run();

    return 0;
}