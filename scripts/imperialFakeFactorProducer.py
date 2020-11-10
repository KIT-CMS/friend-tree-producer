#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Example:
        python scripts/imperialFakeFactorProducer.py \
            --input /ceph/htautau/deeptau_04-27/2018/ntuples/Embedding2018A_ElTauFinalState_inputDoubleMu102XminiAODv1_13TeV_USER_v1/Embedding2018A_ElTauFinalState_inputDoubleMu102XminiAODv1_13TeV_USER_v1.root \
            --output-dir . --rooworkspace-file /ceph/sbrommer/copy_dir/fakefactors_ws_et_mssm_2017.root \
            --start 0 --end 1000 --pipeline et_nominal --config data/config_imperialFakeFactorProducer.yaml
"""
import os
import json
import yaml
import ROOT
from array import array
import argparse
import logging
import numpy as np
logger = logging.getLogger()


def setup_logging(output_file, level=logging.DEBUG):
    logger.setLevel(level)
    formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    file_handler = logging.FileHandler(output_file, "w")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="create friend trees for imperial fake factors form the given RooWorkspace"
    )
    parser.add_argument("--input", required=True, type=str, help="Input file.")

    parser.add_argument("--tree",
                        default="ntuple",
                        type=str,
                        help="Name of the root tree.")
    parser.add_argument(
        "--enable-logging",
        action="store_true",
        help="Enable loggging for debug purposes.",
    )
    parser.add_argument(
        "--cmsswbase",
        default=os.environ["CMSSW_BASE"],
        help="Set path for to local cmssw for submission with Grid-Control",
    )

    parser.add_argument(
        "--first-entry",
        "--first_entry",
        "--start",
        default=0,
        type=int,
        help="Index of first event to process.",
    )
    parser.add_argument(
        "--last-entry",
        "--last_entry",
        "--end",
        default=-1,
        type=int,
        help="Index of last event to process.",
    )

    parser.add_argument(
        "--pipeline",
        "--pipelines",
        "--folder",
        nargs="?",
        default=None,
        type=str,
        help="Directory within rootfile.",
    )
    parser.add_argument("--output-dir",
                        type=str,
                        default=".",
                        help="Tag of output files.")

    parser.add_argument("--dry",
                        action="store_true",
                        default=False,
                        help="dry run")

    return parser.parse_args()


def calc_delta_phi(phi1, phi2):
    dPhi = phi1 - phi2
    if(dPhi > np.pi):
        return dPhi - 2 * np.pi
    if(dPhi < -np.pi):
        return dPhi + 2 * np.pi
    return dPhi


def calculate_met_var_qcd(event):
    met_phi = event.puppimetphi
    met_et = event.puppimetsumet
    hadron_pt = event.pt_2
    hadron_phi = event.phi_2
    delta_phi = calc_delta_phi(met_phi, hadron_phi)
    return met_et / hadron_pt * np.cos(delta_phi)


def calculate_met_var_w(event):
    lepton_pt = event.pt_1
    lepton_eta = event.eta_1
    lepton_phi = event.phi_1
    lepton_mt = event.mt_1_puppi
    met = event.puppimet
    met_phi = event.puppimetphi
    met_eta = 0.0
    met_et = event.puppimetsumet
    hadron_pt = event.pt_2
    hadron_phi = event.phi_2
    # construct met and lepton 4 vectors
    lepton_vec = ROOT.Math.PtEtaPhiMVector(
        lepton_pt, lepton_eta, lepton_phi, lepton_mt)
    met_vec = ROOT.Math.PtEtaPhiEVector(met, met_eta, met_phi, met_et)
    # make the vectorial sum of the two, and use the resulting vector instead
    # the pure met vector
    comb = lepton_vec + met_vec
    delta_phi = calc_delta_phi(comb.phi(), hadron_phi)
    return comb.Et() / hadron_pt * np.cos(delta_phi)


class FakeFactorProducer(object):

    def __init__(self, era, inputfile, outputfile, eventrange, workspace,
                 config, treename, channel, pipelines):
        self.inputfile = ROOT.TFile(inputfile, "read")
        self.era = era
        self.eventrange = eventrange
        workspacefile = ROOT.TFile(workspace, "read")
        self.workspace = workspacefile.Get("w")
        self.outputfile = self.make_outputfile(outputfile)
        self.config = yaml.load(open(config))["rooworkspace"]
        self.variable_mapping = yaml.load(open(config))["map_arguments"]
        self.treename = treename
        self.channel = channel
        self.pipelines = pipelines

        if self.pipelines == "all":
            self.pipelines = []
            for key in self.inputfile.GetListOfKeys():
                if key.GetName().startswith(self.channel[0]):
                    pipelines.append(
                        key.GetName().split(
                            "/")[0].replace(self.channel[0] + "_", "")
                    )
            print("process pipelines: %s" % pipelines)

    def make_outputfile(self, outputfile):
        print outputfile
        outputfile = os.path.abspath(outputfile)
        if not os.path.exists(os.path.dirname(outputfile)):
            os.makedirs(os.path.dirname(outputfile))
        return ROOT.TFile(outputfile, "recreate")

    def getQuantity(self, function, arguments, event, variable_mapping):
        # get quantities from the event
        argset = self.workspace.argSet(",".join(arguments))
        roofunction = self.workspace.function(function)
        for para in arguments:
            if para == 'met_var_qcd':
                value = calculate_met_var_qcd(event)
            elif para == 'met_var_w':
                value = calculate_met_var_w(event)
            elif para == 'os':
                # same sign = 0 / opposite sign = 1
                if event.q_1 == event.q_2:
                    value = 0.0
                else:
                    value = 1.0
            else:
                value = getattr(event, variable_mapping[para])
            logger.debug("Parameter: {} - Value: {}".format(para, value))
            argset.setRealValue(para, value)
        result = roofunction.getVal(argset)
        logger.debug("Result: {} - Value: {}".format(function, result))
        return result

    def run(self):
        self.functions = {}
        for pipeline in self.pipelines:
            pipeline = pipeline.replace(self.channel + "_", "")
            # Prepare data inputs
            input_tree = self.inputfile.Get(
                "%s_%s/%s" % (self.channel, pipeline, self.treename))
            output_root_dir = self.outputfile.mkdir("%s_%s" %
                                                    (self.channel, pipeline))
            output_root_dir.cd()
            output_tree = ROOT.TTree(self.treename, self.treename)
            # Prepare branches
            output_buffer = {}
            # add up and down variation for all uncertainties
            branches = self.config[self.channel]["functions"]["main"].keys()
            for uncertainty in self.config[self.channel]["functions"]["uncertainties"]:
                branches.extend(["{}_up".format(uncertainty),
                                 "{}_down".format(uncertainty)])
            for branch in branches:
                output_buffer[branch] = array("d", [0])
                output_tree.Branch(branch, output_buffer[branch],
                                   "%s/D" % branch)
                output_tree.SetBranchAddress(branch, output_buffer[branch])
            # Fill tree
            if self.eventrange[1] > 0:
                nev = self.eventrange[1] - self.eventrange[0] + 1
                if self.eventrange[1] >= input_tree.GetEntries():
                    raise Exception("The last entry exceeds maximum")
            else:
                nev = input_tree.GetEntries() - self.eventrange[0]
            printout_on = [
                self.eventrange[0] + i * int(nev / 10)
                for i in range(0, 11)
            ]
            for evt_i, event in enumerate(input_tree):
                if self.eventrange is not None:
                    if evt_i < self.eventrange[0]:
                        continue
                    elif (
                            evt_i > self.eventrange[1]
                            and self.eventrange[1] >= 0
                    ):  # latter condition allows to set negative upper limit in order to have it ignored
                        break
                    elif evt_i in printout_on:
                        print("\t ...processing %d %% [Event %d]" %
                              (printout_on.index(evt_i) * 10, evt_i))

                # Evaluating weights
                for branch in branches:
                    output_buffer[branch] = self.getQuantity(
                        branch,
                        self.config[self.channel]["arguments"], event, self.variable_mapping[self.channel])

                output_tree.Fill()
            # Save
            output_tree.Write()
            print("Tree successfully written")


def main(args):
    nickname = os.path.basename(args.input).replace(".root", "")
    channel = args.pipeline.split("_")[0]
    # Get path to cmssw and fakefactordatabase
    cmsswbase = args.cmsswbase

    # Determine era
    with open(
            cmsswbase +
            "/src/HiggsAnalysis/friend-tree-producer/data/input_params/datasets.json"
    ) as json_file:
        datasets = json.load(json_file)
    era = str(datasets[nickname]["year"])

    # Set Workspace Path
    workspacepath = cmsswbase + \
        "/src/HiggsAnalysis/friend-tree-producer/data/imperial_ff/fakefactors_ws_{}_mssm_{}.root".format(
            channel, era)
    # Set config path
    configpath = cmsswbase + \
        "/src/HiggsAnalysis/friend-tree-producer/data/config_imperialFakeFactorProducer.yaml"

    # Create friend tree
    output_path = os.path.join(args.output_dir, nickname)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    outputfile = os.path.join(
        output_path,
        "_".join(
            filter(
                None,
                [
                    nickname,
                    args.pipeline,
                    str(args.first_entry),
                    str(args.last_entry),
                ],
            )) + ".root",
    )
    producer = FakeFactorProducer(
        era=era,
        inputfile=args.input,
        outputfile=outputfile,
        eventrange=[args.first_entry, args.last_entry],
        workspace=workspacepath,
        config=configpath,
        treename=args.tree,
        channel=channel,
        pipelines="all" if args.pipeline is None else [args.pipeline])
    producer.run()


if __name__ == "__main__":
    args = parse_arguments()

    if args.enable_logging:
        setup_logging(
            "imperialFakeFactorProducer%s_%s_%s_%s.log" % (
                os.path.basename(args.input).replace(".root", ""),
                args.folder,
                args.first_entry,
                args.last_entry,
            ),
            logging.WARNING,
        )

    main(args)
