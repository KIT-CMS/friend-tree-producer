#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Example:
        python HiggsAnalysis/friend-tree-producer/scripts/singletauTriggerWeightsProducer.py \
            --input /ceph/htautau/deeptau_04-27/2018/ntuples/Embedding2018A_ElTauFinalState_inputDoubleMu102XminiAODv1_13TeV_USER_v1/Embedding2018A_ElTauFinalState_inputDoubleMu102XminiAODv1_13TeV_USER_v1.root \
            --output-dir . --rooworkspace-file /work/mburkart/Run2Legacy/LegacyCorrectionsWorkspace/output/htt_scalefactors_legacy_2018.root \
            --start 0 --end 11 --pipeline tt_nominal \
"""
import re
import os
import sys
import json
import yaml
import ROOT
import numpy
import copy
from array import array
import six
import argparse
import logging
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
        description="Write single tau trigger weights to friend-trees.")
    parser.add_argument(
            "--input",
           required=True,
           type=str,
           help="Input file."
    )
    parser.add_argument(
            "--tree",
            default="ntuple",
            type=str,
            help="Name of the root tree."
    )
    parser.add_argument(
            "--enable-logging",
            action="store_true",
            help="Enable loggging for debug purposes.")
    parser.add_argument(
            "--cmsswbase",
            default=os.environ['CMSSW_BASE'],
            help="Set path for to local cmssw for submission with Grid-Control"
    )
    parser.add_argument(
            "--first-entry",
            "--first_entry",
            "--start",
            default=0,
            type=int,
            help="Index of first event to process."
    )
    parser.add_argument(
            "--last-entry",
            "--last_entry",
            "--end",
            default=-1,
            type=int,
            help="Index of last event to process."
    )
    parser.add_argument(
            "--pipeline",
            "--pipelines",
            "--folder",
            nargs="?",
            default=None,
            type=str,
            help="Directory within rootfile."
    )
    parser.add_argument(
            "--output-dir",
            type=str,
            default='.',
            help="Tag of output files."
    )
    parser.add_argument(
            "--rooworkspace-file",
            "--rooworkspace_file",
            type=str,
            default="workspace.root",
            help="Path to the file containing the workspace."
    )
    parser.add_argument(
            "--input-friends",
            "--input_friends",
            type=str,
            default=None,
            help="Path to the friend-tree used."
    )
    return parser.parse_args()


def calculate_emb_med_region(event, outbuf, wp="tight"):
    wp = wp.lower()
    # eff_ditau_1 = getattr(event, "crossTriggerEMBEfficiencyWeight_{}_DeepTau_1".format(wp))
    # eff_ditau_2 = getattr(event, "crossTriggerEMBEfficiencyWeight_{}_DeepTau_2".format(wp))
    eff_ditau_1 = outbuf["t_trg_double_t_emb_1"][0]
    eff_ditau_2 = outbuf["t_trg_double_t_emb_2"][0]
    eff_stau_1 = outbuf["t_trg_single_t_emb_1"][0]
    result = (eff_ditau_1 * eff_ditau_2) + eff_stau_1 - (eff_ditau_2 * eff_stau_1)
    return result


def calculate_emb_high_region(event, outbuf, wp="tight"):
    wp = wp.lower()
    # eff_ditau_1 = getattr(event, "crossTriggerEMBEfficiencyWeight_{}_DeepTau_1".format(wp))
    # eff_ditau_2 = getattr(event, "crossTriggerEMBEfficiencyWeight_{}_DeepTau_2".format(wp))
    eff_ditau_1 = outbuf["t_trg_double_t_emb_1"][0]
    eff_ditau_2 = outbuf["t_trg_double_t_emb_2"][0]
    eff_stau_1 = outbuf["t_trg_single_t_emb_1"][0]
    eff_stau_2 = outbuf["t_trg_single_t_emb_2"][0]
    result = eff_ditau_1*eff_ditau_2 + eff_stau_1 + eff_stau_2 \
             - eff_ditau_1*eff_stau_2 - eff_ditau_2*eff_stau_1
    return result


def calculate_mc_med_region(event, outbuf, wp="tight"):
    wp = wp.lower()
    # eff_ditau_1 = getattr(event, "crossTriggerMCEfficiencyWeight_{}_DeepTau_1".format(wp))
    # eff_ditau_2 = getattr(event, "crossTriggerMCEfficiencyWeight_{}_DeepTau_2".format(wp))
    eff_ditau_1 = outbuf["t_trg_double_t_mc_1"][0]
    eff_ditau_2 = outbuf["t_trg_double_t_mc_2"][0]
    eff_stau_1 = outbuf["t_trg_single_t_mc_1"][0]
    eff_stau_and_ditau_1 = outbuf["t_trg_singleANDdouble_t_mc_1"][0]
    result = eff_ditau_1*eff_ditau_2 + eff_stau_1 - eff_ditau_2*eff_stau_and_ditau_1
    return result


def calculate_mc_high_region(event, outbuf, wp="tight"):
    wp = wp.lower()
    # eff_ditau_1 = getattr(event, "crossTriggerMCEfficiencyWeight_{}_DeepTau_1".format(wp))
    # eff_ditau_2 = getattr(event, "crossTriggerMCEfficiencyWeight_{}_DeepTau_2".format(wp))
    eff_ditau_1 = outbuf["t_trg_double_t_mc_1"][0]
    eff_ditau_2 = outbuf["t_trg_double_t_mc_2"][0]
    eff_stau_1 = outbuf["t_trg_single_t_mc_1"][0]
    eff_stau_2 = outbuf["t_trg_single_t_mc_2"][0]
    eff_stau_and_ditau_1 = outbuf["t_trg_singleANDdouble_t_mc_1"][0]
    eff_stau_and_ditau_2 = outbuf["t_trg_singleANDdouble_t_mc_2"][0]
    result = eff_ditau_1*eff_ditau_2 + eff_stau_1 + eff_stau_2 \
             - eff_ditau_1*eff_stau_and_ditau_2 - eff_ditau_2*eff_stau_and_ditau_1 - eff_stau_1*eff_stau_2 \
             + eff_stau_and_ditau_1*eff_stau_and_ditau_2
    return result


def calculate_data_med_region(event, outbuf, wp="tight", is_emb=False):
    wp = wp.lower()
    eff_ditau_1 = getattr(event, "crossTriggerDataEfficiencyWeight_{}_DeepTau_1".format(wp))
    eff_ditau_2 = getattr(event, "crossTriggerDataEfficiencyWeight_{}_DeepTau_2".format(wp))
    if is_emb:
        sf_ditau_1 = eff_ditau_1/getattr(event, "crossTriggerCorrectedMCEfficiencyWeight_{}_DeepTau_1".format(wp))
        sf_ditau_2 = eff_ditau_2/getattr(event, "crossTriggerCorrectedMCEfficiencyWeight_{}_DeepTau_2".format(wp))
    else:
        sf_ditau_1 = eff_ditau_1/getattr(event, "crossTriggerMCEfficiencyWeight_{}_DeepTau_1".format(wp))
        sf_ditau_2 = eff_ditau_2/getattr(event, "crossTriggerMCEfficiencyWeight_{}_DeepTau_2".format(wp))
    sf_ditau = sf_ditau_1*sf_ditau_2
    eff_ditau_1_mc = outbuf["t_trg_double_t_mc_1"][0]
    eff_ditau_2_mc = outbuf["t_trg_double_t_mc_2"][0]
    eff_stau_1 = outbuf["t_trg_single_t_wstar_sf_1"][0]*outbuf["t_trg_single_t_mc_1"][0]
    eff_stau_and_ditau_1 = outbuf["t_trg_single_t_wstar_sf_1"][0]*outbuf["t_trg_singleANDdouble_t_mc_1"][0]
    result = sf_ditau*eff_ditau_1_mc*eff_ditau_2_mc + eff_stau_1 - sf_ditau*eff_ditau_2_mc*eff_stau_and_ditau_1
    return result


def calculate_data_high_region(event, outbuf, wp="tight", is_emb=False):
    wp = wp.lower()
    eff_ditau_1 = getattr(event, "crossTriggerDataEfficiencyWeight_{}_DeepTau_1".format(wp))
    eff_ditau_2 = getattr(event, "crossTriggerDataEfficiencyWeight_{}_DeepTau_2".format(wp))
    if is_emb:
        sf_ditau_1 = eff_ditau_1/getattr(event, "crossTriggerCorrectedMCEfficiencyWeight_{}_DeepTau_1".format(wp))
        sf_ditau_2 = eff_ditau_2/getattr(event, "crossTriggerCorrectedMCEfficiencyWeight_{}_DeepTau_2".format(wp))
    else:
        sf_ditau_1 = eff_ditau_1/getattr(event, "crossTriggerMCEfficiencyWeight_{}_DeepTau_1".format(wp))
        sf_ditau_2 = eff_ditau_2/getattr(event, "crossTriggerMCEfficiencyWeight_{}_DeepTau_2".format(wp))
    sf_ditau = sf_ditau_1*sf_ditau_2
    eff_ditau_1_mc = outbuf["t_trg_double_t_mc_1"][0]
    eff_ditau_2_mc = outbuf["t_trg_double_t_mc_2"][0]
    eff_stau_1 = outbuf["t_trg_single_t_wstar_sf_1"][0]*outbuf["t_trg_single_t_mc_1"][0]
    eff_stau_2 = outbuf["t_trg_single_t_wstar_sf_2"][0]*outbuf["t_trg_single_t_mc_2"][0]
    eff_stau_and_ditau_1 = outbuf["t_trg_single_t_wstar_sf_1"][0]*outbuf["t_trg_singleANDdouble_t_mc_1"][0]
    eff_stau_and_ditau_2 = outbuf["t_trg_single_t_wstar_sf_2"][0]*outbuf["t_trg_singleANDdouble_t_mc_2"][0]
    result = sf_ditau*eff_ditau_1_mc*eff_ditau_2_mc + eff_stau_1 + eff_stau_2 \
             - sf_ditau*eff_ditau_1_mc*eff_stau_and_ditau_2 - sf_ditau*eff_ditau_2_mc*eff_stau_and_ditau_1 - eff_stau_1*eff_stau_2 \
             + sf_ditau*eff_ditau_1_mc*eff_ditau_2_mc*eff_stau_and_ditau_1*eff_stau_and_ditau_2
    return result


class SingletauWeightsProducer(object):

    additional_variables_emb = [
        "singleTauTriggerEMBWeightHighPt",
        "singleTauTriggerEMBWeightMedPt",
        "singleTauTriggerDataWeightHighPt",
        "singleTauTriggerDataWeightMedPt",
        "singleTauTriggerSFWeightHighPt",
        "singleTauTriggerSFWeightMedPt"
        ]
    additional_variables = [
        "singleTauTriggerMCWeightHighPt",
        "singleTauTriggerMCWeightMedPt",
        "singleTauTriggerDataWeightHighPt",
        "singleTauTriggerDataWeightMedPt",
        "singleTauTriggerSFWeightHighPt",
        "singleTauTriggerSFWeightMedPt",
        ]


    def __init__(self, era, inputfile, outputfile, eventrange, workspace,
                 config, treename, channel, pipelines, friend_file):
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
        self.friend_file = friend_file
        return

    def make_outputfile(self, outputfile):
        outputfile = os.path.abspath(outputfile)
        if not os.path.exists(os.path.dirname(outputfile)):
            os.makedirs(os.path.dirname(outputfile))
        return ROOT.TFile(outputfile, "recreate")

    def get_quantity(self, function, arguments, event, variable_mapping):
        # get quantities from the event.
        argset = self.workspace.argSet(",".join(variable_mapping[arg] for arg in arguments))
        roofunction = self.workspace.function(function)
        for par in arguments:
            value = getattr(event, par)
            logger.debug("Parameter: {} - Value: {}".format(par, value))
            argset.setRealValue(variable_mapping[par], value)
        result = roofunction.getVal(argset)
        logger.debug("Result: {} - Value: {}".format(function, result))
        return result

    def run(self):
        nickname = os.path.basename(self.inputfile.GetName()).replace(".root", "")
        for pipeline in self.pipelines:
            pipeline = pipeline.replace(self.channel+ "_", "")

            # Prepare inputs from ntuple
            input_tree = self.inputfile.Get(
                    "{}_{}/{}".format(self.channel, pipeline, self.treename))
            if self.friend_file is not None:
                input_tree.AddFriend("{}_{}/{}".format(self.channel, pipeline, self.treename),
                                     self.friend_file)

            # Only load branches necessary for the writeout of weights.
            input_tree.SetBranchStatus("*", 0)
            for var in self.variable_mapping.keys():
                input_tree.SetBranchStatus(var, 1)
            input_tree.SetBranchStatus("crossTrigger*", 1)

            output_root_dir = self.outputfile.mkdir("%s_%s" %
                                                    (self.channel, pipeline))
            output_root_dir.cd()
	    output_tree = ROOT.TTree(self.treename, self.treename)
            # Prepare branches
            output_buffer = {}
            branches = self.config.keys()
            add_variables = self.additional_variables_emb if "Embedding" in nickname \
                            else self.additional_variables
            for branch in branches + add_variables:
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
                    output_buffer[branch][0] = self.get_quantity(
                        self.config[branch]["WorkspaceWeightNames"],
                        self.config[branch]["WorkspaceObjectArguments"], event, self.variable_mapping)

                wp = "medium"
                if "Embedding" in nickname:
                    output_buffer["singleTauTriggerEMBWeightHighPt"][0] = calculate_emb_high_region(event, output_buffer, wp=wp)
                    output_buffer["singleTauTriggerEMBWeightMedPt"][0] = calculate_emb_high_region(event, output_buffer, wp=wp)
                    output_buffer["singleTauTriggerDataWeightHighPt"][0] = calculate_data_high_region(event, output_buffer, wp=wp, is_emb=True)
                    output_buffer["singleTauTriggerDataWeightMedPt"][0] = calculate_data_med_region(event, output_buffer, wp=wp, is_emb=True)
                    output_buffer["singleTauTriggerSFWeightHighPt"][0] = output_buffer["singleTauTriggerDataWeightHighPt"][0]/output_buffer["singleTauTriggerEMBWeightHighPt"][0]
                    output_buffer["singleTauTriggerSFWeightMedPt"][0] = output_buffer["singleTauTriggerDataWeightMedPt"][0]/output_buffer["singleTauTriggerEMBWeightMedPt"][0]
                else:
                    output_buffer["singleTauTriggerMCWeightHighPt"][0] = calculate_mc_high_region(event, output_buffer, wp=wp)
                    output_buffer["singleTauTriggerMCWeightMedPt"][0] = calculate_mc_med_region(event, output_buffer, wp=wp)
                    output_buffer["singleTauTriggerDataWeightHighPt"][0] = calculate_data_high_region(event, output_buffer, wp=wp)
                    output_buffer["singleTauTriggerDataWeightMedPt"][0] = calculate_data_med_region(event, output_buffer, wp=wp)
                    output_buffer["singleTauTriggerSFWeightHighPt"][0] = output_buffer["singleTauTriggerDataWeightHighPt"][0]/output_buffer["singleTauTriggerMCWeightHighPt"][0]
                    output_buffer["singleTauTriggerSFWeightMedPt"][0] = output_buffer["singleTauTriggerDataWeightMedPt"][0]/output_buffer["singleTauTriggerMCWeightMedPt"][0]
                output_tree.Fill()

            # Save
            output_tree.Write()
            print("Tree successfully written")


def main(args):
    nickname = os.path.basename(args.input).replace(".root", "")
    channel = args.pipeline.split("_")[0]
    cmsswbase = args.cmsswbase

    # Determine era
    with open(cmsswbase +
              "/src/HiggsAnalysis/friend-tree-producer/data/input_params/datasets.json"
              ) as json_file:
        datasets = json.load(json_file)
    era = str(datasets[nickname]["year"])

    # Set workspace path
    # workspacepath = cmsswbase + \
    #         "/src/HiggsAnalysis/friend-tree-producer/data/{}.root"
    workspacepath = args.rooworkspace_file
    configpath = cmsswbase + \
            "/src/HiggsAnalysis/friend-tree-producer/data/config_singletauTriggerWeightsProducer.yaml"

    # Create friend tree
    output_path = os.path.join(args.output_dir, nickname)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    outputfile = os.path.join(
            output_path,
            "_".join(
                filter(None,
                       [nickname,
                        args.pipeline,
                        str(args.first_entry),
                        str(args.last_entry),
                        ],
                       )) + ".root",
                )
    producer = SingletauWeightsProducer(
            era=era,
            inputfile=args.input,
            outputfile=outputfile,
            eventrange=[args.first_entry, args.last_entry],
            workspace=workspacepath,
            config=configpath,
            treename=args.tree,
            channel=channel,
            pipelines="all" if args.pipeline is None else [args.pipeline],
            friend_file=args.input_friends)
    producer.run()
    return


if __name__ == "__main__":
    args = parse_arguments()

    if args.enable_logging:
        setup_logging(
                "SingletauWeightsProducer%s_%s_%s_%s.log" % (
                    os.path.basename(args.input).replace(".root", ""),
                    args.pipeline,
                    args.first_entry,
                    args.last_entry,
                    ),
                logging.DEBUG,
                )

    main(args)
