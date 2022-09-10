import argparse
import copy
import gc
import json
import logging
import os
import queue
import shutil
import subprocess
from multiprocessing import Lock, Process, Queue, current_process

import ROOT as r

logger = logging.getLogger("job_managment")
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")


def parse_arguments():
    # args are executable cores custom_workdir_path mode
    parser = argparse.ArgumentParser(
        description="Run the executable on the remote machine"
    )
    parser.add_argument("--executable", type=str, help="The executable to run")
    parser.add_argument(
        "--custom_workdir_path", type=str, help="The path to the remote machine"
    )
    parser.add_argument("--mode", type=str, help="The mode to run the executable")
    parser.add_argument("--cores", type=int, help="The number of cores")
    return parser.parse_args()


def write_trees_to_files(tasks_to_accomplish, tasks_that_are_done):
    while True:
        try:
            info = tasks_to_accomplish.get_nowait()
        except queue.Empty:
            break
        nick = info[0]
        collection_path = info[1]
        db = info[2]
        era = info[3]
        channel = info[4]
        # from the nick, remove the last _#number
        foldername = "_".join(nick.split("_")[:-1])
        # the outputpath follows the pattern era/channel/nick/files.root
        nick_path = os.path.join(collection_path, era, foldername, channel)
        if not os.path.exists(nick_path):
            os.makedirs(nick_path)
        outputfilepath = os.path.join(nick_path, nick + ".root")
        ##################
        # dump one db to a json file
        print(
            "For {} combining {} files".format(
                nick, len([db[friend]["files"] for friend in db.keys()])
            )
        )
        with open("jsondump.json", "w") as configdump:
            json.dump(db, configdump)
        chain = r.TChain("ntuple")
        friendchains = {}
        leavelist = []
        # prepare the friend chains
        for friend in db.keys():
            if friend != "nominal":
                friendchains[friend] = r.TChain("ntuple")

        for friend in db.keys():
            if friend == "nominal":
                for rootfile_path in db[friend]["files"]:
                    chain.Add(rootfile_path)
                    # rfile = r.TFile.Open(rootfile_path)
                    # columns.extend([b.GetName() for b in rfile.Get("ntuple").GetListOfLeaves()])
            else:
                for rootfile_path in db[friend]["files"]:
                    # if the file is not available, skip it
                    if os.path.exists(rootfile_path):
                        rfile = r.TFile.Open(rootfile_path)
                        leavelist.append(
                            copy.deepcopy(rfile.Get("ntuple").GetListOfLeaves())
                        )
                        rfile.Close()
                        r.SetOwnership(rfile, True)
                        if len(leavelist[-1]) != 0:
                            friendchains[friend].Add(rootfile_path)
                        else:
                            pass
        # add the friend chains to the base tree
        for friend in friendchains.keys():
            if friendchains[friend].GetEntries() != 0:
                chain.AddFriend(friendchains[friend], friend)
        # now create an TDataFrame from the chain
        df = r.RDataFrame(chain)
        chain_numentries = df.Count().GetValue()
        if chain_numentries == 0:
            logger.fatal("Chain does not contain any events.")
            raise Exception
        opt = r.RDF.RSnapshotOptions()
        opt.fMode = "RECREATE"
        df.Snapshot("ntuple", outputfilepath, ".*", opt)
        # now delete all the stuff, thanks root ....
        r.SetOwnership(chain, True)
        r.SetOwnership(df, True)
        chain.Reset()
        for friend in friendchains.keys():
            r.SetOwnership(friendchains[friend], True)
            friendchains[friend].Reset()
        del chain
        del df
        del leavelist
        gc.collect()
        print(
            "Found {} events for {}, written to outputfile".format(
                chain_numentries, nick
            )
        )
        tasks_that_are_done.put(nick)
    return True


def collect_outputs(executable, cores, custom_workdir_path, mode):
    if custom_workdir_path:
        workdir_path = os.path.join(
            custom_workdir_path, executable.replace(".py", "") + "_workdir"
        )
    else:
        workdir_path = os.path.join(
            os.environ["CMSSW_BASE"], "src", executable.replace(".py", "") + "_workdir"
        )
    jobdb_path = os.path.join(
        workdir_path, "condor_" + executable.replace(".py", "") + ".json"
    )
    gc_path = os.path.join(
        workdir_path, "grid_control_{}.conf".format(executable.replace(".py", ""))
    )
    jobdb_file = open(jobdb_path, "r")
    jobdb = json.loads(jobdb_file.read())
    datasetdb = {}
    collection_path = os.path.join(
        workdir_path, executable.replace(".py", "") + "_collected"
    )
    if not os.path.exists(collection_path):
        os.makedirs(collection_path)
    for jobnumber in sorted([int(k) for k in jobdb.keys()]):
        for subjobnumber in range(len(jobdb[str(jobnumber)])):
            nick = str(
                jobdb[str(jobnumber)][subjobnumber]["input"]
                .split("/")[-1]
                .replace(".root", "")
            )
            channel = str(jobdb[str(jobnumber)][subjobnumber]["channel"])
            folder = str(jobdb[str(jobnumber)][subjobnumber]["folder"])
            era = str(jobdb[str(jobnumber)][subjobnumber]["era"])
            first = jobdb[str(jobnumber)][subjobnumber]["first_entry"]
            last = jobdb[str(jobnumber)][subjobnumber]["last_entry"]
            first_last_combination = (first, last)
            filename = "_".join([nick, folder, era, channel, str(first), str(last)]) + ".root"
            filepath = os.path.join(workdir_path, "{}_{}_{}".format(era, channel, nick), filename)
            if nick not in datasetdb:
                datasetdb[nick] = {}
            if era not in datasetdb[nick]:
                datasetdb[nick][era] = {}
            if channel not in datasetdb[nick][era]:
                datasetdb[nick][era][channel] = {}
            if folder not in datasetdb[nick][era][channel].keys():
                datasetdb[nick][era][channel][folder] = {}
                datasetdb[nick][era][channel][folder]["files"] = []
                datasetdb[nick][era][channel][folder]["first_last"] = []

            datasetdb[nick][era][channel][folder]["first_last"].append(
                first_last_combination
            )
            datasetdb[nick][era][channel][folder]["files"].append(filepath)
            datasetdb[nick][era][channel][folder]["workdir"] = workdir_path

    if mode == "local":
        commands = []
        nicks = sorted(datasetdb)
        for nick in nicks:
            for era in datasetdb[nick]:
                for channel in datasetdb[nick][era]:
                    # create the folders for each nick
                    # from the nick, remove the last _#number
                    foldername = "_".join(nick.split("_")[:-1])
                    nick_path = os.path.join(collection_path, era, foldername, channel)
                    if not os.path.exists(nick_path):
                        os.makedirs(nick_path)
                    commands.append(
                        [
                            nick,
                            collection_path,
                            datasetdb[nick][era][channel],
                            era,
                            channel,
                        ]
                    )
        # pool = Pool(cores)
        # pool.map(write_trees_to_files, commands)
        # r.EnableImplicitMT(cores)
        # for command in commands:
        #     write_trees_to_files(command)
    else:
        raise NotImplementedError("Mode {} is not implemented".format(mode))

    number_of_task = len(commands)
    print("will run {} tasks using {} processes".format(number_of_task, cores))
    number_of_processes = cores
    tasks_to_accomplish = Queue()
    tasks_that_are_done = Queue()
    processes = []

    for command in commands:
        tasks_to_accomplish.put(command)

    # creating processes
    for w in range(number_of_processes):
        p = Process(
            target=write_trees_to_files, args=(tasks_to_accomplish, tasks_that_are_done)
        )
        processes.append(p)
        p.start()

    # completing process
    for p in processes:
        p.join()

    # print the output
    while not tasks_that_are_done.empty():
        print(tasks_that_are_done.get())


# main function
if __name__ == "__main__":
    args = parse_arguments()
    collect_outputs(args.executable, args.cores, args.custom_workdir_path, args.mode)
