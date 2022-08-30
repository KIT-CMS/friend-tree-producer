import os, json, logging
import argparse
logger = logging.getLogger("job_managment")
from streampaths import *
import ROOT as r
import shutil
import subprocess
from multiprocessing import Pool

r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")


def parse_arguments():
    # args are executable cores custom_workdir_path mode
    parser = argparse.ArgumentParser(description='Run the executable on the remote machine')
    parser.add_argument('--executable', type=str, help='The executable to run')
    parser.add_argument('--custom_workdir_path', type=str, help='The path to the remote machine')
    parser.add_argument('--mode', type=str, help='The mode to run the executable')
    parser.add_argument('--cores', type=int, help='The number of cores')
    return parser.parse_args()



def write_trees_to_files(info):
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
    # hadd version does normally not work :-(
    # hadd_cmd = ["hadd -f", outputfile]
    # # unpack the filenames to total paths
    # for filepath in db["files"]:
    #     full_filepath = os.path.join(collection_path, filepath)
    #     hadd_cmd.append(full_filepath)
    # try:
    #     print(" ".join(hadd_cmd))
    #     # dump the hadd_cmd to a file
    #     with open("hadd_cmd.txt", "w") as hadd_file:
    #         hadd_file.write(" ".join(hadd_cmd))
    #     # run the hadd command
    #     exit()
    #     output = subprocess.check_output(" ".join(hadd_cmd), shell=True)
    #     # use decode function to convert to string
    #     # print('###############')
    #     # print('Output:', output.decode("utf-8"))
    # except subprocess.CalledProcessError as error:
    #     print('Error code:', error.returncode, '. Output:', error.output.decode("utf-8"))
    ##################
    # dump one db to a json file
    with open("jsondump.json", "w") as configdump:
        json.dump(db, configdump)
    columns = []
    chain = r.TChain("ntuple")
    friendchains = {}
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
                rfile = r.TFile.Open(rootfile_path)
                if len(rfile.Get("ntuple").GetListOfLeaves()) != 0:
                    friendchains[friend].Add(rootfile_path)
                else:
                    pass
                rfile.Close()
                # rfile = r.TFile.Open(rootfile_path)
                # columns.extend([b.GetName() for b in rfile.Get("ntuple").GetListOfLeaves()])
    # for first_last in db["first_last"]:
    #     print(first_last)
    #     nominal_file = "_".join([nick, "nominal", era, channel, str(first_last[0]), str(first_last[1])]) + ".root"
    #     nominal_file_path = os.path.join(db["workdir"], "{}_{}_{}".format(era, channel, nick), nominal_file)
    #     chain.AddFile(nominal_file_path)
    #     # open the file to gat all the branches
    #     rootfile = r.TFile.Open(nominal_file_path)
    #     columns.add(set([b.GetName() for b in rootfile.Get("ntuple").GetListOfLeaves()]))
    #     # now add all the friends to the friend chains
    #     for friend in db["folders"]:
    #         if friend != "nominal":
    #             friend_file = "_".join([nick, friend, era, channel, str(first_last[0]), str(first_last[1])]) + ".root"
    #             friend_file_path = os.path.join(db["workdir"], "{}_{}_{}".format(era, channel, nick), friend_file)
    #             friendchains[friend].AddFile(friend_file_path)
    #             rootfile = r.TFile.Open(friend_file_path)
    #             columns.add(set([b.GetName() for b in rootfile.Get("ntuple").GetListOfLeaves()]))
    # add the friend chains to the base tree
    for friend in friendchains.keys():
        chain.AddFriend(friendchains[friend], friend)
    # now create an TDataFrame from the chain
    df = r.RDataFrame(chain)
    chain_numentries = df.Count().GetValue()
    if chain_numentries == 0:
        logger.fatal("Chain does not contain any events.")
        raise Exception
    print("Found {} events for {}".format(chain_numentries, nick))
    opt = r.RDF.RSnapshotOptions()
    opt.fMode = "RECREATE"
    df.Snapshot("ntuple", outputfilepath, ".*", opt)

def collect_outputs(executable, cores, custom_workdir_path, mode):
    if custom_workdir_path:
        workdir_path = os.path.join(custom_workdir_path, executable.replace('.py', '') + "_workdir")
    else:
        workdir_path = os.path.join(
            os.environ["CMSSW_BASE"], "src", executable.replace('.py', '') + "_workdir"
        )
    jobdb_path = os.path.join(workdir_path, "condor_" + executable.replace('.py', '') + ".json")
    gc_path = os.path.join(workdir_path, "grid_control_{}.conf".format(executable.replace('.py', '')))
    jobdb_file = open(jobdb_path, "r")
    jobdb = json.loads(jobdb_file.read())
    datasetdb = {}
    collection_path = os.path.join(workdir_path, executable.replace('.py', '') + "_collected")
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
            # print("Adding nick: {}, channel: {}, folder: {}, era: {}, first: {}, last: {}".format(nick, channel, folder, era, first, last))
            if nick not in datasetdb:
                datasetdb[nick] = {}
            if era not in datasetdb[nick]:
                datasetdb[nick][era] = {}
            if channel not in datasetdb[nick][era]:
                datasetdb[nick][era][channel] = {}
            if folder not in datasetdb[nick][era][channel].keys():
                # print("{} not in {}".format(folder, datasetdb[nick][era][channel].keys()))
                datasetdb[nick][era][channel][folder] = {}
                datasetdb[nick][era][channel][folder]["files"] = []
                datasetdb[nick][era][channel][folder]["first_last"] = []
                # print("{} is now in {}".format(folder, datasetdb[nick][era][channel].keys()))

            # elif mode == "xrootd":
            #     with open(gc_path, "r") as gc_file:
            #         for line in gc_file.readlines():
            #             if "se path" in line:
            #                 filepath = (
            #                     server_xrootd["GridKA"]
            #                     + "/store/"
            #                     + line.split("/store/")[1].strip("\n")
            #                     + "/"
            #                     + filename
            #                 )
            #                 break

            datasetdb[nick][era][channel][folder]["first_last"].append(first_last_combination)
            datasetdb[nick][era][channel][folder]["files"].append(filepath)
            datasetdb[nick][era][channel][folder]["workdir"] = workdir_path
            # print("Added nick: {}, channel: {}, folder: {}, era: {}, first: {}, last: {}".format(nick, channel, folder, era, first, last))


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
                    # print("Writing {} to {}".format(nick, nick_path))
                    if not os.path.exists(nick_path):
                        os.makedirs(nick_path)
                    commands.append([nick, collection_path, datasetdb[nick][era][channel], era, channel, ])
                    # write_trees_to_files([nick, collection_path, datasetdb[nick][era][channel], era, channel])
        print("Using {} cores to collect outputs".format(cores))
        pool = Pool(cores)
        pool.map(write_trees_to_files, commands)
        # for command in commands:
        #     write_trees_to_files(command)
    else:
        raise NotImplementedError("Mode {} is not implemented".format(mode))
        # pool.map(
        #     write_trees_to_files,
        #     zip(nicks, [collection_path] * len(nicks), [datasetdb] * len(nicks)),
        # )

    # elif mode == "xrootd":  # it did not complete when using Pool in xrootd mode
    #     for nick in nicks:
    #         write_trees_to_files([nick, collection_path, datasetdb])


# main function
if __name__ == "__main__":
    args = parse_arguments()
    collect_outputs(args.executable, args.cores, args.custom_workdir_path, args.mode)