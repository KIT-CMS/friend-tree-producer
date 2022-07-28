import os, json, logging

logger = logging.getLogger("job_managment")
from streampaths import *
import ROOT as r
import shutil
import subprocess
from multiprocessing import Pool

r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")


def write_trees_to_files(info):
    nick = info[0]
    collection_path = info[1]
    db = info[2]
    era = info[3]
    channel = info[4]
    # from the nick, remove the last _#number
    nick = nick.replace("_{}".format(nick.split("_")[-1]), "")
    print "Copying trees for %s" % nick
    print("Nick: {}, collection_path: {}, db: {}, era {}, channel {}".format(nick, collection_path, db, era, channel))
    # the outputpath follows the pattern era/channel/nick/files.root
    nick_path = os.path.join(collection_path, era, channel, nick)
    if not os.path.exists(nick_path):
        os.makedirs(nick_path)
    outputfile = os.path.join(nick_path, nick + ".root")
    # now merge all files in db["files"] into the output file using hadd
    if len(db["files"]) > 1:
        print "Merging %s" % nick
        hadd_cmd = ["hadd", outputfile]
        # unpack the filenames to total paths
        for filepath in db["files"]:
            full_filepath = os.path.join(collection_path, filepath)
            hadd_cmd.append(full_filepath)
        print(" ".join(hadd_cmd))
        subprocess.call(hadd_cmd)
    else:
        print "Copying %s" % nick
        try:
            shutil.copy(db["files"][0], outputfile)
        except IOError:
            print("Could not copy %s to %s" % (db["files"][0], outputfile))


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
            filename = "_".join([nick, folder, str(first), str(last)]) + ".root"
            filepath = os.path.join(workdir_path, nick, filename)
            if nick not in datasetdb:
                datasetdb[nick] = {}
                if era not in datasetdb[nick]:
                    datasetdb[nick][era] = {}
                    if channel not in datasetdb[nick][era]:
                        datasetdb[nick][era][channel] = {}
                        datasetdb[nick][era][channel]["files"] = []
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
            datasetdb[nick][era][channel]["files"].append(filepath)


    if mode == "local":
        nicks = sorted(datasetdb)
        for nick in nicks:
            for era in datasetdb[nick]:
                for channel in datasetdb[nick][era]:
                    write_trees_to_files([nick, collection_path, datasetdb[nick][era][channel], era, channel])
        # pool = Pool(cores)
        # pool.map(
        #     write_trees_to_files,
        #     zip(nicks, [collection_path] * len(nicks), [datasetdb] * len(nicks)),
        # )

    # elif mode == "xrootd":  # it did not complete when using Pool in xrootd mode
    #     for nick in nicks:
    #         write_trees_to_files([nick, collection_path, datasetdb])
