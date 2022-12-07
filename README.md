# friend-tree-producer
Software setup to produce additional quantities on top of ntuples in sync format with batch-system support

## Setting up the software
To setup the environment, please download the checkout and execute the checkout script as follows:

```bash
wget https://raw.githubusercontent.com/KIT-CMS/friend-tree-producer/ul_crown_ffs/scripts/checkout.sh
bash checkout.sh
```
## General setup of friend tree executables
Several external software relies on `C++` or `python` executables to create friend trees.
The `C++` executables should be placed into the [bin](https://github.com/KIT-CMS/friend-tree-producer/tree/master/bin)
folder with corresponding entry in the [BuildFile.xml](https://github.com/KIT-CMS/friend-tree-producer/tree/master/bin/BuildFile.xml),
whereas the `python` executables should be placed into the [scripts](https://github.com/KIT-CMS/friend-tree-producer/tree/master/scripts) folder.
After compilation within CMSSW, these executables should be accessible within the command line.

A common set of command line options should be used for these executables in order to process them equally through the batch system job management.
This set of command line options should be:

 * `--input`: path to the input ntuple to be used
 * `--folder`: systematics folder to be accessed within the file defined with the option above
 * `--tree`: tree name to be accessed within the folder defined with the option above
 * `--first_entry`: first entry to be processed within the tree
 * `--last_entry`: last entry to be processed within the tree
 * `--input-friends`: (optional) list of friend tree files

Asuming the absolute path to the input file is `/path/to/the/<input>.root`, the path to the output file starting from the current directory should read:
`<input>/<input>_<folder>_<first_entry>_<last_entry>.root`.

### Example command with SVFit executable

```bash
SVFit --input DY1JetsToLLM50_RunIIFall17MiniAODv2_PU2017_13TeV_MINIAOD_madgraph-pythia8_v1.root --folder mt_nominal --last_entry 999 --tree ntuple --first_entry 0
```
### Example command with FakeFactors executable
```bash
FakeFactors --input /ceph/olavoryk/smhtt/2022_07_v6/ntuples/2018/SingleMuon_Run2018A-UL2018/mt/SingleMuon_Run2018A-UL2018_0.root --ff_file data/fake_factors_2018_mt_nmssm.json --channel mt --era 2018 --last_entry 99 --folder nominal
```

## Job management for condor batch systems
The main script to submit jobs is [job_management.py](https://github.com/KIT-CMS/friend-tree-producer/blob/master/scripts/job_management.py). Following options are available:

 * `--executable`: Executable to be used for friend tree creation ob the batch system. Currently, only the choice for `SVFit` and `MELA` available.
 * `--batch_cluster`: Batch system cluster to be used. Currently available choices: `naf`, `etp` and `lxplus`. The templates for the `.jdl` files can be found in the [data](https://github.com/KIT-CMS/friend-tree-producer/tree/master/data) folder.
 * `--command`: Command to be done by the job manager. The `submit` command preprares a submission `.jdl` file for the desired condor batch system. The `collect` command merges the produced outputs to a single output file.
 * `--input_ntuples_directory`: Directory where the input files can be found. The file structure in the directory should match `*/*.root` wildcard.
 * `--friend_ntuples_directories`: List of directories where the friend files can be found. The file structure in the directory should match the one of the base ntuples. Channel dependent parts of the path can be inserted like /commonpath/{et:et_folder,mt:mt_folder,tt:tt_folder}/commonpath. If channel dependecies are given, this option is only forwarded to job executables for the respective channels.
 * `--events_per_job`: Event to be processed by each job.
 * `--walltime`: This option should be only set, if it is required by the batch cluster you are using. Currently, for the `etp` cluster.
 * `--cores`: Number of cores to be used for the collect command.

Please use also the `--help` command for this script to see the details of its execution.

### Example command with MELA executable for job submission
In order to submit the jobs to a batch system, you can use the following command as an example to start with:

```bash
job_management.py --executable MELA --input_ntuples_directory Full_2017_test_mt_11_05_2019/ --batch_cluster etp --command submit --events_per_job 100000 --walltime 3600
```

This command will create (if needed) a directory `MELA_executable` and will make the following printout:

```bash
To run the condor submission, execute the following:

cd <PWD>/MELA_workdir; condor_submit <PWD>/MELA_workdir/condor_MELA.jdl

```

The expression `<PWD>` stands here in this `README.md` for the current directory you use.

### Example command with MELA executable to collect the job outputs
To collect the outputs of the command before, you can execute the following command:

```bash
job_management.py --executable MELA --input_ntuples_directory Full_2017_test_mt_11_05_2019/ --batch_cluster etp --command collect --events_per_job 100000 --walltime 3600 --cores 10
```

This command will create another folder at `<PWD>/MELA_workdir/MELA_collected/` and put there the merged outputs matching a `*/*.root` structure.

### Check and resubmit failed condor jobs
The collect command does not end successfully, if some condor jobs did not finish successfully. To create a configuration to resubmit the crashed jobs, run 

```bash
job_management.py --executable MELA --command check --custom_workdir_path "/path/to/workdir" --logfile path/to/logfile 
```
This will create a `*resubmit.jdl` file which only submits jobs that did not end successfully
