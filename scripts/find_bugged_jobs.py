import os
import time
import argparse
from rich import print
from rich.progress import (
    BarColumn,
    TimeRemainingColumn,
    Progress,
)


def get_sample(jobid):
    sample = ""
    with open(submission_config, "r") as f:
        for line in f:
            line = line.rstrip()
            if ' -eq {} ];'.format(jobid) in line:
                job_command = next(f)
                sample = [x for x in job_command.split(" ")
                          if "root://" in x][0]
                return sample.split("/")[-2]


parser = argparse.ArgumentParser(
    description=
    'Small script to merge artus outputs from local or xrootd resources using multiprocessing.'
)
parser.add_argument('--workdir',
                    required=True,
                    help='Path to the friend tree workdir')
parser.add_argument('--submit-file',
                    required=True,
                    help='Path to the shell script for the submission')
parser.add_argument('--error-string',
                    type=str,
                    help='String that should be checked for each job',
                    default='Error in <TNetXNGFile::TNetXNGFile>')
parser.add_argument('--stdout',
                    action='store_true',
                    help='if set, the stdout and the stderr file of the job is checked')

args = parser.parse_args()

progress = Progress(
    "[progress.description]{task.description}",
    BarColumn(),
    "{task.completed} of {task.total}",
    "[progress.percentage]{task.percentage:>3.0f}%",
    TimeRemainingColumn(),
)

basepath = args.workdir
if "gc_workdir" not in basepath:
    basepath = os.path.join(basepath, "gc_workdir/output")
submission_config = args.submit_file
samplelist = []
with progress:
    totaljobs = len(next(os.walk(basepath))[1])
    task = progress.add_task("[red]Check folders...",
                             total=totaljobs + 1,
                             expand=True)
    print("Processing {} logfiles ".format(totaljobs + 1))
    for subdir, dirs, files in os.walk(basepath):
        for i, filename in enumerate(files):
            filepath = subdir + os.sep + filename
            if args.stdout:
                checkstring = ".std"
            else:
                checkstring = ".stderr"
            if filepath.endswith(checkstring):
                if args.error_string in open(filepath).read():
                    print("Error found in log {}".format(filepath))
                    jobid = filepath.split("/")[-2].split("_")[1]
                    samplename = get_sample(jobid)
                    if samplename not in samplelist:
                        samplelist.append(samplename)
                        print("Buggy sample: {}".format(samplename))
        progress.update(task, advance=1)
if (len(samplelist) > 0):
    print(
        "[bold red]Error![/bold red]"
    )
    print("[bold red] Missing samples are: [/bold red]")
    samplelist.sort()
    print(set(samplelist))
else:
    print(
        "[bold green] Well done ! All Samples should be correct ! [/bold green]"
    )
