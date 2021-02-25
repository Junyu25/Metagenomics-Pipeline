'''
Copyright {2020} Junyu Chen

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''

import os
#import shutil 
import argparse
import subprocess
import pandas as pd
#from Bio import SeqIO
#from shutil import copyfile
from itertools import repeat
from multiprocessing import Pool, freeze_support

## Generate manifest Table
def manifestGen(InDir, OutDir):
    SampleID = []
    df = pd.DataFrame()
    for subdir, dirs, files in os.walk(InDir):
        for file in files:
            filePath = os.path.join(subdir, file)
            if file.endswith(r1_end) and os.path.getsize(filePath) > 0:
                R1 = os.path.join(subdir, file)
                R2 = os.path.join(subdir, file.replace(r1_end, r2_end))
                if os.path.exists(R2) and os.path.getsize(filePath) > 0:
                    SampleID = file.replace(r1_end, "")
                    df = df.append({"SampleID": SampleID, "R1":R1, "R2":R2},ignore_index=True)
    return df # return pair end dataframe

## Kneadata pair end tirmming
'''
kneaddata -i R1 -i R2 -o OutDir -db dbPath --output-prefix prefix --threads 16 --remove-intermediate-output --run-fastqc-start --run-fastqc-end --cat-final-output
'''
def RunKneaddata(R1, R2, db, trimmomaticPath, prefix, OutDir, threads):
    cmd = "kneaddata -i " + R1 + " -i " + R2 + " --reference-db " + db + " -o " + OutDir + " --output-prefix " + prefix + \
        " --thread " + str(threads) + " --trimmomatic " + trimmomaticPath + \
        " --remove-intermediate-output --run-fastqc-start --run-fastqc-end --cat-final-output"
    subprocess.call(cmd, shell=True)
## Run Kneaddata in parallel
def RunKneaddataParallel(R1List, R2List, db, trimmomaticPath, prefixList, OutDir, threads, jobs):
    pool = Pool(processes = jobs)
    pool.starmap(RunKneaddata, zip(R1List, R2List, repeat(db), repeat(trimmomaticPath), prefixList, repeat(OutDir), repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()
## Get summary table of Kneadata result
def sumKneaddata(kneadataDir, OutDir):
    #Clean up the output directory (helps downstream commands) by moving the discarded sequences to a subfolder:
    contam_seq_dir = os.path.join(kneadataDir, "contam_seq")
    if os.path.exists(contam_seq_dir) == 0:
        os.makedirs(contam_seq_dir, 0o777, True)
    #shutil.move(source, destination)
    subprocess.call("mv " + kneadataDir + "/*_contam*.fastq" + " " + contam_seq_dir, shell=True)
    #You can produce a logfile summarizing the kneadData output with this command:
    cmd = "kneaddata_read_count_table --input " +  kneadataDir + " --output " + os.path.join(OutDir, "kneaddata_read_counts.txt")
    subprocess.call(cmd, shell=True)


## Argument Parser
parser = argparse.ArgumentParser(description="Kneaddata pre clean reads for Humann3")
parser.add_argument('-i', '--input', dest='InDir', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OutDir', type=str, required=True,
                    help="the output path of reads")
parser.add_argument('-s', '--scriptPath', dest='scriptPath', type=str, required=False, default="/home/junyuchen/Lab/Metagenomics-Pipeline/Scripts",
                    help="the script path")
parser.add_argument('-d', '--database', dest='Kneaddata_db', type=str, required=False, default="/home/LDlab/Databases/GRCh38_PhiX_bowtie2_index/GRCh38_PhiX",
                    help="the database path")
parser.add_argument('-tp', '--trimmomaticPath', dest='trimmomaticPath', type=str, required=False, default="/home/junyuchen/Biosoft/anaconda3/share/trimmomatic-0.39-1",
                    help="the trimmomaticPath path")
parser.add_argument('-j', '--jobs', dest='jobs', type=str,  required=False, default='4',
                    help="the number of jobs run in parallel")
parser.add_argument('-t', '--threads', dest='threads', type=str, required=False, default='6',
                    help="the number of threads run for a job")
parser.add_argument('-F', '--sepF', dest='r1_end', type=str, required=False, default='_1.clean.fq.gz',
                    help="It is the surfix to recognize the forward info, default='_1.clean.fq.gz'.")
parser.add_argument('-R', '--sepR', dest='r2_end', type=str, required=False, default='_2.clean.fq.gz',
                    help="It is the surfix to recognize the reverse info, default='_2.clean.fq.gz'.")
args = parser.parse_args()

## Argument
InDir = os.path.abspath(args.InDir)
OutDir = os.path.abspath(args.OutDir)
scriptPath = os.path.abspath(args.scriptPath)
Kneaddata_db = os.path.abspath(args.Kneaddata_db)
trimmomaticPath = os.path.abspath(args.trimmomaticPath)
jobs = int(args.jobs)
threads = int(args.threads)
r1_end = str(args.r1_end)
r2_end = str(args.r2_end)

## process manifest
df = manifestGen(InDir, OutDir)
prefixList = df["SampleID"].tolist()
R1List = df["R1"].tolist()
R2List = df["R2"].tolist()
df.to_csv(os.path.join(OutDir, "SampleTable.csv"), index = None)


## processing...
RunKneaddataParallel(R1List, R2List, Kneaddata_db, trimmomaticPath, prefixList, kneadataDir, threads, jobs)
sumKneaddata(kneadataDir, OutDir)