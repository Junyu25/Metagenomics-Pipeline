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

#db = "/thinker/storage/org/SIAT/LD/kraken_database/minikraken_8GB_20200312"


'''
kraken2 --gzip-compressed --paired /thinker/globe/org/SIAT/LD/biodata/LD_lab/JunyuChen/data/PRJEB7949/ERR695627/ERR695627_1.fastq.gz /thinker/globe/org/SIAT/LD/biodata/LD_lab/JunyuChen/data/PRJEB7949/ERR695627/ERR695627_2.fastq.gz --db /thinker/storage/org/SIAT/LD/kraken_database/minikraken_8GB_20200312 --threads 20 --report sum.txt
'''
#RunKraken2
def RunKraken2(prefix, R1, R2, database, OutDir, threads):
    
    #cmd = "spades.py --isolate -1 " + R1 + " -2 " + R2 + " -o " + OutDir
    cmd = "kraken2 --paired " + R1 + " " + R2 + " -db " + database + " --threads " + str(threads) + " --report " + os.path.join(OutDir, prefix + ".txt")
    subprocess.call(cmd, shell=True)
#Run Kraken2 in parallel
def RunKraken2Parallel(prefixList, R1List, R2List, database, OutDir, threads, jobs):
    #numOfprocess = len(R1List)
    #pool = Pool(processes=numOfprocess)
    pool = Pool(processes=jobs)
    pool.starmap(RunKraken2, zip(prefixList, R1List, R2List, repeat(database), repeat(OutDir), repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()


## Argument Parser
parser = argparse.ArgumentParser(description="kraken2")
parser.add_argument('-i', '--input', dest='InDir', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OutDir', type=str, required=True,
                    help="the output path of kraken2 result")
parser.add_argument('-d', '--database', dest='kraken2_db', type=str, required=False, default="/home/junyuchen/3-Resources/databases/k2_standard_20201202",
                    help="the database path")
parser.add_argument('-j', '--jobs', dest='jobs', type=str,  required=False, default='4',
                    help="the number of jobs run in parallel")
parser.add_argument('-t', '--threads', dest='threads', type=str, required=False, default='6',
                    help="the number of threads run for a job")
parser.add_argument('-F', '--sepF', dest='r1_end', type=str, required=False, default='_paired_1.fastq',
                    help="It is the surfix to recognize the kneaddata forward info, default='_paired_1.fastq'.")
parser.add_argument('-R', '--sepR', dest='r2_end', type=str, required=False, default='_paired_2.fastq',
                    help="It is the surfix to recognize the kneaddata reverse info, default='_paired_2.fastq'.")
args = parser.parse_args()

## Argument
InDir = os.path.abspath(args.InDir)
OutDir = os.path.abspath(args.OutDir)
kraken2_db = os.path.abspath(args.kraken2_db)
jobs = int(args.jobs)
threads = int(args.threads)
r1_end = str(args.r1_end)
r2_end = str(args.r2_end)

## init output dir
if os.path.exists(OutDir) == 0:
    os.makedirs(OutDir, 0o777, True)
kraken2Dir = os.path.join(OutDir, "kraken2_out")
if os.path.exists(kraken2Dir) == 0:
    os.makedirs(kraken2Dir, 0o777, True)

## process manifest
df = manifestGen(InDir, OutDir)
prefixList = df["SampleID"].tolist()
R1List = df["R1"].tolist()
R2List = df["R2"].tolist()
df.to_csv(os.path.join(OutDir, "SampleTable.csv"), index = None)

RunKraken2Parallel(prefixList, R1List, R2List, kraken2_db, OutDir, threads, jobs)