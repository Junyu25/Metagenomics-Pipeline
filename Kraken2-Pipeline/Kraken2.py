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
import argparse
import subprocess
import pandas as pd
from itertools import repeat
from multiprocessing import Pool, freeze_support

#db = "/thinker/storage/org/SIAT/LD/kraken_database/minikraken_8GB_20200312"
#jobs
#threads


#Run Kraken2 in parallel
def RunKraken2Parallel(R1List, R2List, db, jobs, threads, outFileList):
    #numOfprocess = len(R1List)
    #pool = Pool(processes=numOfprocess)
    pool = Pool(processes=jobs)
    pool.starmap(RunKraken2, zip(R1List, R2List, repeat(db), repeat(threads), outFileList))
    pool.close()
    pool.join()
    pool.terminate()

'''
kraken2 --gzip-compressed --paired /thinker/globe/org/SIAT/LD/biodata/LD_lab/JunyuChen/data/PRJEB7949/ERR695627/ERR695627_1.fastq.gz /thinker/globe/org/SIAT/LD/biodata/LD_lab/JunyuChen/data/PRJEB7949/ERR695627/ERR695627_2.fastq.gz --db /thinker/storage/org/SIAT/LD/kraken_database/minikraken_8GB_20200312 --threads 20 --report sum.txt
'''
#RunKraken2
def RunKraken2(R1, R2, db, threads, OutFile):
    
    #cmd = "spades.py --isolate -1 " + R1 + " -2 " + R2 + " -o " + OutDir
    cmd = "kraken2 --gzip-compressed --paired " + R1 + " " + R2 + " -db " + db + " --threads " + str(threads) + " --report " + os.path.join(ouputDir, OutFile + ".txt")
    subprocess.call(cmd, shell=True)

#jobs
#threads
parser = argparse.ArgumentParser(description='Kraken2')
parser.add_argument('-i', '--input', dest='fileDir', type=str, required=True,
                    help="the path of the file path table")
parser.add_argument('-o', '--output', dest='OpDir', type=str, required=True,
                    help="the output path of manifest")
parser.add_argument('-j', '--jobs', dest='jobs', type=str, required=True,
                    help="the number of jobs run in parallel")
parser.add_argument('-t', '--threads', dest='threads', type=str, required=True,
                    help="the number of threads run for a job")
parser.add_argument('-db', '--database', dest='database', type=str, required=False, default="/thinker/storage/org/SIAT/LD/kraken_database/minikraken_8GB_20200312",
                    help="the path of kraken database")                    
args = parser.parse_args()

inputDir = str(args.fileDir)
ouputDir = os.path.abspath(args.OpDir)
jobs = int(args.jobs)
threads = int(args.threads)
db = os.path.abspath(args.database)

df = pd.read_table(inputDir)
outFileList = df["#SampleID"].tolist()
R1List = df["forward-absolute-filepath"].tolist()
R2List = df["reverse-absolute-filepath"].tolist()

RunKraken2Parallel(R1List, R2List, db, jobs, threads, outFileList)