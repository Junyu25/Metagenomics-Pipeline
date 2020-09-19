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
def RunMegahitParallel(R1List, R2List, jobs, outFileList):
    #numOfprocess = len(R1List)
    #pool = Pool(processes=numOfprocess)
    pool = Pool(processes=jobs)
    pool.starmap(RunMegahit, zip(R1List, R2List,  outFileList))
    pool.close()
    pool.join()
    pool.terminate()

'''
kraken2 --gzip-compressed --paired /thinker/globe/org/SIAT/LD/biodata/LD_lab/JunyuChen/data/PRJEB7949/ERR695627/ERR695627_1.fastq.gz /thinker/globe/org/SIAT/LD/biodata/LD_lab/JunyuChen/data/PRJEB7949/ERR695627/ERR695627_2.fastq.gz --db /thinker/storage/org/SIAT/LD/kraken_database/minikraken_8GB_20200312 --threads 20 --report sum.txt
'''
#RunKraken2
def RunMegahit(R1, R2, OutFile):
    
    #cmd = "spades.py --isolate -1 " + R1 + " -2 " + R2 + " -o " + OutDir
    cmd = "megahit -1 " + R1 + " -2 " + R2 + " -o  " + os.path.join(outputDir, OutFile)
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
args = parser.parse_args()

inputDir = str(args.fileDir)
outputDir = os.path.abspath(args.OpDir)
jobs = int(args.jobs)


df = pd.read_table(inputDir)
outFileList = df["#SampleID"].tolist()
R1List = df["forward-absolute-filepath"].tolist()
R2List = df["reverse-absolute-filepath"].tolist()

RunMegahitParallel(R1List, R2List, jobs, outFileList)