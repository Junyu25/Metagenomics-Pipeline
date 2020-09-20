import os
import argparse
import subprocess
import pandas as pd
from itertools import repeat
from multiprocessing import Pool, freeze_support

#Run Spades on a directory
def RunSpadesDirectory(inputDir, ouputDir):
    R1 = ""
    R2 = ""
    outputFilePath = ""
    SpadesFilePath = ""
    SpadesOutDir = ""
    R1List = []
    R2List = []
    outFileList = []
    SpadesFileList = []
    SpadesOutList = []
    for subdir, dirs, files in os.walk(inputDir):
        for file in files:
            if file.endswith("_1_kneaddata_paired_1.fastq"):
                R1 = os.path.join(subdir, file)
                R2 = os.path.join(subdir, file[:-7]+"2.fastq")
                print(R1)
                print(R2)
                R1List.append(R1)
                R2List.append(R2)
                sampleStr = os.path.splitext(file)[0].replace("_1_kneaddata_paired_1", "")
                outputFilePath = os.path.join(ouputDir, sampleStr)
                outFileList.append(outputFilePath)
                SpadesOutDir = os.path.join(ouputDir, sampleStr)
                SpadesOutList.append(SpadesOutDir)
                #make out dir for every run
                os.makedirs(os.path.join(ouputDir, sampleStr), 0o777, True)
    RunSpadesParallel(R1List, R2List, SpadesOutList)

#Run Spades in parallel
def RunSpadesParallel(R1List, R2List, outFileList):
    #numOfprocess = len(R1List)
    #pool = Pool(processes=numOfprocess)
    pool = Pool(processes=jobs)
    pool.starmap(RunSpades, zip(R1List, R2List, outFileList))
    pool.close()
    pool.join()
    pool.terminate()

#SPAdes Assembling
def RunSpades(R1, R2, OutFile):
    OutDir = os.path.join(outputDir, OutFile)
    os.makedirs(OutDir, 0o777, True)
    #cmd = "spades.py --isolate -1 " + R1 + " -2 " + R2 + " -o " + OutDir
    cmd = "spades.py --meta -1 " + R1 + " -2 " + R2 + " -o " + OutDir + " -t 16"
    subprocess.call(cmd, shell=True)


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

RunSpadesParallel(R1List, R2List, outFileList)