import os
import argparse
import subprocess
import pandas as pd
from itertools import repeat
from multiprocessing import Pool, freeze_support
'''
def RunBWAindexParallel(InFileList):
    pool = Pool(processes=jobs)
    pool.starmap(RunBWAindex, zip(InFileList))
    pool.close()
    pool.join()
    pool.terminate()
'''    
def RunBWAParallel(indexList, R1List, R2List, jobs, threads, outFileList):
    pool = Pool(processes=jobs)
    pool.starmap(RunBWA, zip(repeat(indexList),R1List, R2List, repeat(threads), outFileList))
    pool.close()
    pool.join()
    pool.terminate()

def RunBBMapParallel(InFileList, OutFileList):
    pool = Pool(processes=jobs)
    pool.starmap(RunBBMap, zip(InFileList, OutFileList))
    pool.close()
    pool.join()
    pool.terminate()

def RunCatRPKMParallel(InFileList, OutFileList):
    pool = Pool(processes=jobs)
    pool.starmap(RunCatRPKM, zip(InFileList, OutFileList))
    pool.close()
    pool.join()
    pool.terminate()

def RunCatCPMParallel(InFileList, OutFileList):
    pool = Pool(processes=jobs)
    pool.starmap(RunCatCPM, zip(InFileList, OutFileList))
    pool.close()
    pool.join()
    pool.terminate()

def RunBWAindex(InFile):
    cmd = "bwa index " + InFile
    subprocess.call(cmd, shell=True)

def RunBWA(index, R1, R2, threads, OutFile):
    cmd = "bwa mem -t " + str(threads) + " "+ index + " " + R1 + " " + R2 + " > " + os.path.join(outputDir, OutFile + ".sam")
    subprocess.call(cmd, shell=True)

def RunBBMap(InFile, OutFile):
    cmd = "pileup.sh in=" + InFile + " out=" + os.path.join(outputDir, OutFile + "_coverage.txt") + " rpkm=" + os.path.join(outputDir, OutFile + "_rpkm.txt")
    subprocess.call(cmd, shell=True)

def RunCatRPKM(InFile, OutFile):
    cmd = "cat " + InFile + " | grep -v '#' | cut -f1, 6 > " + os.path.join(outputDir, OutFile + "_gene_abundance_rpkm.txt")
    subprocess.call(cmd, shell=True)

def RunCatCPM(InFile, OutFile):
    cmd = "cat " + InFile + " | grep -v '#' | cut -f1,2,5 > " + os.path.join(outputDir, OutFile + "_reads.input.txt")
    subprocess.call(cmd, shell=True)


parser = argparse.ArgumentParser(description='Calculate CAZy gene coverage')
parser.add_argument('-i', '--input', dest='fileDir', type=str, required=True,
                    help="the path of the bwa index .fasta file path")
parser.add_argument('-d', '--kneaddata', dest='kneaddata', type=str, required=True,
                    help="the path of the keandata file path")
parser.add_argument('-o', '--output', dest='OpDir', type=str, required=True,
                    help="the output path of result")
parser.add_argument('-j', '--jobs', dest='jobs', type=str, required=True,
                    help="the number of jobs run in parallel")
parser.add_argument('-t', '--threads', dest='threads', type=str, required=True,
                    help="the number of threads run for a job")               
args = parser.parse_args()

index = os.path.abspath(args.fileDir)
inputDir = os.path.abspath(args.kneaddata)
outputDir = os.path.abspath(args.OpDir)
jobs = int(args.jobs)
threads = int(args.threads)


R1 = ""
R2 = ""
OutFile = ""
R1List = []
R2List = []
OutFileList = []
for files in os.listdir(inputDir):
    if files.endswith("_kneaddata_paired_1.fastq"):
        R1 = os.path.join(inputDir, files)
        R2 = os.path.join(inputDir, files[:-25]+"_kneaddata_paired_2.fastq")
        OutFile = files.replace("_kneaddata_paired_1.fastq", "")
        R1List.append(R1)
        R2List.append(R2)
        OutFileList.append(OutFile)
RunBWAParallel(index,R1List, R2List, jobs, threads, OutFileList)

BBMapIn = ""
BBMapOut = ""
BBMapInList = []
BBMapOutList = []
for files in os.listdir(outputDir):
    if files.endswith(".sam"):
        BBMapIn = os.path.join(outputDir, files)
        BBMapOut = files.replace(".sam", "")
        BBMapInList.append(BBMapIn)
        BBMapOutList.append(BBMapOut)
RunBBMapParallel(BBMapInList, BBMapOutList)

for file in os.listdir(outputDir):
    if file.endswith("_rpkm.txt"):
        InFile = os.path.join(outputDir, file)
        OutFile = files.replace("_rpkm.txt", "")
        RunCatRPKM(InFile, OutFile)
        RunCatCPM(InFile, OutFile)