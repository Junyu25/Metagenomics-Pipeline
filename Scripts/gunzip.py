import os
import argparse
import subprocess
import multiprocessing
import pandas as pd
from shutil import copyfile
from itertools import repeat
from multiprocessing import Pool, freeze_support



def RunKneadDataParallel(F_list, R_list):
    #numOfprocess = len(R1List)
    #pool = Pool(processes=numOfprocess)
    pool = Pool(processes=20)
    arguments = zip(F_list, R_list)
    pool.map(RunKneadData, iterable=arguments)
    pool.close()
    pool.join()
    pool.terminate()

def RunGunzip(F_list, R_list):
    cmd = "gunzip -k "
    subprocess.call(cmd, shell=True)

def RunCpF(F_list, ouputDir):
    print(F_list)
    cmd = "cp " + F_list + " " + ouputDir 
    subprocess.call(cmd, shell=True)
def RunCpR(R_list, ouputDir):
    print(F_list)
    cmd = "cp " + R_list + " " + ouputDir 
    subprocess.call(cmd, shell=True)

def RunCpFParallel(F_list, ouputDir):
    #numOfprocess = len(R1List)
    #pool = Pool(processes=numOfprocess)
    pool = Pool(processes=20)
    pool.starmap(RunCpF, zip(F_list, repeat(ouputDir)))
    pool.close()
    pool.join()
    pool.terminate()
    
parser = argparse.ArgumentParser(description='RunKneadData')
parser.add_argument('-i', '--input', dest='fileDir', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OpDir', type=str, required=True,
                    help="the output path of reads")
args = parser.parse_args()

inputfile = str(args.fileDir)
ouputDir = os.path.abspath(args.OpDir)

#print(inputfile)
df = pd.read_table(inputfile)
F_list = df["forward-absolute-filepath"].tolist() 
R_list = df["reverse-absolute-filepath"].tolist() 
#print(F_list)
RunCpFParallel(F_list, ouputDir)
