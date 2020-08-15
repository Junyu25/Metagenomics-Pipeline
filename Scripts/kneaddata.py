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

def RunKneadData(F_list, R_list):
    cmd = "kneaddata -i " + F_list + " -i " + R_list + " -o /home/junyuchen/Lab/MetaBGC/1520/1522-Out/pair \
    -db /home/junyuchen/Databases/GRCh38_PhiX_bowtie2_index/GRCh38_PhiX \
    --trimmomatic /home/junyuchen/Biosoft/anaconda3/envs/humann2/share/trimmomatic-0.39-1/ \
    -t 4 --trimmomatic-options \"SLIDINGWINDOW:4:20 MINLEN:50\" \
    --bowtie2-options \"--very-sensitive --dovetail\" --remove-intermediate-output"
    subprocess.call(cmd, shell=True)


parser = argparse.ArgumentParser(description='RunKneadData')
parser.add_argument('-i', '--input', dest='fileDir', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OpDir', type=str, required=True,
                    help="the output path of reads")
args = parser.parse_args()

inputfile = str(args.fileDir)
ouputDir = os.path.abspath(args.OpDir)

print(inputfile)
df = pd.read_table(inputfile)
F_list = df["forward-absolute-filepath"].tolist() 
R_list = df["reverse-absolute-filepath"].tolist() 
print(F_list)

RunKneadDataParallel(F_list, R_list)