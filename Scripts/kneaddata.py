import os
import argparse
import subprocess
import multiprocessing
import pandas as pd
from shutil import copyfile
from itertools import repeat
from multiprocessing import Pool, freeze_support



def RunKneadDataParallel(F_list, R_list, outputDir, refdb_path, trimmomatic_path, nnodes, njobs):
    #numOfprocess = len(R1List)
    #pool = Pool(processes=numOfprocess)
    pool = Pool(processes=njobs)
    #arguments = zip(F_list, R_list)
    #pool.map(RunKneadData, iterable=arguments)
    pool.starmap(RunKneadData, zip(F_list, R_list, repeat(outputDir), repeat(refdb_path), repeat(trimmomatic_path), repeat(nnodes)))
    pool.close()
    pool.join()
    pool.terminate()

def RunKneadData(F_list, R_list, outputDir, refdb_path, trimmomatic_path, nnodes):
    cmd = "kneaddata -i " + F_list + " -i " + R_list + " -o " + outputDir + " -db " + refdb_path + " --trimmomatic " + trimmomatic_path + " -t " + str(nnodes) + " --run-bmtagger --remove-intermediate-output"
    subprocess.call(cmd, shell=True)


parser = argparse.ArgumentParser(description='Run KneadData in Parallel')
parser.add_argument('-i', '--input', dest='fileDir', type=str, required=True,
                    help="the path of the PathTable")
parser.add_argument('-o', '--output', dest='OpDir', type=str, required=True,
                    help="the output path of reads")
parser.add_argument('-n', '--node', dest='node', type=str, required=False, default='2',
                        help="the number of nodes to request")
parser.add_argument('-j', '--jobs', dest='jobs', type=str, required=False, default='8',
                    help="The number of jobs run parallell in humann2 step. Default is '8'; It is bounded by the total number of memory available. Each job should have 16GB memory")
parser.add_argument('-t', '--tpath', dest='tp', type=str, required=False, default='/home/LDlab/anaconda3/envs/humann2/share/trimmomatic-0.39-1/',
                    help="The path of the trimmomatic ref files. Default is /home/LDlab/anaconda3/envs/humann2/share/trimmomatic-0.39-1/")
parser.add_argument('-d', '--dpath', dest='dp', type=str, required=False, default='/home/yxtan/ref_databases/GRCh38_PhiX_bowtie2_index/GRCh38_PhiX',
                    help="The path of the ref database. Default is /home/yxtan/ref_databases/GRCh38_PhiX_bowtie2_index/GRCh38_PhiX")

args = parser.parse_args()

inputfile = os.path.abspath(args.fileDir)
outputDir = os.path.abspath(args.OpDir)
nnodes = int(args.node)
njobs = int(args.jobs)
trimmomatic_path = os.path.abspath(args.tp)
refdb_path = os.path.abspath(args.dp)


print(inputfile)
df = pd.read_table(inputfile)
F_list = df["forward-absolute-filepath"].tolist() 
R_list = df["reverse-absolute-filepath"].tolist() 
print(F_list)

RunKneadDataParallel(F_list, R_list, outputDir, refdb_path, trimmomatic_path, nnodes, njobs)