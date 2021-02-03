import os
import argparse
import subprocess
import pandas as pd
#from Bio import SeqIO
from shutil import copyfile
from itertools import repeat
from multiprocessing import Pool, freeze_support

## Kneadata pair end tirmming
'''
kneaddata -i R1 -i R2 -o OutDir -db dbPath --output-prefix prefix --threads 16 --remove-intermediate-output --run-fastqc-start --run-fastqc-end --cat-final-output
'''
def RunKneadata(R1, R2, prefix, OutDir, threads):
    cmd = "kneaddata -i " + R1 + " -i " + R2 + " -o " + OutDir + " --output-prefix " + prefix + " --thread " + str(threads) + \
    " --remove-intermediate-output --run-fastqc-start --run-fastqc-end --cat-final-output"
    subprocess.call(cmd, shell=True)
## Run RunKneadata in parallel
def RunKneadataParallel(R1List, R2List, prefixList, OutDir, threads, jobs):
    pool = Pool(processes = jobs)
    pool.starmap(RunKneadata, zip(R1List, R2List, prefixList, repeat(OutDir), repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()

def sumKneaddata(kneadataDir, OutDir):
    cmd = "kneaddata_read_count_table --input " +  kneadataDir + " --output " os.path.join(OutDir, "kneaddata_read_counts.txt")
    subprocess.call(cmd, shell=True)

#os.system('perl %s/concat_paired_end.pl -p %s --no_R_match -o cat_reads kneaddata_out/*_paired_*.fastq' % (script_path, nnodes)) 

# set database dir or add it to parameter

def RunHumann3(fastq, prefix, OutDir, threads):
    cmd = "humann --input " + fastq + " --output " + os.path.join(OutDir, prefix) + " "