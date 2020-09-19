import os
import sys
import argparse
import pandas as pd

def PreData(inputDir, outputDir, r1, r2):
    for subdir, dirs, files in os.walk(inputDir):
        f1 = ''
        r2 = ''
        for fastq in files:
            if fastq.endswith(r1):
                #global f1 
                f1 = os.path.join(subdir, fastq)
                #print(f1)
                if os.path.exists(os.path.join(subdir, fastq.replace(r1, r2))):
                    r2 = os.path.join(subdir, fastq.replace(r1, r2))
                    #print(r2)
                    print(fastq.split("_")[0])
                    path = path.append({'#SampleID':fastq.split("_")[0], "forward-absolute-filepath":str(f1), "reverse-absolute-filepath":str(r2)}, ignore_index=True)
    path.to_csv(os.path.join(outputDir, "PathTable.tsv"), sep="\t", index=False, encoding = "utf-8")

parser = argparse.ArgumentParser(description='Kraken2')
parser.add_argument('-i', '--input', dest='fileDir', type=str, required=True,
                    help="the path of the file path table")
parser.add_argument('-o', '--output', dest='OpDir', type=str, required=True,
                    help="the output path of manifest")
parser.add_argument('-r1', '--r1', dest='r1_prefix', type=str, required=False, default='_1.clean.fq.gz',
                    help="the r1 prefix")
parser.add_argument('-r2', '--r2', dest='r2_prefix', type=str, required=False, default='_2.clean.fq.gz',
                    help="the r2 prefix")                    
args = parser.parse_args()

inputDir = os.path.abspath(args.fileDir)
outputDir = os.path.abspath(args.OpDir)
r1 = str(args.r1_prefix)
r2 = str(args.r2_prefix)

PreData(inputDir, outputDir, r1, r2)