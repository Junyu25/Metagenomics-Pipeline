import os
import sys
import pandas as pd

path = pd.DataFrame()
#!input the full path of project
project = sys.argv[1]
f1 = ''
r2 = ''
#project = "/home/junyuchen/Lab/Meta-Analysis/DataDownload/data/PRJEB23207"
for subdir, dirs, files in os.walk(project):
    for fastq in files:
        if fastq.endswith("1.fastq.gz"):
            #global f1 
            f1 = os.path.join(subdir, fastq)
            print(f1)
        elif fastq.endswith("2.fastq.gz"):
            #global r2 
            r2 = os.path.join(subdir, fastq)
            print(r2)
    path = path.append({'#SampleID':str(subdir.split("/")[-1]), "forward-absolute-filepath":str(f1), "reverse-absolute-filepath":str(r2)}, ignore_index=True)

path.to_csv("PathTable.tsv", sep="\t", index=False, encoding = "utf-8")
    