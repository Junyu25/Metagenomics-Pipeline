import os
import sys
import pandas as pd

path = pd.DataFrame()
#!input the full path of project
project = sys.argv[1]
#project = "/home/junyuchen/Lab/Meta-Analysis/DataDownload/data/PRJEB23207"
for run in os.listdir(project):
    #panda run column
    print(run)
    
    for fastq in os.listdir(os.path.join(project, run)):
        if fastq.endswith("1.fastq.gz"):
            f1 = os.path.join(project, run, fastq)
            print(f1)
            #add to pandas run column
        elif fastq.endswith("2.fastq.gz"):
            r2 = os.path.join(project, run, fastq)
            print(r2)
    path = path.append({'#SampleID':str(run), "forward-absolute-filepath":str(f1), "reverse-absolute-filepath":str(r2)}, ignore_index=True)

path.to_csv(project+"/PathTable.tsv", sep="\t", index=False, encoding = "utf-8")
    