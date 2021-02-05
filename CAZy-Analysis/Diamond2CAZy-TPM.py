import os
import argparse
import subprocess
import pandas as pd
from itertools import repeat
from multiprocessing import Pool, freeze_support

'''
diamond blastx --query /home/junyuchen/Lab/Liuhongbin/result-2020-08/cat_reads/LN1-200410-P7-LD200408-0003-GTCCTAAG.fastq --evalue 1.0 --threads 16 --max-target-seqs 1 --outfmt 6 --db /home/junyuchen/Lab/Custom-DataBase/CAZy/CAZy --out /home/junyuchen/Lab/Custom-DataBase/CAZy/Test/Test2.tsv
'''
def RunDiamondDirect(inputDir, ouputDir):
    fastaList = []
    outFileList = []
    rpkmList = []
    for subdir, dirs, files in os.walk(inputDir):
        fasta = ""
        for file in files:
            if file.endswith(".fastq"):
                fasta = os.path.join(subdir, file)
                fastaList.append(fasta)
                OutFile = os.path.join(ouputDir, file + "_blastp.tsv")
                outFileList.append(OutFile)
                #rpkmList.append(os.path.join(ouputDir, file + "_tpm.tsv"))
                #RunDiamond(fasta, db, threads, OutFile)
                calRPKM(OutFile, os.path.join(ouputDir, file + "_rpk.tsv"), os.path.join(ouputDir, file + "_tpm.tsv"))
    #RunDiamondParallel(fastaList, db, jobs, threads, outFileList)
    #calRPKMParallel(outFileList, rpkmList)


def RunDiamondParallel(fastaList, db, jobs, threads, outFileList):
    pool = Pool(processes=jobs)
    pool.starmap(RunDiamond, zip(fastaList, repeat(db), repeat(threads), outFileList))
    pool.close()
    pool.join()
    pool.terminate()

def RunDiamond(fasta, db, threads, OutFile):
    cmd = "diamond blastx -q " + fasta  + " -o " + OutFile + " --evalue 1.0 --max-target-seqs 1 --outfmt 6 --db " + db + " -p " + str(threads) 
    subprocess.call(cmd, shell=True)



def calRPKMParallel(BlastTsvFileList, outFileList):
    pool = Pool(processes=jobs)
    pool.starmap(calRPKM, zip(BlastTsvFileList, outFileList))
    pool.close()
    pool.join()
    pool.terminate()

def calRPKM(BlastTsvFile, rpk_out, OutFile):
    df = pd.read_table(BlastTsvFile, header=None)
    fileName = os.path.split(BlastTsvFile)[1].replace(".fastq_blastp.tsv", "")
    df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    idMap = pd.read_csv(idmap)
    df1 = df.loc[(df["length"] >= 25) & (df["pident"] >= 35)]
    df2 = df1["sseqid"].value_counts()
    idCount = pd.DataFrame()
    idCount["ID"] = df2.index
    df4 = pd.DataFrame(df2.tolist())
    idCount["Count"] = df4[0]
    result = pd.merge(idCount, idMap, on='ID')
    rpk = result["Count"]/((result["Len"]/1000))
    cpm = rpk * (1/sum(rpk)) * 1000000
    out = pd.DataFrame()
    out_rpk = pd.DataFrame()
    out["# Gene Family"] = result["ID"]
    out[fileName] = cpm
    #out rpk
    out_rpk["# Gene Family"] = result["ID"]
    out_rpk[fileName] = rpk
    
    out_rpk.to_csv(rpk_out, index=None, sep="\t")
    out.to_csv(OutFile, index=None, sep="\t")



parser = argparse.ArgumentParser(description='Run Diamond')
parser.add_argument('-i', '--input', dest='fileDir', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OpDir', type=str, required=True,
                    help="the output path of reads")
parser.add_argument('-d', '--database', dest='database', type=str,  required=False, default='/home/junyuchen/Lab/Custom-DataBase/CAZy/CAZyDB_0.9.25',
                    help="the database path")
parser.add_argument('-m', '--mapping', dest='mapping', type=str,  required=False, default='/home/junyuchen/Lab/Liuhongbin/CAZyidMapping.csv',
                    help="the ID Mapping file of len and ID of CAZy database path")
parser.add_argument('-j', '--jobs', dest='jobs', type=str,  required=False, default='4',
                    help="the number of jobs run in parallel")
parser.add_argument('-t', '--threads', dest='threads', type=str, required=False, default='6',
                    help="the number of threads run for a job")
args = parser.parse_args()

inputDir = str(args.fileDir)
ouputDir = os.path.abspath(args.OpDir)
db = os.path.abspath(args.database)
idmap = os.path.abspath(args.mapping)
jobs = int(args.jobs)
threads = int(args.threads)

RunDiamondDirect(inputDir, ouputDir)