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
                outFileList.append(os.path.join(ouputDir, file + "_blastp.tsv"))
                rpkmList.append(os.path.join(ouputDir, file + "_rpkm.tsv"))
    RunDiamondParallel(fastaList, db, jobs, threads, outFileList)
    calRPKMParallel(outFileList, rpkmList)


def RunDiamondParallel(fastaList, db, jobs, threads, outFileList):
    pool = Pool(processes=jobs)
    pool.starmap(RunDiamond, zip(fastaList, repeat(db), repeat(threads), outFileList))
    pool.close()
    pool.join()
    pool.terminate()

def RunDiamond(fasta, db, threads, OutFile):
    cmd = "diamond blastx --query " + fasta + " --evalue 1.0 --max-target-seqs 1 --outfmt 6 --db " + db + " --out " + OutFile + " " + str(threads) 
    subprocess.call(cmd, shell=True)



def calRPKMParallel(BlastTsvFileList, outFileList):
    pool = Pool(processes=jobs)
    pool.starmap(calRPKM, zip(BlastTsvFileList, outFileList))
    pool.close()
    pool.join()
    pool.terminate()

def calRPKM(BlastTsvFile, OutFile):
    df = pd.read_table(BlastTsvFile, header=None)
    df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    idMap = pd.read_csv(idmap, header=None)
    df1 = df.loc[(df["length"] >= 25) & (df["pident"] >= 35)]
    df2 = df1["sseqid"].value_counts()
    idCount = pd.DataFrame()
    idCount["ID"] = df2.index
    df4 = pd.DataFrame(df2.tolist())
    idCount["Count"] = df4[0]
    result = pd.merge(idCount, idMap, on='ID')
    rpkm = result["Count"]/((result["Len"]/1000) * (sum(result["Count"]))/1000000)
    result["rpkm"] = rpkm
    result.to_csv(OutFile)



parser = argparse.ArgumentParser(description='Run Diamond')
parser.add_argument('-i', '--input', dest='fileDir', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OpDir', type=str, required=True,
                    help="the output path of reads")
parser.add_argument('-d', '--database', dest='database', type=str,  required=False, default='/home/junyuchen/Lab/Custom-DataBase/CAZy/CAZy',
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