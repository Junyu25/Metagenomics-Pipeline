import os
import shutil 
import argparse
import subprocess
import pandas as pd
#from Bio import SeqIO
#from shutil import copyfile
from itertools import repeat
from multiprocessing import Pool, freeze_support

## Generate manifest Table
def manifestGen(InDir, OutDir):
    SampleID = []
    df = pd.DataFrame()
    for subdir, dirs, files in os.walk(InDir):
        for file in files:
            filePath = os.path.join(subdir, file)
            if file.endswith(r1_end) and os.path.getsize(filePath) > 0:
                R1 = os.path.join(subdir, file)
                R2 = os.path.join(subdir, file.replace(r1_end, r2_end))
                if os.path.exists(R2) and os.path.getsize(filePath) > 0:
                    SampleID = file.replace(r1_end, "")
                    df = df.append({"SampleID": SampleID, "R1":R1, "R2":R2},ignore_index=True)
    return df # return pair end dataframe

## Kneadata pair end tirmming
'''
kneaddata -i R1 -i R2 -o OutDir -db dbPath --output-prefix prefix --threads 16 --remove-intermediate-output --run-fastqc-start --run-fastqc-end --cat-final-output
'''
def RunKneadata(R1, R2, db, trimmomaticPath, prefix, OutDir, threads):
    cmd = "kneaddata -i " + R1 + " -i " + R2 + " --reference-db " + db + " -o " + OutDir + " --output-prefix " + prefix + \
        " --thread " + str(threads) + " --trimmomatic " + trimmomaticPath + \
        " --remove-intermediate-output --run-fastqc-start --run-fastqc-end --cat-final-output"
    subprocess.call(cmd, shell=True)
## Run RunKneadata in parallel
def RunKneadataParallel(R1List, R2List, db, trimmomaticPath, prefixList, OutDir, threads, jobs):
    pool = Pool(processes = jobs)
    pool.starmap(RunKneadata, zip(R1List, R2List, repeat(db), repeat(trimmomaticPath), prefixList, repeat(OutDir), repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()

def sumKneaddata(kneadataDir, OutDir):
    #Clean up the output directory (helps downstream commands) by moving the discarded sequences to a subfolder:
    contam_seq_dir = os.path.join(kneadataDir, "contam_seq")
    if os.path.exists(contam_seq_dir) == 0:
        os.makedirs(contam_seq_dir, 0o777, True)
    #shutil.move(source, destination)
    subprocess.call("mv " + kneadataDir + "/*_contam*.fastq" + " " + contam_seq_dir, shell=True)
    #You can produce a logfile summarizing the kneadData output with this command:
    cmd = "kneaddata_read_count_table --input " +  kneadataDir + " --output " + os.path.join(OutDir, "kneaddata_read_counts.txt")
    subprocess.call(cmd, shell=True)

#os.system('perl %s/concat_paired_end.pl -p %s --no_R_match -o cat_reads kneaddata_out/*_paired_*.fastq' % (script_path, nnodes)) 
def CatReads(kneadataDir, OutDir, threads):
    cmd = "perl " + scriptPath+"/concat_paired_end.pl -p " + threads + " --no_R_match -o " + cat_reads_dir + " " + kneadataDir + "/*_paired_*.fastq"
    subprocess.call(cmd, shell=True)
def parseCatReads(InDir):
    prefixList = []
    fastqList = []
    for file in os.listdir(cat_reads_dir):
        filePath = os.path.join(InDir, file)
        if file.endswith(".fastq") and os.path.getsize(filePath) > 0:
            prefixList.append(file.replace(".fastq", ""))
            fastqList.append(filePath)
    return fastqList, prefixList

# set database dir or add it to parameter
def RunHumann3(fastq, prefix, OutDir, threads):
    cmd = "humann --input " + fastq + " --output " + os.path.join(OutDir, prefix) + " --threads " + str(threads)
def RunHumann3Parallel(fastqList, prefixList, OutDir, threads, jobs):
    pool = Pool(processes = jobs)
    pool.starmap(RunHumann3, zip(fastqList, prefixList, repeat(OutDir), repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()
#humann3_final_out
def SumHumann3(humann3Dir, OutDir):
    humann3_final_out_dir = os.path.join(OutDir, "humann3_final_out")
    if os.path.exists(humann3_final_out_dir) == 0:
        os.makedirs(humann3_final_out_dir, 0o777, True)
    # Merge individual sample data together
    subprocess.call("humann_join_tables -s --input " + humann3Dir + " --file_name pathabundance --output " + os.path.join(humann3_final_out_dir, "humann3_pathabundance.tsv"), shell=True)
    subprocess.call("humann_join_tables -s --input " + humann3Dir + " --file_name pathcoverage --output " + os.path.join(humann3_final_out_dir, "humann3_pathcoverage.tsv"), shell=True)
    subprocess.call("humann_join_tables -s --input " + humann3Dir + " --file_name genefamilies --output " + os.path.join(humann3_final_out_dir, "humann3_genefamilies.tsv"), shell=True)
    #Re-normalize gene family and pathway abundances (so that all samples are in units of copies per million).
    subprocess.call("humann_renorm_table --input " + os.path.join(humann3_final_out_dir, "humann3_pathabundance.tsv") + " --units cpm --output " + os.path.join(humann3_final_out_dir, "humann3_pathabundance_cpm.tsv"), shell=True)
    subprocess.call("humann_renorm_table --input " + os.path.join(humann3_final_out_dir, "humann3_pathcoverage.tsv") + " --units cpm --output " + os.path.join(humann3_final_out_dir, "humann3_pathcoverage_cpm.tsv"), shell=True)
    subprocess.call("humann_renorm_table --input " + os.path.join(humann3_final_out_dir, "humann3_genefamilies.tsv") + " --units cpm --output " + os.path.join(humann3_final_out_dir, "humann3_genefamilies_cpm.tsv"), shell=True)
    #Split HUMAnN2 output abundance tables in stratified and unstratified tables (stratified tables include the taxa associated with a functional profile).
    subprocess.call("humann_split_stratified_table --input " + os.path.join(humann3_final_out_dir, "humann3_pathabundance_cpm.tsv") + " --output " + humann3_final_out_dir, shell=True)
    subprocess.call("humann_split_stratified_table --input " + os.path.join(humann3_final_out_dir, "humann3_pathcoverage_cpm.tsv") + " --output " + humann3_final_out_dir, shell=True)
    subprocess.call("humann_split_stratified_table --input " + os.path.join(humann3_final_out_dir, "humann3_genefamilies_cpm.tsv") + " --output " + humann3_final_out_dir, shell=True)

def SumMetaphlan3(metaphlan3Dir, humann3Dir, OutDir):
    subprocess.call("cp " + humann3Dir + "/*/*/*metaphlan_bugs_list.tsv" + metaphlan3Dir, shell=True)
    subprocess.call("merge_metaphlan_tables.py " + metaphlan3Dir + "/*metaphlan_bugs_list.tsv > " + os.path.join(OutDir, "metaphlan3_merged.txt"), shell=True)
    subprocess.call("sed -i 's/_metaphlan_bugs_list//g' " + os.path.join(OutDir, "metaphlan3_merged.txt"), shell=True)

parser = argparse.ArgumentParser(description="Kneaddata pre clean reads for Humann3")
parser.add_argument('-i', '--input', dest='InDir', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OutDir', type=str, required=True,
                    help="the output path of reads")
parser.add_argument('-s', '--scriptPath', dest='scriptPath', type=str, required=False, default="/home/junyuchen/Lab/Metagenomics-Pipeline/Scripts",
                    help="the script path")
parser.add_argument('-d', '--database', dest='Kneaddata_db', type=str, required=False, default="/home/LDlab/Databases/GRCh38_PhiX_bowtie2_index/GRCh38_PhiX",
                    help="the database path")
parser.add_argument('-tp', '--trimmomaticPath', dest='trimmomaticPath', type=str, required=False, default="/home/junyuchen/Biosoft/anaconda3/share/trimmomatic-0.39-1",
                    help="the trimmomaticPath path")
parser.add_argument('-j', '--jobs', dest='jobs', type=str,  required=False, default='4',
                    help="the number of jobs run in parallel")
parser.add_argument('-t', '--threads', dest='threads', type=str, required=False, default='6',
                    help="the number of threads run for a job")
parser.add_argument('-F', '--sepF', dest='r1_end', type=str, required=False, default='_1.clean.fq.gz',
                    help="It is the surfix to recognize the forward info, default='_1.clean.fq.gz'.")
parser.add_argument('-R', '--sepR', dest='r2_end', type=str, required=False, default='_2.clean.fq.gz',
                    help="It is the surfix to recognize the reverse info, default='_2.clean.fq.gz'.")
args = parser.parse_args()


InDir = os.path.abspath(args.InDir)
OutDir = os.path.abspath(args.OutDir)
scriptPath = os.path.abspath(args.scriptPath)
Kneaddata_db = os.path.abspath(args.Kneaddata_db)
trimmomaticPath = os.path.abspath(args.trimmomaticPath)
jobs = int(args.jobs)
threads = int(args.threads)
r1_end = str(args.r1_end)
r2_end = str(args.r2_end)

## init out dir
if os.path.exists(OutDir) == 0:
    os.makedirs(OutDir, 0o777, True)
kneadataDir = os.path.join(OutDir, "kneadata_out")
if os.path.exists(kneadataDir) == 0:
    os.makedirs(kneadataDir, 0o777, True)
cat_reads_dir = os.path.join(OutDir, "cat_reads")
if os.path.exists(cat_reads_dir) == 0:
    os.makedirs(cat_reads_dir, 0o777, True)
humann3Dir = os.path.join(OutDir, "humann3_out")
if os.path.exists(humann3Dir) == 0:
    os.makedirs(humann3Dir, 0o777, True)
metaphlan3Dir = os.path.join(OutDir, "metaphlan3_out")
if os.path.exists(metaphlan3Dir) == 0:
    os.makedirs(metaphlan3Dir, 0o777, True)
## process manifest
df = manifestGen(InDir, OutDir)
df.to_csv(os.path.join(OutDir, "SampleTable.csv"), index = None)
prefixList = df["SampleID"].tolist()
R1List = df["R1"].tolist()
R2List = df["R2"].tolist()

RunKneadataParallel(R1List, R2List, Kneaddata_db, trimmomaticPath, prefixList, kneadataDir, threads, jobs)
sumKneaddata(kneadataDir, OutDir)
CatReads(kneadataDir, OutDir, threads)
fastqList, prefixList = parseCatReads(cat_reads_dir)
RunHumann3Parallel(fastqList, prefixList, OutDir, threads, jobs)
SumHumann3(humann3Dir, OutDir)
SumMetaphlan3(metaphlan3Dir, humann3Dir, OutDir)