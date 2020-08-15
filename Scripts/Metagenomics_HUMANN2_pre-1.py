#coding:utf-8
from __future__ import print_function
import glob,os,re,argparse
import pandas as pd

def m2_run_kneaddata(nnodes, F_list, R_list):

    os.system('mkdir kneaddata_out')
    os.system("parallel -j %s --link 'kneaddata -i {1} -i {2} -o kneaddata_out/ \
    -db /home/junyuchen/Databases/GRCh38_PhiX_bowtie2_index/GRCh38_PhiX --trimmomatic /home/junyuchen/Biosoft/anaconda3/envs/humann2/share/trimmomatic-0.39-1/ \
    -t 4 --trimmomatic-options \"SLIDINGWINDOW:4:20 MINLEN:50\" \
    --bowtie2-options \"--very-sensitive --dovetail\" --remove-intermediate-output' \
     ::: %s ::: %s" % (nnodes, F_list, R_list))
    #Clean up the output directory (helps downstream commands) by moving the discarded sequences to a subfolder:
    os.system('mkdir -p kneaddata_out/contam_seq')
    os.system('mv kneaddata_out/*_contam*.fastq kneaddata_out/contam_seq')
    #You can produce a logfile summarizing the kneadData output with this command:
    os.system('kneaddata_read_count_table --input kneaddata_out --output kneaddata_read_counts.txt')
    #2.2 Concatenate unstitched output 
    #os.system('perl %s/concat_paired_end.pl -p %s --no_R_match -o cat_reads kneaddata_out/*_paired_*.fastq' % (script_path, nnodes)) 

rd_dir = "/home/junyuchen/Lab/MetaBGC/1520/1522-pair.tsv"
nnodes = 12
print(rd_dir)
df = pd.read_csv(rd_dir, sep='\t')
F_list = df["forward-absolute-filepath"].tolist() 
R_list = df["reverse-absolute-filepath"].tolist() 
print(R_list[1])
print(F_list[1])
m2_run_kneaddata(nnodes, F_list, R_list)






