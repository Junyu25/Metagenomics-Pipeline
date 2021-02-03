#coding:utf-8
from __future__ import print_function
import glob,os,re,argparse
import pandas as pd


def m1_fastqc(ALL_list,nnodes):
    """
    :ALL_list: the string of all the paths of fastqc files '/'.
    :return: None
    """     
    if not os.path.isdir('fastqc_out'):
        print('Input folder is not exist; mkdir now.')
        os.system('mkdir fastqc_out')
        os.system('fastqc -t %s %s -o fastqc_out/' % (nnodes, ALL_list))
    else:
        print('fastqc step is skipped.')

def m2_run_kneaddata(nnodes, F_list, R_list, script_path):
    """
    :nnodes: the number of nodes to parallel
    :F_list: the string of all the paths of forward fastqc files '/'.
    :R_list: the string of all the paths of paired reverse fastqc files '/'.
    :script_path: the path of custome script, default is /home/yxtan/HUMANN2_SOP_scripts/
    :return: None
    """
    #2.1 Running KneadData
    #Use kneaddata to run pre-processing tools. First Trimmomatic is run to remove low quality sequences. Then Bowtie2 is run to screen out contaminant sequences. Below we are screening out reads that map to the human or PhiX genomes. Note KneadData is being run below on all unstitched FASTQ pairs with parallel, you can see our quick tutorial on this tool here. For a detailed breakdown of the options in the below command see this page. The forward and reverse reads will be specified by "_1" and "_2" in the output files, ignore the "R1" in each filename. Note that the \ characters at the end of each line are just to split the command over multiple lines to make it easier to read.
    #the order of quotes are extremely important here.
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

def m2_run_kneaddata_single(nnodes, F_list,script_path):
    """
    :nnodes: the number of nodes to parallel
    :F_list: the string of all the paths of forward fastqc files '/'.
    :script_path: the path of custome script, default is /home/yxtan/HUMANN2_SOP_scripts/
    :return: None
    """
    #2.1 Running KneadData
    #Use kneaddata to run pre-processing tools. First Trimmomatic is run to remove low quality sequences. Then Bowtie2 is run to screen out contaminant sequences. Below we are screening out reads that map to the human or PhiX genomes. Note KneadData is being run below on all unstitched FASTQ pairs with parallel, you can see our quick tutorial on this tool here. For a detailed breakdown of the options in the below command see this page. The forward and reverse reads will be specified by "_1" and "_2" in the output files, ignore the "R1" in each filename. Note that the \ characters at the end of each line are just to split the command over multiple lines to make it easier to read.
    #the order of quotes are extremely important here.
    os.system('mkdir kneaddata_out')
    os.system("parallel -j %s 'kneaddata -i {1} -o kneaddata_out/ \
    -db /home/junyuchen/Databases/GRCh38_PhiX_bowtie2_index/GRCh38_PhiX --trimmomatic /home/junyuchen/Biosoft/anaconda3/envs/humann2/share/trimmomatic-0.39-1/ \
    -t 4 --trimmomatic-options \"SLIDINGWINDOW:4:20 MINLEN:50\" \
    --bowtie2-options \"--very-sensitive --dovetail\" --remove-intermediate-output' \
     ::: %s " % (nnodes, F_list))
    #Clean up the output directory (helps downstream commands) by moving the discarded sequences to a subfolder:
    os.system('mkdir -p kneaddata_out/contam_seq')
    os.system('mv kneaddata_out/*_contam*.fastq kneaddata_out/contam_seq')
    #You can produce a logfile summarizing the kneadData output with this command:
    os.system('kneaddata_read_count_table --input kneaddata_out --output kneaddata_read_counts.txt')
    #2.2 Concatenate unstitched output 
    #os.system('mkdir cat_reads/')
    #os.system('mv kneaddata_out/*kneaddata.fastq cat_reads/') 




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='rd', type=str, required=True,
                        help="the tabular-table contains paths of the raw data")
    parser.add_argument('-n', '--node', dest='node', type=str, required=False, default='20',
                        help="the number of nodes to request")
    parser.add_argument('-sp', '--spath',dest='sp', type=str, required=False, default='/home/junyuchen/Lab/Meta-Analysis/Scripts/',
                        help="path of the custom scripts")
    parser.add_argument('-e', '--pair', dest='pair', type=str, required=False, default='True',
                        help="Is it pair-end seq data? Default is 'True'; Any other strings will be considered False")
    parser.add_argument('-j', '--jobs', dest='jobs', type=str, required=False, default='8',
                        help="The number of jobs run parallell in humann2 step. Default is '8'; It is bounded by the total number of memory available. Each job should have 16GB memory")
    args = parser.parse_args()
    print('Usage example with minimum parameters: python /home/yxtan/HUMANN2_SOP_scripts/Metagenomics_HUMANN2.py -i sample_table.txt -n 4')
    rd_dir = os.path.abspath(args.rd)
    script_path = os.path.abspath(args.sp)
    pair_end = args.pair
    nnodes = args.node
    njobs = args.jobs
    
    #需要检查输入的参数是否正确，主要是路径是否存在    
    if not os.path.isfile(rd_dir):
        print('Input sample table is not exist; Exit now.')
        exit(0)
    if not os.path.isdir(script_path):
        print('The folder of custom scripts is not exist; Exit now.')
        exit(0)

    #1. First Steps
    #1.1 Generate the list of samples
    #rd_dir = "example_data_file.txt"
    print(rd_dir)
    df = pd.read_csv(rd_dir, sep='\t')
    F_list = df["forward-absolute-filepath"].tolist() 
    R_list = df["reverse-absolute-filepath"].tolist() 
    print(R_list[1])
    print(F_list[1])
    #print(Counter)
    #1.2 Inspect read quality
    #if pair_end == 'True':
    #    m1_fastqc(ALL_list,nnodes)
    #else:
    #    m1_fastqc(F_list,nnodes)
    
    #2. Read Quality-Control and Contaminant Screens and connect to a long read
    #if pair_end == 'True':
    m2_run_kneaddata(nnodes, F_list, R_list, script_path)
    #else:
    #    m2_run_kneaddata_single(nnodes, F_list,script_path)
    
    #3. Determine Functions with HUMAnN2
    #m3_humann2(nnodes, script_path,njobs)






