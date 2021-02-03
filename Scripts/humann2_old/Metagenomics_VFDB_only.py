#coding:utf-8
from __future__ import print_function
import glob,os,re,argparse

"""
#   Copyright {2019} Yuxiang Tan
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

#This script will run HUMANN2 automatically，mostly following the microbiome_helper workflow
#The output folder will be the current working folder (it is pointed by the workfd in qsub)
#Usage rule: python /home/yxtan/HUMANN2_SOP_scripts/Metagenomics_HUMANN2.py -i sample_table.txt -n 4'
#For single-end data, the sample table must have 3 columns as same as the pair-end data, but the third column will not be used and could be filled by any strings.

#This script starts from raw fqs before joined. No required rules for sample names, but information of each sample must be list in the tabular sample table.
#Library file requirement:
#HUMANN2 must be installed by anaconda, depending databases should be installed following the homepage instruction.
#kneaddata must be installed by anaconda


"""

def m1_fastqc(ALL_list):
    """
    :ALL_list: the string of all the paths of fastqc files '/'.
    :return: None
    """     
    if not os.path.isdir('fastqc_out'):
        print('Input folder is not exist; mkdir now.')
        os.system('mkdir fastqc_out')
        os.system('fastqc -t 4 %s -o fastqc_out/' % ALL_list)
    else:
        print('fastqc step is skipped.')

def m2_run_kneaddata(nnodes, F_list,R_list,script_path):
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
    os.system("parallel -j %s --link 'kneaddata -i {1} -i {2} -o kneaddata_out/ \
    -db /home/yxtan/ref_databases/GRCh38_PhiX_bowtie2_index/GRCh38_PhiX --trimmomatic /home/yxtan/miniconda2/envs/humann2/share/trimmomatic-0.39-1/ \
    -t 4 --trimmomatic-options \"SLIDINGWINDOW:4:20 MINLEN:50\" \
    --bowtie2-options \"--very-sensitive --dovetail\" --remove-intermediate-output' \
     ::: %s ::: %s" % (nnodes, F_list, R_list))
    #Clean up the output directory (helps downstream commands) by moving the discarded sequences to a subfolder:
    os.system('mkdir kneaddata_out/contam_seq')
    os.system('mv kneaddata_out/*_contam*.fastq kneaddata_out/contam_seq')
    #You can produce a logfile summarizing the kneadData output with this command:
    os.system('kneaddata_read_count_table --input kneaddata_out --output kneaddata_read_counts.txt')
    #2.2 Concatenate unstitched output 
    os.system('perl %s/concat_paired_end.pl -p %s --no_R_match -o cat_reads kneaddata_out/*_paired_*.fastq' % (script_path, nnodes)) 

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
    os.system("parallel -j %s 'kneaddata -i {1} -o kneaddata_out/ \
    -db /home/yxtan/ref_databases/GRCh38_PhiX_bowtie2_index/GRCh38_PhiX --trimmomatic /home/yxtan/miniconda2/envs/humann2/share/trimmomatic-0.39-1/ \
    -t 4 --trimmomatic-options \"SLIDINGWINDOW:4:20 MINLEN:50\" \
    --bowtie2-options \"--very-sensitive --dovetail\" --remove-intermediate-output' \
     ::: %s " % (nnodes, F_list))
    #Clean up the output directory (helps downstream commands) by moving the discarded sequences to a subfolder:
    os.system('mkdir kneaddata_out/contam_seq')
    os.system('mv kneaddata_out/*_contam*.fastq kneaddata_out/contam_seq')
    #You can produce a logfile summarizing the kneadData output with this command:
    os.system('kneaddata_read_count_table --input kneaddata_out --output kneaddata_read_counts.txt')
    #2.2 Concatenate unstitched output 
    os.system('mkdir cat_reads/')
    os.system('mv kneaddata_out/*kneaddata.fastq cat_reads/') 


#dowanlaod databases following http://huttenhower.sph.harvard.edu/humann2, humann2 will update the database path to $DIR
# DIR="/home/yxtan/HUMANN2_SOP_scripts/dateabases/"
# humann2_databases --download chocophlan full $DIR
# humann2_databases --download uniref uniref90_ec_filtered_diamond $DIR
# humann2_databases --download uniref uniref90_diamond $DIR
# humann2_databases --download uniref uniref50_diamond $DIR
# humann2_databases --download utility_mapping full $DIR
def m3_humann2(nnodes, script_path,njobs):
    """
    :nnodes: the number of nodes to parallel
    :script_path: the path of custome script, default is /home/yxtan/HUMANN2_SOP_scripts/
    :return: None
    """
    #3.1 Run HUMAnN2
    #Run humann2 with parallel to calculate abundance of UniRef90 gene families and MetaCyc pathways. If you are processing environment data (e.g. soil samples) the vast majority of the reads may not map using this approach. Instead, you can try mapping against the UniRef50 database (which you can point to with the --protein-database option).
    
    os.system('mkdir humann2_vf_out')
    os.system('humann2_config --update database_folders protein /home/junyuchen/Lab/Custom-DataBase/VFDB/DataBase/VFDB-2019-12-16/Diamod-VFDB_setB_pro')
    import math
    nthread=str(int(math.ceil(float(int(nnodes))/int(njobs))))
    #this is the most time consuming step
    os.system("parallel -j %s 'humann2 --threads %s --bypass-nucleotide-index --nucleotide-database /home/junyuchen/Lab/Custom-DataBase/VFDB/DataBase/VFDB-2019-12-20 --remove-temp-output --input {} --id-mapping /home/junyuchen/Lab/Custom-DataBase/VFDB/id-mapping-VF.tsv --output humann2_vf_out/{/.}' ::: cat_reads/*fastq" % (njobs,nthread))
    #3.2 Merge individual sample data together
    os.system('mkdir c')
    os.system('humann2_join_tables -s --input humann2_vf_out/ --file_name pathabundance --output humann2_vf_out/humann2_pathabundance.tsv')
    os.system('humann2_join_tables -s --input humann2_vf_out/ --file_name pathcoverage --output humann2_vf_out/humann2_pathcoverage.tsv')
    os.system('humann2_join_tables -s --input humann2_vf_out/ --file_name genefamilies --output humann2_vf_out/humann2_genefamilies.tsv')
    #3.3 Table output normalization
    #Re-normalize gene family and pathway abundances (so that all samples are in units of copies per million).
    os.system('humann2_renorm_table --input humann2_vf_out/humann2_pathabundance.tsv --units cpm --output humann2_vf_out/humann2_pathabundance_cpm.tsv')
    os.system('humann2_renorm_table --input humann2_vf_out/humann2_genefamilies.tsv --units cpm --output humann2_vf_out/humann2_genefamilies_cpm.tsv')
    #3.4 Separate out taxonomic contributions
    #Split HUMAnN2 output abundance tables in stratified and unstratified tables (stratified tables include the taxa associated with a functional profile).
    os.system('humann2_split_stratified_table --input humann2_vf_out/humann2_pathabundance_cpm.tsv --output humann2_final_out')
    os.system('humann2_split_stratified_table --input humann2_vf_out/humann2_genefamilies_cpm.tsv --output humann2_final_out')
    os.system('humann2_split_stratified_table --input humann2_vf_out/humann2_pathcoverage.tsv --output humann2_final_out')
    #3.5 Format STAMP function file
    #Convert unstratified HUMAnN2 abundance tables to STAMP format by changing header-line. These commands remove the comment character and the spaces in the name of the first column. Trailing descriptions of the abundance datatype are also removed from each sample's column name.
    os.system("sed 's/_Abundance-RPKs//g' humann2_vf_out/humann2_genefamilies_cpm_unstratified.tsv | sed 's/# Gene Family/GeneFamily/' > humann2_vf_out/humann2_genefamilies_cpm_unstratified.spf")
    os.system("sed 's/_Abundance//g' humann2_vf_out/humann2_pathabundance_cpm_unstratified.tsv | sed 's/# Pathway/Pathway/' > humann2_vf_out/humann2_pathabundance_cpm_unstratified.spf")
    #3.6 Extract MetaPhlAn2 taxonomic compositions
    #Since HUMAnN2 also runs MetaPhlAn2 as an initial step, we can use the output tables already created to get the taxonomic composition of our samples. First we need to gather all the output MetaPhlAn2 results per sample into a single directory and then merge them into a single table using MetaPhlAn2's merge_metaphlan_tables.py command. After this file is created we can fix the header so that each column corresponds to a sample name without the trailing "_metaphlan_bugs_list" description. Note that MetaPhlAn2 works best for well-characterized environments, like the human gut, and has low sensitivity in other environments.
    #os.system('mkdir metaphlan2_out')
    #os.system('cp humann2_out/*/*/*metaphlan_bugs_list.tsv metaphlan2_out/')
    #os.system('merge_metaphlan_tables.py metaphlan2_out/*metaphlan_bugs_list.tsv > metaphlan2_merged.txt')
    #os.system("sed -i 's/_metaphlan_bugs_list//g' metaphlan2_merged.txt")
    #3.7 Format STAMP taxonomy file
    #Lastly we can convert this MetaPhlAn2 abundance table to STAMP format
    #os.system('perl %s/metaphlan_to_stamp.pl metaphlan2_merged.txt > metaphlan2_merged.spf' % (script_path))



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='rd', type=str, required=True,
                        help="the tabular-table contains paths of the raw data")
    parser.add_argument('-n', '--node', dest='node', type=str, required=False, default='4',
                        help="the number of nodes to request")
    parser.add_argument('-sp', '--spath',dest='sp', type=str, required=False, default='/home/yxtan/HUMANN2_SOP_scripts/',
                        help="path of the custom scripts")
    parser.add_argument('-e', '--pair', dest='pair', type=str, required=False, default='True',
                        help="Is it pair-end seq data? Default is 'True'; Any other strings will be considered False")
    parser.add_argument('-j', '--jobs', dest='jobs', type=str, required=False, default='8',
                        help="The number of jobs run parallell in humann2 step. Default is '8'; It is bounded by the total number of memory available. Each job should have 16GB memory")
    args = parser.parse_args()
    print('Usage example with minimum parameters: python /home/yxtan/HUMANN2_SOP_scripts/Metagenomics_HUMANN2_only.py -i sample_table.txt -n 4')
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
    ID_list = []
    Counter=0
    File = open(rd_dir)
    for line in File:
        if line[0] == '#':
            continue
        line = line.strip()
        Fields = line.split("\t")
        print(Fields[0])
        if Counter==0:
            ID_list.append(Fields[0])
            F_list =Fields[1]
            R_list =Fields[2]
            ALL_list ='%s %s' % (Fields[1], Fields[2])
            Counter+=1
        else:
            if len(Fields) == 3:
                ID_list.append(Fields[0])
                F_list ='%s %s' % (F_list, Fields[1])
                R_list ='%s %s' % (R_list, Fields[2])
                ALL_list ='%s %s %s' % (ALL_list, Fields[1], Fields[2])
            elif len(Fields) == 1:
                if Fields[0] == '':
                    print("warning: there is an empty row, if more than 1 warning like this, please check the sample table")
                else:
                    print('The input sample table has less columns than it should be; Exit now.')
                    exit(0)
            else:
                print('The input sample table has more columns than it should be; Exit now.')
                exit(0)
    
    File.close()
    
    #1.2 Inspect read quality
    # if pair_end == 'True':
    #     m1_fastqc(ALL_list)
    # else:
    #     m1_fastqc(F_list)
    
    #2. Read Quality-Control and Contaminant Screens and connect to a long read
    # if pair_end == 'True':
    #     m2_run_kneaddata(nnodes, F_list,R_list,script_path)
    # else:
    #     m2_run_kneaddata_single(nnodes, F_list,script_path)
    
    #3. Determine Functions with HUMAnN2
    m3_humann2(nnodes, script_path,njobs)







