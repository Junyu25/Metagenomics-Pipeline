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

#This script will run QIIME2 automatically，mostly following the microbiome_helper workflow
#The output folder will be the current working folder (it is pointed by the workfd in qsub)
#Usage rule: python /home/yxtan/QIIME2_workflow.py -i raw_data_folder -n 20 -r reference_classification -fl len_of_forward_end -rl len_of_reverse_end -md metadata_path  -fp 'ACTCCTACGGGAGGCAGCA' -rp 'GGACTACHVGGGTWTCTAAT'

#This script starts from raw fqs before joined. The sample name should not contain R1 or R2, which will conflict with the R1\R2 annotation of 1st and 2nd ends. 
#Library file requirement:
#QIIME2 must be installed by anaconda
#cutadapt must be installed by anaconda

"""

#modules
def m1_fastqc(ALL_list, nnodes):
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

def m1a_cutadapt_pair(ALL_list,nnodes,script_path,fw_primer,rev_primer, F_list, R_list, file_list_out):
    """
    :rd_dir: the dir where the raw data are '/'.
    :nnodes: the number of nodes to parallel
    :script_path: the path of custome script, default is /home/yxtan/qiime2_custom_scripts/
    :return: None
    """    
    #1.3 Inspect read quality
    m1_fastqc(ALL_list, nnodes)
    #1.4 Trim primers with cutadapt
    os.system('mkdir primer_trimmed_fastqs')
    os.system('parallel -k --xapply -j %s \
    "cutadapt --pair-filter any \
    --no-indels \
    --discard-untrimmed \
    -g %s \
    -G %s \
    -o primer_trimmed_fastqs/{1/} \
    -p primer_trimmed_fastqs/{2/} \
    {1} {2} \
    > primer_trimmed_fastqs/{1/}_cutadapt_log.txt" \
    ::: %s ::: %s' % (nnodes, fw_primer, rev_primer, F_list, R_list))

    #1.4.2Inspect read quality after trimming
    if not os.path.isdir('primer_trimmed_fastqs/fastqc_out'):
        os.system('mkdir primer_trimmed_fastqs/fastqc_out')
        os.system('fastqc -t %s %s -o primer_trimmed_fastqs/fastqc_out/' % (nnodes, file_list_out))
    else:
        print('fastqc for primer_trimmed is skipped.')

    os.system('python %s/parse_cutadapt_logs.py -i primer_trimmed_fastqs/*log.txt' % (script_path))
    os.system('rm primer_trimmed_fastqs/*log.txt')

def m1a_cutadapt_single(nnodes,script_path,fw_primer, F_list, file_list_out_single):
    """
    :rd_dir: the dir where the raw data are '/'.
    :nnodes: the number of nodes to parallel
    :script_path: the path of custome script, default is /home/yxtan/qiime2_custom_scripts/
    :return: None
    """    
    #1.3 Inspect read quality
    m1_fastqc(F_list, nnodes)
    #1.4 Trim primers with cutadapt
    os.system('mkdir primer_trimmed_fastqs')
    os.system('parallel -k -j %s \
    "cutadapt --no-indels \
    --discard-untrimmed \
    -g %s \
    -o primer_trimmed_fastqs/{1/} \
    {1} \
    > primer_trimmed_fastqs/{1/}_cutadapt_log.txt" \
    ::: %s ' % (nnodes, fw_primer, F_list))
    #1.4.2Inspect read quality after trimming
    if not os.path.isdir('primer_trimmed_fastqs/fastqc_out'):
        os.system('mkdir primer_trimmed_fastqs/fastqc_out')
        os.system('fastqc -t %s %s -o primer_trimmed_fastqs/fastqc_out/' % (nnodes, file_list_out_single))
    else:
        print('fastqc for primer_trimmed is skipped.')

    os.system('python %s/parse_cutadapt_logs.py -i primer_trimmed_fastqs/*log.txt' % (script_path))
    os.system('rm primer_trimmed_fastqs/*log.txt')
    print(file_list_out_single)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='rd', type=str, required=True,
                        help="the tabular-table contains paths of the raw data")
    parser.add_argument('-n', '--node', dest='node', type=str, required=False, default='4',
                        help="the number of nodes to request")
    parser.add_argument('-sp', '--spath',dest='sp', type=str, required=False, default='/home/yxtan/qiime2_custom_scripts/',
                        help="path of the custom scripts")
    parser.add_argument('-fp', '--primerf', dest='fp', type=str, required=False, default='ACTCCTACGGGAGGCAGCA',
                        help="the seq of forward primer")
    parser.add_argument('-rp', '--primerr', dest='rp', type=str, required=False, default='GGACTACHVGGGTWTCTAAT',
                        help="the seq of reverse primer")
    parser.add_argument('-e', '--pair', dest='pair', type=str, required=False, default='True',
                        help="Is it pair-end seq data? Default is 'True'; Any other strings will be considered False")
    args = parser.parse_args()
    print('Usage example with minimum parameters: python /home/yxtan/qiime2_custom_scripts/QIIME2_fq_preprocess.py -i sample_table.txt')
    rd_dir = os.path.abspath(args.rd)
    script_path = os.path.abspath(args.sp)
    pair_end = args.pair
    nnodes = int(args.node)
    fw_primer = args.fp
    rev_primer = args.rp

    #########要增加standard out的check point。
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
            file1_list = Fields[1].split("/")
            file2_list = Fields[2].split("/")
            file_list_out ='primer_trimmed_fastqs/%s primer_trimmed_fastqs/%s' % (file1_list[len(file1_list)-1], file2_list[len(file2_list)-1])
            file_list_out_single ='primer_trimmed_fastqs/%s' % (file1_list[len(file1_list)-1])
            Counter+=1
        else:
            if len(Fields) == 3:
                ID_list.append(Fields[0])
                F_list ='%s %s' % (F_list, Fields[1])
                R_list ='%s %s' % (R_list, Fields[2])
                ALL_list ='%s %s %s' % (ALL_list, Fields[1], Fields[2])
                file1_list = Fields[1].split("/")
                file2_list = Fields[2].split("/")
                file_list_out ='%s primer_trimmed_fastqs/%s primer_trimmed_fastqs/%s' % (file_list_out, file1_list[len(file1_list)-1], file2_list[len(file2_list)-1])
                file_list_out_single ='%s primer_trimmed_fastqs/%s' % (file_list_out_single, file1_list[len(file1_list)-1])
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
    
    #run fastq and cutadapt
    if pair_end == 'True':
        m1a_cutadapt_pair(ALL_list,nnodes,script_path,fw_primer,rev_primer, F_list, R_list, file_list_out)
    else:
        m1a_cutadapt_single(nnodes,script_path,fw_primer, F_list, file_list_out_single)