#coding:utf-8
from __future__ import print_function
import glob,os,re,argparse

from pandas import DataFrame as df


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
#Usage rule: python /home/yxtan/qiime2_custom_scripts/QIIME2_workflow.py -i raw_data_folder -n 20 -r reference_classification -fl len_of_forward_end -rl len_of_reverse_end -md metadata_path  -fp 'ACTCCTACGGGAGGCAGCA' -rp 'GGACTACHVGGGTWTCTAAT'

#This script starts from raw fqs before joined. The sample name should not contain R1 or R2, which will conflict with the R1\R2 annotation of 1st and 2nd ends. 
#Library file requirement:
#QIIME2 must be installed by anaconda
#cutadapt must be installed by anaconda

"""

#modules
def import_manifest(pair_end, rd_dir):
    """
    :rd_dir: the path of sample table of the raw data '/'.
    :pair_end: define pair-end data as input or not
    :return: None
    """
    if pair_end == 'True':
        os.system('qiime tools import --type SampleData[PairedEndSequencesWithQuality] \
                       --input-path %s \
                       --output-path reads_qza/reads_trimmed.qza \
                       --input-format PairedEndFastqManifestPhred33V2' % (rd_dir))
    else:
        ####################single end
        os.system('qiime tools import --type SampleData[SequencesWithQuality] \
                       --input-path %s \
                       --output-path reads_qza/reads_trimmed.qza \
                       --input-format SingleEndFastqManifestPhred33V2' % (rd_dir))    

def m1b_import_data(pair_end,rd_dir):
    """
    :rd_dir: the path of sample table of the raw data '/'.
    :pair_end: define pair-end data as input or not
    :return: None
    """
    #1.5 Import primer trimmed FASTQs as QIIME2 artifact 最大的逻辑差异，主要就是统一了格式
    os.system('mkdir reads_qza')
    #如果是没有demultiplex的数据，那需要casava格式，实际上QIIME2在读入的时候会同时进行casava的分样品。。。。
    #如果是公司已经分开成具体每个样品了，那又是另外一种manifest格式。。。（这里就牵涉到要提供manifest，所以如果有给manifest（实际上可以上面cuteadapt那一步直接生成，所以只需要提供True or False），那就是这个格式，输入FALSE，那就是casava。
    import_manifest(pair_end, rd_dir)
    #1.6 Generating table of sample depths at each step
    #os.system('%s/qiime2_fastq_lengths.py [QZA1] [QZA2] ... [QZAN] --proc %s -o read_counts.tsv' % (script_path, nnodes)) 

def m2a_data_preprocess(m1_out, script_path, nnodes):
    """
    :m1_out: the out put from m1 step as input here'/'.
    :script_path: the path of custome script, default is /home/yxtan/qiime2_custom_scripts/
    :nnodes: the number of nodes to parallel
    :return: None
    """ 
    #2:比对分析，有两条路线A是deblur会快些（需要用vsearch），B是dada2貌似接受度高一些。
    #2A.1 Join paired-end reads
    if pair_end == 'True':
        os.system('qiime vsearch join-pairs --i-demultiplexed-seqs %s \
                     --o-joined-sequences reads_qza/reads_trimmed_joined.qza' % m1_out)
        #2A.2 Filter out low-quality reads.
        os.system('qiime quality-filter q-score-joined --i-demux reads_qza/reads_trimmed_joined.qza \
                    --o-filter-stats reads_qza/filt_stats.qza \
                    --o-filtered-sequences reads_qza/reads_trimmed_joined_filt.qza')
    else:
        #2A.2 Filter out low-quality reads.
        os.system('qiime quality-filter q-score --i-demux %s \
                    --o-filter-stats reads_qza/filt_stats.qza \
                    --o-filtered-sequences reads_qza/reads_trimmed_joined_filt.qza' % m1_out)
    ##########这里可能需要做个QC，检验一下read的损失高不高，高的话可能得预警
    if pair_end == 'True':
        os.system('python %s/qiime2_fastq_lengths.py %s reads_qza/reads_trimmed_joined.qza reads_qza/reads_trimmed_joined_filt.qza --proc %s -o reads_qza/read_counts.tsv' % (script_path , m1_out , nnodes)) 
    else:
        os.system('python %s/qiime2_fastq_lengths.py %s reads_qza/reads_trimmed_joined_filt.qza --proc %s -o reads_qza/read_counts.tsv' % (script_path , m1_out , nnodes)) 

def m2b_deblur(m2a_out, nnodes, trim_len_1):
    """
    :m2a_out: the out put from m2a step as input here'/'.
    :nnodes: the number of nodes to parallel
    :trim_len_1: the trim_len setting, default is 0 for auto detection
    :return: None
    """
    #2A.3 Running deblur
    #deblur有个假设是所有read都是同样长度，否则会报错，所以肯定是要用trim-len的，但是比这个短的read全部都会被扔掉，所以我觉得就用2%短的就好了，不然损失会太大。
    os.system('qiime demux summarize --i-data %s --o-visualization reads_trimmed_joined_filt_summary.qzv' % m2a_out)
    os.system('qiime tools export --input-path reads_trimmed_joined_filt_summary.qzv --output-path demux_sum/')
    trim_len_2 = int(os.popen("grep -A1 '<th>2%</th>' demux_sum/quality-plot.html | grep '<td>' | cut -d '>' -f2 | cut -f1 -d' '").readlines()[0].strip())
    if trim_len_1 < trim_len_2:
        trim_len = trim_len_2
        print('trim_length is %s, automatic detected' % trim_len_2)
    else:
        trim_len = trim_len_1
        print('trim_length is %s, by user input' % trim_len_1)
    os.system('qiime deblur denoise-16S --i-demultiplexed-seqs %s \
                     --p-trim-length %s \
                     --p-sample-stats \
                     --p-jobs-to-start %s \
                     --p-min-reads 1 \
                     --output-dir deblur_output' % (m2a_out,trim_len, nnodes))

def m2b_dada2(m1_out, nnodes, pair_end, len_R1, trunc_R1, trunc_R2, minOL):
    """
    :m1_out: the out put from m1 step as input here'/'.
    :nnodes: the number of nodes to parallel
    :pair_end: control whether the data is pair-end or single-end
    :trunc_R1:the length to be trunc in dada2 after cutadapt on the 1st end, only use it when you want to trim the low quality reads at the 3 end. The default is doing no trimming. Or you can use an extreme large number to activate using the full length and filter the short reads.
    :trunc_R2:the length to be trunc in dada2 after cutadapt on the 2st end, only use it when you want to trim the low quality reads at the 3 end. The default is doing no trimming. Or you can use an extreme large number to activate using the full length and filter the short reads.
    :minOL: the number of min overlap allowed for dada2
    :return: None
    """
    #2B. Running DADA2 workflow (Option 2):
    #2B.1 Running DADA2
    #############dada2对单端的是不是又不一样？
    max_ee = len_R1//100+1
    if max_ee < 2:
        max_ee = 2
    #其实trunc-len可以用0，就不会trime了，但是那样的话，结果可能会长度参差不齐，比较好的方式就是用理论长度或者cuteadapt以后的fastqc的结果
    if os.path.isdir('dada2_output'):
        os.system('rm -rf dada2_output')
    if pair_end == 'True':
        print('min-overlap is %s' % minOL)
        os.system('qiime dada2 denoise-paired --i-demultiplexed-seqs %s \
                               --p-trunc-len-f %s \
                               --p-trunc-len-r %s \
                               --p-max-ee %s \
                               --p-n-threads %s \
                               --p-min-overlap %s \
                               --output-dir dada2_output' % (m1_out, trunc_R1, trunc_R2, max_ee, nnodes,minOL))
    else:
        os.system('qiime dada2 denoise-single --i-demultiplexed-seqs %s \
                       --p-trunc-len %s \
                       --p-max-ee %s \
                       --p-n-threads %s \
                       --output-dir dada2_output' % (m1_out, trunc_R1, max_ee, nnodes))
    #2B.2 Convert stats QZA to TSV
    #You should take a look at the read count table to see how many reads were retained at each step of the DADA2 pipeline:
    os.system('qiime tools export --input-path dada2_output/denoising_stats.qza --output-path dada2_output')

#now no filter will be used in the pipeline. Filters will be used in the downstream analysis.
def m3_filter_ASVs(m2b_out,Freq):    
    """
    :m2b_out: the out put from m2b step as input here'/'.
    :Freq:default is -1, then it will be auto detected by the Mean frequency of sample and divided by 1000. User can input, but features smaller than it will be filtered.
    :return: None
    """
    #后面的每一个大模块都同时对两个方法的结果适用，因此都应该模块化
    #3.1 Summarize output table
    os.system('qiime feature-table summarize --i-table %s/table.qza --o-visualization %s/table_summary.qzv' % (m2b_out, m2b_out))
    os.system('qiime tools export --input-path %s/table_summary.qzv --output-path %s/table_sum/' % (m2b_out, m2b_out))
    
    #to calculate percentage，need to tansfer to relative abundance by https://docs.qiime2.org/2019.7/plugins/available/feature-table/relative-frequency/
    os.system('qiime feature-table relative-frequency --i-table %s/table.qza --o-relative-frequency-table %s/table_relative.qza' % (m2b_out, m2b_out))
    
    #just cp, so the two file will be redundance, but just for the convinient of modification of downstream names.
    os.system('cp %s/table.qza %s/table_filt.qza' % (m2b_out, m2b_out))
    os.system('cp %s/representative_sequences.qza %s/rep_seqs_filt.qza' % (m2b_out, m2b_out))
    
    #Finally, you can make a new summary of the filtered abundance table:
    os.system('qiime feature-table summarize --i-table %s/table_filt.qza --o-visualization %s/table_filt_summary.qzv' % (m2b_out, m2b_out))
    os.system('qiime tools export --input-path %s/table_filt_summary.qzv --output-path %s/table_filt_sum/' % (m2b_out, m2b_out))

#this function will not be used in this pipeline
def m3_filter_ASVs_MH(m2b_out,Freq):    
    """
    :m2b_out: the out put from m2b step as input here'/'.
    :Freq:default is -1, then it will be auto detected by the Mean frequency of sample and divided by 1000. User can input, but features smaller than it will be filtered.
    :return: None
    """
    #后面的每一个大模块都同时对两个方法的结果适用，因此都应该模块化
    #3.1 Summarize output table
    os.system('qiime feature-table summarize --i-table %s/table.qza --o-visualization %s/table_summary.qzv' % (m2b_out, m2b_out))
    os.system('qiime tools export --input-path %s/table_summary.qzv --output-path %s/table_sum/' % (m2b_out, m2b_out))
    #from data/index.html, extract Mean frequency, divided by 1000, and round, get Freq
    freq_0001 = int(float(os.popen("grep -A1 'Mean frequency' %s/table_sum/index.html | grep '<td>' | cut -d '>' -f2 | cut -f1 -d'<'" % m2b_out).readlines()[0].strip().replace(',',''))//1000)
    if Freq <0:
        Freq = freq_0001
        print('The auto detect frequency is %s' % Freq)
    #3.2 Filter out rare ASVs
    ############Based on the above summary visualization you can choose a cut-off for how frequent a variant needs to be (and optionally how many samples need to have the variant) for it to be retained. One possible choice would be to remove all ASVs that have a frequency of less than 0.1% of the mean sample depth. This cut-off excludes ASVs that are likely due to MiSeq bleed-through between runs (reported by Illumina to be 0.1% of reads). To calculate this cut-off you would identify the mean sample depth in the above visualization, multiply it by 0.001, and round to the nearest integer.（这一步其实应该可以自动化，提取计算值而已）
    os.system('qiime feature-table filter-features --i-table %s/table.qza \
                                        --p-min-frequency %s \
                                        --p-min-samples 1 \
                                        --o-filtered-table %s/table_filt_mh.qza' % (m2b_out, Freq, m2b_out))
    #Since the ASVs will be in a separate FASTA file you can exclude the low frequency ASVs from the sequence file with this command:
    os.system('qiime feature-table filter-seqs --i-data %s/representative_sequences.qza \
                                    --i-table %s/table_filt_mh.qza \
                                    --o-filtered-data %s/rep_seqs_filt_mh.qza' % (m2b_out, m2b_out, m2b_out))
    #这里目前缺少了一个从relative去filter total的方式，可能需要创建一个新的plugin，基于relative的命令集
    
    #Finally, you can make a new summary of the filtered abundance table:
    os.system('qiime feature-table summarize --i-table %s/table_filt_mh.qza --o-visualization %s/table_filt_mh_summary.qzv' % (m2b_out, m2b_out))
    os.system('qiime tools export --input-path %s/table_filt_mh_summary.qzv --output-path %s/table_filt_mh_sum/' % (m2b_out, m2b_out))

def m4_phylogeny(m3_out,m2b_out):    
    """
    :m2a_out: the out put from m2a step as input here'/'.
    :m3_out: the out put from m3 step as input here'/'.
    :return: None
    """
    #4. Build quick phylogeny with FastTree（图像化）
    #4.1 Making multiple-sequence alignment
    os.system('mkdir %s/tree_out' % m2b_out)
    os.system('qiime alignment mafft --i-sequences %s \
                          --p-n-threads %s \
                          --o-alignment %s/tree_out/rep_seqs_filt_aligned.qza' % (m3_out, nnodes, m2b_out))
    #4.2 Filtering multiple-sequence alignment
    os.system('qiime alignment mask --i-alignment %s/tree_out/rep_seqs_filt_aligned.qza \
    --o-masked-alignment %s/tree_out/rep_seqs_filt_aligned_masked.qza' % (m2b_out, m2b_out))
    #4.3 Running FastTree
    os.system('qiime phylogeny fasttree --i-alignment %s/tree_out/rep_seqs_filt_aligned_masked.qza \
                             --p-n-threads %s \
                             --o-tree %s/tree_out/rep_seqs_filt_aligned_masked_tree' % (m2b_out, nnodes, m2b_out))
    #4.4 Add root to tree
    os.system('qiime phylogeny midpoint-root --i-tree %s/tree_out/rep_seqs_filt_aligned_masked_tree.qza \
                                  --o-rooted-tree %s/tree_out/rep_seqs_filt_aligned_masked_tree_rooted.qza' % (m2b_out, m2b_out))

def m5_rarefaction(m3_table_filt_out,m4_root_tree_out,m2b_out,max_depth,metadata_category):    
    """
    :m3_table_filt_out: the table filtered output from m3 step as input here'/'.
    :m4_root_tree_out: the rooted tree output from m1 step as input here'/'.
    :m2b_out: the folder path output from m2 step as input here'/'.
    :max_depth: default is 0 and auto detection will be done. User can specify a input, however, if the input is bigger than the max_depth, the program will collapse.
    :metadata_category: only the category columns in the metadata will be process.
    :return: None
    """
    #extract max_depth first
    max_depth_auto = int(float(os.popen("grep -A1 'Minimum frequency' %s/table_filt_sum/index.html | grep '<td>' | cut -d '>' -f2 | cut -f1 -d'<'" % m2b_out).readlines()[0].strip().replace(',','')))
    if max_depth <=0:
        max_depth = max_depth_auto
        print('The auto detect max_depth for rarefaction (in fact, the min depth across samples) is %s' % max_depth_auto)
    #5. Generate rarefaction curves
    #A key quality control step is to plot rarefaction curves for all of your samples to determine if you performed sufficient sequencing. 
    #grouped plot
    os.system('qiime diversity alpha-rarefaction --i-table %s \
                                      --p-max-depth %s \
                                      --p-steps 20 \
                                      --i-phylogeny %s \
                                      --m-metadata-file %s \
                                      --o-visualization %s/rarefaction_curves.qzv' % (m3_table_filt_out, max_depth, m4_root_tree_out, metadata_category, m2b_out))
    #each sample plot（without using metadata info）
    os.system('qiime diversity alpha-rarefaction --i-table %s \
                                      --p-max-depth %s \
                                      --p-steps 20 \
                                      --i-phylogeny %s \
                                      --o-visualization %s/rarefaction_curves_eachsample.qzv' % (m3_table_filt_out, max_depth, m4_root_tree_out, m2b_out))

def m6_taxonomy(m3_out,m2b_out,classified_ref,nnodes,m3_table_filt_out,metadata_category):    
    """
    :m3_out: the out put from m3 step as input here'/'.
    :m3_table_filt_out: the table filtered output from m3 step as input here'/'.
    :classified_ref: the reference based classification model as primery input'/'.
    :m2b_out: the folder path output from m2 step as input here'/'.
    :nnodes: the number of nodes to parallel
    :metadata_category: only the category columns in the metadata will be process
    :return: None
    """
    #6. Assign taxonomy
    #You can assign taxonomy to your ASVs using a Naive-Bayes approach implemented in the scikit learn Python library and the SILVA database. 
    if os.path.isdir('%s/diversity' % m2b_out):
        os.system('rm -rf %s/taxa' % m2b_out)
    os.system('qiime feature-classifier classify-sklearn --i-reads %s \
                                              --i-classifier %s \
                                              --p-n-jobs %s \
                                              --output-dir %s/taxa ' % (m3_out, classified_ref, nnodes, m2b_out))
    os.system('qiime tools export --input-path %s/taxa/classification.qza --output-path %s/taxa' % (m2b_out, m2b_out))
    #the interactive stacked bar-charts of the taxonomic abundances across samples
    os.system('qiime taxa barplot --i-table %s \
                       --i-taxonomy %s/taxa/classification.qza \
                       --m-metadata-file %s \
                       --o-visualization %s/taxa/taxa_barplot.qzv' % (m3_table_filt_out, m2b_out, metadata_category, m2b_out))
    meta_f = open(metadata_category)
    L1=meta_f.readline()
    col_ID=L1.rstrip().split('\t')[1:]
    meta_f.close()
    #For each category，run grouped analysis
    for i_ID in col_ID:
        os.system('qiime feature-table group --i-table %s \
                              --p-axis sample \
                              --p-mode sum \
                              --m-metadata-file %s \
                              --m-metadata-column %s \
                              --o-grouped-table %s/feature_table_filt_%s.qza' % (m3_table_filt_out, metadata_category, i_ID, m2b_out, i_ID))
        os.system('qiime taxa barplot --i-table %s/feature_table_filt_%s.qza \
                       --i-taxonomy %s/taxa/classification.qza \
                       --m-metadata-file %s \
                       --o-visualization %s/taxa/taxa_barplot_%s.qzv' % (m2b_out, i_ID, m2b_out, metadata_category, m2b_out, i_ID))

def m7_diversity(m3_table_filt_out,m4_root_tree_out, nnodes,m2b_out,metadata_category,min_depth):    
    """
    :m3_table_filt_out: the table filtered output from m3 step as input here'/'.
    :m4_root_tree_out: the rooted tree output from m1 step as input here'/'.
    :m2b_out: the folder path output from m2 step as input here'/'.
    :nnodes: the number of nodes to parallel.
    :metadata_category: only the category columns in the metadata will be process
    :min_depth: default is 0 for auto detection; otherwise, user can input, but samples with depth below this cut-off will be excluded.
    :return: None
    """
    min_depth_auto = int(float(os.popen("grep -A1 'Minimum frequency' %s/table_filt_sum/index.html | grep '<td>' | cut -d '>' -f2 | cut -f1 -d'<'" % m2b_out).readlines()[0].strip().replace(',','')))
    if min_depth <=0:
        min_depth = min_depth_auto
        print('The auto detect min_depth_auto for diversity calculation is %s' % min_depth_auto)
    #7. Calculating diversity metrics and generating ordination plots
    if os.path.isdir('%s/diversity' % m2b_out):
        os.system('rm -rf %s/diversity' % m2b_out)
    os.system('qiime diversity core-metrics-phylogenetic --i-table %s \
                                              --i-phylogeny %s \
                                              --p-sampling-depth %s \
                                              --m-metadata-file %s \
                                              --p-n-jobs %s \
                                              --output-dir %s/diversity' % (m3_table_filt_out, m4_root_tree_out, min_depth, metadata_category, nnodes, m2b_out))
    os.system('qiime diversity alpha-group-significance --i-alpha-diversity %s/diversity/shannon_vector.qza \
                                             --m-metadata-file %s \
                                             --o-visualization %s/diversity/shannon_compare_groups.qzv' % (m2b_out, metadata_category, m2b_out))
    os.system('qiime diversity alpha-group-significance \
  --i-alpha-diversity %s/diversity/faith_pd_vector.qza \
  --m-metadata-file %s \
  --o-visualization %s/faith-pd-group-significance.qzv' % (m2b_out, metadata_category, m2b_out))
    os.system('qiime diversity alpha-group-significance \
  --i-alpha-diversity %s/diversity/evenness_vector.qza \
  --m-metadata-file %s \
  --o-visualization %s/evenness-group-significance.qzv' % (m2b_out, metadata_category, m2b_out))

def m8_ANCOM(m3_table_filt_out, m2b_out,metadata_category):    
    """
    :m3_table_filt_out: the table filtered output from m3 step as input here'/'.
    :m2b_out: the folder path output from m2 step as input here'/'.
    :metadata_category: only the category columns in the metadata will be process.
    :return: None
    """
    #8. Identifying deferentially abundant features with ANCOM
    #ANCOM is one method to test for difference in the relative abundance of features between sample groupings. It is a compositional approach that makes no assumptions about feature distributions. However, it requires that all features have non-zero abundances so a pseudocount of 1 first needs to be added:
    os.system('qiime composition add-pseudocount --i-table %s \
                                      --o-composition-table %s/table_filt_pseudocount.qza' % (m3_table_filt_out, m2b_out))
    #############note that CATEGORY is a placeholder for the text label of your category of interest from the metadata file:
    meta_f = open(metadata_category)
    L1=meta_f.readline()
    col_ID=L1.rstrip().split('\t')[1:]
    meta_f.close()
    #For each category，run grouped analysis
    for i_ID in col_ID:
        if os.path.isdir('%s/ancom_%s_output' % (m2b_out, i_ID)):
            os.system('rm -rf %s/ancom_%s_output' % (m2b_out, i_ID))
        os.system('qiime composition ancom --i-table %s/table_filt_pseudocount.qza \
                            --m-metadata-file %s \
                            --m-metadata-column %s \
                            --output-dir %s/ancom_%s_output' % (m2b_out, metadata_category, i_ID, m2b_out, i_ID))

def m9_exporting(m3_out,m2b_out,m3_table_filt_out, script_path, m4_out):    
    """
    :m3_out: the out put from m3 step as input here'/'.
    :m3_table_filt_out: the table filtered output from m3 step as input here'/'.
    :m2b_out: the folder path output from m2 step as input here'/'.
    :nnodes: the number of nodes to parallel
    :script_path: the path of custome script, default is /home/yxtan/qiime2_custom_scripts/
    :return: None
    """
    #9. Exporting the final abundance, profile and sequence files
    #Lastly, to get the BIOM file (with associated taxonomy), the STAMP profile file (.spf) and FASTA file (one per ASV) for your dataset to plug into other programs you can use the commands below. Note that there are a couple of file manipulations included here since the conversion scripts expect certain slightly different headers than those produced by QIIME2 by default and certain SILVA taxonomy labeling problems need to be fixed for STAMP:
    # For FASTA:
    os.system('qiime tools export --input-path %s --output-path %s/output_exported' % (m3_out, m2b_out))
    # For BIOM w/taxonomy:
    os.system("sed -i -e '1 s/Feature/#Feature/' -e '1 s/Taxon/taxonomy/' %s/taxa/taxonomy.tsv" % (m2b_out))
    if os.path.isdir('%s/output_exported' % (m2b_out)):
        os.system('rm -rf %s/output_exported' % (m2b_out))
    os.system('qiime tools export --input-path %s --output-path %s/output_exported' % (m3_table_filt_out, m2b_out))
    os.system('biom add-metadata -i %s/output_exported/feature-table.biom -o %s/output_exported/feature-table_w_tax.biom --observation-metadata-fp %s/taxa/taxonomy.tsv --sc-separated taxonomy' % (m2b_out, m2b_out, m2b_out))
    os.system('biom convert -i %s/output_exported/feature-table_w_tax.biom -o %s/output_exported/feature-table_w_tax.txt --to-tsv --header-key taxonomy' % (m2b_out, m2b_out))
    #For STAMP profile:
    os.system('python %s/biom_to_stamp.py -m taxonomy %s/output_exported/feature-table_w_tax.biom > %s/output_exported/feature-table_w_tax.spf' % (script_path , m2b_out , m2b_out))
    os.system('python %s/fix_spf.py -i %s/output_exported/feature-table_w_tax.spf -o %s/output_exported/feature-table_w_tax_final.spf' % (script_path , m2b_out , m2b_out))
    #export tree.nwk
    os.system('qiime tools export --input-path %s --output-path %s/output_exported/rooted_tree' % (m4_out, m2b_out))

def m10_picrust2(m2b_out,tax_tb,rep_seq,nnodes, metadata_category,script_path):
    """
    :m2b_out: the out put from m2 step as path here'/'.
    :tax_tb: the ASV count table, such as table.qza
    :rep_seq: the ASV seq file, such as representative_sequences.qza
    :nnodes: the number of nodes to parallel
    :metadata_category: only the category columns in the metadata will be process.
    :script_path: the path of custome script, default is /home/yxtan/qiime2_custom_scripts/
    :return: None
    """
    #10. do picrust2 analysis and export outputs
    os.system('qiime picrust2 full-pipeline \
       --i-table %s \
       --i-seq %s \
       --output-dir %s/picrust2_output \
       --p-threads %s \
       --verbose' % (tax_tb, rep_seq, m2b_out,nnodes))
    m10_picrust2_sub(m2b_out,"ec_metagenome", nnodes, metadata_category,script_path)
    m10_picrust2_sub(m2b_out,"ko_metagenome", nnodes, metadata_category,script_path)
    m10_picrust2_sub(m2b_out,"pathway_abundance", nnodes, metadata_category,script_path)


def m10_picrust2_sub(m2b_out,out_type, nnodes, metadata_category,script_path):
    """
    :m1_out: the out put from m1 step as input here'/'.
    :nnodes: the number of nodes to parallel
    :return: None
    """
    os.system('qiime feature-table summarize \
       --i-table %s/picrust2_output/%s.qza \
       --o-visualization %s/picrust2_output/%s.qzv'% (m2b_out,out_type,m2b_out,out_type))
    #get min depth and generate PCoA plot
    os.system('qiime tools export --input-path %s/picrust2_output/%s.qzv --output-path %s/picrust2_output/%s/' % (m2b_out,out_type,m2b_out,out_type))
    min_freq = int(float(os.popen("grep -A1 'Minimum frequency' %s/picrust2_output/%s/index.html | grep '<td>' | cut -d '>' -f2 | cut -f1 -d'<'" % (m2b_out,out_type)).readlines()[0].strip().replace(',','')))   
    print('The auto detect min_freq for %s calculation is %s' % (out_type,min_freq))
    os.system('qiime diversity core-metrics \
       --i-table %s/picrust2_output/%s.qza \
       --p-sampling-depth %s \
       --m-metadata-file %s \
       --output-dir %s/picrust2_output/%s_diversity_out \
       --p-n-jobs %s' % (m2b_out,out_type, min_freq, metadata_category, m2b_out, out_type,nnodes))
    
    #transfer outputs into more readable formats
    os.system('qiime tools export --input-path %s/picrust2_output/%s.qza --output-path %s/picrust2_output/%s/' % (m2b_out,out_type,m2b_out,out_type))
    os.system('biom convert \
       -i %s/picrust2_output/%s/feature-table.biom \
       -o %s/picrust2_output/%s/feature-table.biom.tsv \
       --to-tsv'% (m2b_out,out_type,m2b_out,out_type))
    #for pathway analysis, there is no taxa info to add in
    #os.system('biom add-metadata -i %s/picrust2_output/%s/feature-table.biom -o %s/picrust2_output/%s/feature-table_w_tax.biom --observation-metadata-fp %s/taxa/taxonomy.tsv --sc-separated taxonomy' % (m2b_out,out_type, m2b_out, out_type,m2b_out))
    #os.system('biom convert -i %s/picrust2_output/%s/feature-table_w_tax.biom -o %s/picrust2_output/%s/feature-table_w_tax.txt --to-tsv --header-key taxonomy' % (m2b_out,out_type, m2b_out,out_type))
    #For STAMP profile, but it does not work because of lack of metainfo:
    #os.system('python %s/biom_to_stamp.py -m taxonomy %s/picrust2_output/%s/feature-table_w_tax.biom > %s/picrust2_output/%s/feature-table_w_tax.spf' % (script_path , m2b_out ,out_type, m2b_out,out_type))
    #os.system('python %s/fix_spf.py -i %s/picrust2_output/%s/feature-table_w_tax.spf -o %s/picrust2_output/%s/feature-table_w_tax_final.spf' % (script_path , m2b_out ,out_type, m2b_out,out_type))



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='rd', type=str, required=True,
                        help="the tabular-table contains paths of the raw data")
    parser.add_argument('-n', '--node', dest='node', type=str, required=False, default='4',
                        help="the number of nodes to request")
    parser.add_argument('-r', '--ref', dest='rf', type=str, required=True,
                        help="paht of the reference based classification model file ")
    parser.add_argument('-sp', '--spath',dest='sp', type=str, required=False, default='/home/yxtan/qiime2_custom_scripts/',
                        help="path of the custom scripts")
    parser.add_argument('-fl', '--forwardlen', dest='fl', type=int, required=True,
                        help="the length of forward primer")
    parser.add_argument('-rl', '--reverselen', dest='rl', type=int, required=True,
                        help="the length of reverse primer")
    parser.add_argument('-md', '--metadata', dest='md', type=str, required=True,
                        help="the path of the metadata file, it must contain only the categorical columns, or the plotting part of the workflow will give alarm")
    parser.add_argument('-fr', '--freq', dest='fr', type=str, required=False, default=-1,
                        help="the frequence used to do rare ASV filter, default is using auto-dection of 0.001 from mean frequency")
    parser.add_argument('-ad', '--deptha', dest='ad', type=str, required=False, default=0,
                        help="the depth used for rarefaction, generally should used the maximum depth from the summary plot, which by default will be auto detected. User can manually input the desired number bigger than 0.")
    parser.add_argument('-id', '--depthi', dest='id', type=str, required=False, default=0,
                        help="the depth used for PCoA, generally should used the minimum depth from the summary plot, which by default will be auto detected. User can manually input the desired number bigger than 0.However, samples with depth below this cut-off will be excluded")
    parser.add_argument('-fp', '--primerf', dest='fp', type=str, required=False, default='ACTCCTACGGGAGGCAGCA',
                        help="the seq of forward primer")
    parser.add_argument('-rp', '--primerr', dest='rp', type=str, required=False, default='GGACTACHVGGGTWTCTAAT',
                        help="the seq of reverse primer")
    parser.add_argument('-a', '--aligner', dest='a', type=str, required=False, default="both",
                        help="the user can choose what aligners to use, default is using both deblur and dada2. Otherwise,you can type in only 'deblur' or 'dada2' as option. Other typeins will all be consider both")
    parser.add_argument('-tl', '--trimelen', dest='tl', type=int, required=False, default=0,
                        help="the user defined trim length for deblur, which must be bigger than the auto detection. The default is using auto-detect 2 percent length after quality filtering.")
    parser.add_argument('-tr1', '--truncR1', dest='tr1', type=int, required=False, default=0,
                        help="the length to be trunc in dada2 after cutadapt on the 1st end, only use it when you want to trim the low quality reads at the 3 end. The default is doing no trimming. Or you can use an extreme large number to activate using the full length and filter the short reads.")
    parser.add_argument('-tr2', '--truncR2', dest='tr2', type=int, required=False, default=0,
                        help="the length to be trunc in dada2 after cutadapt on the 2nd end, only use it when you want to trim the low quality reads at the 3 end. The default is doing no trimming. Or you can use an extreme large number to activate using the full length and filter the short reads.")
    parser.add_argument('-e', '--pair', dest='pair', type=str, required=False, default='True',
                        help="Is it pair-end seq data? Default is 'True'; Any other strings will be considered False")
    parser.add_argument('-mo', '--minOL', dest='minO', type=str, required=False, default='12',
                        help="The min overlap allowed for dada2, default is 12")
    args = parser.parse_args()
    print('Usage example with minimum parameters: python /home/yxtan/qiime2_custom_scripts/QIIME2_workflow.py -i sample_table.txt  -fl 249 -rl 237 -r /home/yxtan/SILVA_132_SSURef_Nr99/SILVA_132_SSURef_Nr99_tax_silva_B337af_B787cr.qza -md /home/yxtan/test2/metadata_categorical.txt')
    rd_dir = os.path.abspath(args.rd)
    script_path = os.path.abspath(args.sp)
    pair_end = args.pair
    nnodes = int(args.node)
    classified_ref=os.path.abspath(args.rf)
    fw_primer = args.fp
    rev_primer = args.rp
    f_len = args.fl
    r_len = args.rl
    trunc_R1 = args.tr1
    trunc_R2 = args.tr2
    trim_len = args.tl
    fp_len = len(fw_primer)
    rp_len = len(rev_primer)
    trunc_f_len = f_len - fp_len
    trunc_r_len = r_len - rp_len
    Freq = args.fr
    max_depth=args.ad
    min_depth=args.id
    metadata_category=os.path.abspath(args.md)
    aligner = args.a
    minOL = int(args.minO)
    #########要增加standard out的check point。
    #需要检查输入的参数是否正确，主要是路径是否存在    
    if not os.path.isfile(rd_dir):
        print('Input sample table is not exist; Exit now.')
        exit(0)
    if not os.path.isfile(classified_ref):
        print('The ref based classification model is not exist; Exit now.')
        exit(0)
    if not os.path.isfile(metadata_category):
        print('The ref based metadata file is not exist; Exit now.')
        exit(0)
    if not os.path.isdir(script_path):
        print('The folder of custom scripts is not exist; Exit now.')
        exit(0)
    
    ##task
    #m1b
    m1b_import_data(pair_end,rd_dir)
    
    ##task
    #m2a
    ###input
    m2a_out="reads_qza/reads_trimmed_joined_filt.qza"
    m1_out="reads_qza/reads_trimmed.qza"
    #calculate trunc_R1和trunc_R2
    if trunc_R1 ==0:
        trunc_R1 = 0
        print('No trimming for dada2.')
    elif trunc_R1 > trunc_f_len:
        trunc_R1 = trunc_f_len
        print('All first end were trimmed by user input length: %s.' % trunc_R1)
    else:
        trunc_f_len = trunc_R1
        print('All first end were trimmed by user input length: %s.' % trunc_R1)
    
    if trunc_R2 ==0:
        trunc_R2 = 0
        print('No trimming for dada2.')
    elif trunc_R2 > trunc_r_len:
        trunc_R2 = trunc_r_len
        print('All second end were trimmed by user input length: %s.' % trunc_R2)
    else:
        trunc_r_len = trunc_R2
        print('All second end were trimmed by user input length: %s.' % trunc_R2)
    
    if aligner == "deblur":
        #用dada2的话，其实根本不需要跑m2a这一步，但是有一些统计结果，所以m2a好像也没有省略掉的必要,但是对于大文件来说，这一步实在是太慢了，主要是解压和压缩。。所以能省则省
        m2a_data_preprocess(m1_out, script_path, nnodes)
        m2b_deblur(m2a_out, nnodes, trim_len)
    elif aligner == "dada2":
        m2b_dada2(m1_out, nnodes, pair_end, trunc_f_len, trunc_R1, trunc_R2, minOL)
    else:
        #用dada2的话，其实根本不需要跑m2a这一步，但是有一些统计结果，所以m2a好像也没有省略掉的必要,但是对于大文件来说，这一步实在是太慢了，主要是解压和压缩。。所以能省则省
        m2a_data_preprocess(m1_out, script_path, nnodes)
        m2b_deblur(m2a_out, nnodes, trim_len)
        m2b_dada2(m1_out, nnodes, pair_end, trunc_f_len, trunc_R1, trunc_R2, minOL)
    
    ##task
    #m3
    ###input
    m2b_deblur_out="deblur_output"
    m2b_dada2_out="dada2_output"
    if aligner == "deblur":
        m3_filter_ASVs(m2b_deblur_out,Freq)
    elif aligner == "dada2":
        m3_filter_ASVs(m2b_dada2_out,Freq)
    else:
        m3_filter_ASVs(m2b_deblur_out,Freq)
        m3_filter_ASVs(m2b_dada2_out,Freq)
    
    ##task
    #m4
    ###input
    m3_out_deblur = '%s/rep_seqs_filt.qza' % m2b_deblur_out
    m3_out_dada2 = '%s/rep_seqs_filt.qza' % m2b_dada2_out
    if aligner == "deblur":
        m4_phylogeny(m3_out_deblur,m2b_deblur_out)
    elif aligner == "dada2":
        m4_phylogeny(m3_out_dada2,m2b_dada2_out)
    else:
        m4_phylogeny(m3_out_deblur,m2b_deblur_out)
        m4_phylogeny(m3_out_dada2,m2b_dada2_out)
    
    print('step 4 finished')
    ##task
    #m5 to m9 are all parallel, no dependency among them
    ###input
    m3_table_filt_deblur = '%s/table_filt.qza' % m2b_deblur_out
    m3_table_filt_dada2 = '%s/table_filt.qza' % m2b_dada2_out
    m4_root_tree_deblur = '%s/tree_out/rep_seqs_filt_aligned_masked_tree_rooted.qza' % m2b_deblur_out
    m4_root_tree_dada2 = '%s/tree_out/rep_seqs_filt_aligned_masked_tree_rooted.qza' % m2b_dada2_out
    #######这里一定要用category，我感觉metadata其实基本上都是用category，时间序列那行感觉好像用不上。
    if aligner == "deblur":
        m5_rarefaction(m3_table_filt_deblur,m4_root_tree_deblur,m2b_deblur_out,max_depth,metadata_category)
    elif aligner == "dada2":
        m5_rarefaction(m3_table_filt_dada2,m4_root_tree_dada2,m2b_dada2_out,max_depth,metadata_category)
    else:
        m5_rarefaction(m3_table_filt_deblur,m4_root_tree_deblur,m2b_deblur_out,max_depth,metadata_category)
        m5_rarefaction(m3_table_filt_dada2,m4_root_tree_dada2,m2b_dada2_out,max_depth,metadata_category)
    
    ##task
    #m6 
    ###input
    m3_out_deblur = '%s/rep_seqs_filt.qza' % m2b_deblur_out
    m3_out_dada2 = '%s/rep_seqs_filt.qza' % m2b_dada2_out
    m3_table_filt_deblur = '%s/table_filt.qza' % m2b_deblur_out
    m3_table_filt_dada2 = '%s/table_filt.qza' % m2b_dada2_out
    if aligner == "deblur":
        m6_taxonomy(m3_out_deblur,m2b_deblur_out,classified_ref,nnodes,m3_table_filt_deblur,metadata_category)
    elif aligner == "dada2":
        m6_taxonomy(m3_out_dada2,m2b_dada2_out,classified_ref,nnodes,m3_table_filt_dada2,metadata_category)
    else:
        m6_taxonomy(m3_out_deblur,m2b_deblur_out,classified_ref,nnodes,m3_table_filt_deblur,metadata_category)
        m6_taxonomy(m3_out_dada2,m2b_dada2_out,classified_ref,nnodes,m3_table_filt_dada2,metadata_category)
    
    ##task
    #m7 Calculating diversity metrics and generating ordination plots
    ###input
    m3_table_filt_deblur = '%s/table_filt.qza' % m2b_deblur_out
    m3_table_filt_dada2 = '%s/table_filt.qza' % m2b_dada2_out
    m4_root_tree_deblur = '%s/tree_out/rep_seqs_filt_aligned_masked_tree_rooted.qza' % m2b_deblur_out
    m4_root_tree_dada2 = '%s/tree_out/rep_seqs_filt_aligned_masked_tree_rooted.qza' % m2b_dada2_out
    if aligner == "deblur":
        m7_diversity(m3_table_filt_deblur,m4_root_tree_deblur, nnodes,m2b_deblur_out,metadata_category,min_depth)
    elif aligner == "dada2":
        m7_diversity(m3_table_filt_dada2,m4_root_tree_dada2, nnodes,m2b_dada2_out,metadata_category,min_depth)
    else:
        m7_diversity(m3_table_filt_deblur,m4_root_tree_deblur, nnodes,m2b_deblur_out,metadata_category,min_depth)
        m7_diversity(m3_table_filt_dada2,m4_root_tree_dada2, nnodes,m2b_dada2_out,metadata_category,min_depth)
    
    ##task
    #m8 Identifying deferentially abundant features with ANCOM
    ###input
    m3_table_filt_deblur = '%s/table_filt.qza' % m2b_deblur_out
    m3_table_filt_dada2 = '%s/table_filt.qza' % m2b_dada2_out
    if aligner == "deblur":
        m8_ANCOM(m3_table_filt_deblur, m2b_deblur_out,metadata_category)
    elif aligner == "dada2":
        m8_ANCOM(m3_table_filt_dada2, m2b_dada2_out,metadata_category)
    else:
        m8_ANCOM(m3_table_filt_deblur, m2b_deblur_out,metadata_category)
        m8_ANCOM(m3_table_filt_dada2, m2b_dada2_out,metadata_category)
    
    ##task
    #m9 Exporting the final abundance, profile and sequence files
    ###input
    m3_out_deblur = '%s/rep_seqs_filt.qza' % m2b_deblur_out
    m3_out_dada2 = '%s/rep_seqs_filt.qza' % m2b_dada2_out
    m3_table_filt_deblur = '%s/table_filt.qza' % m2b_deblur_out
    m3_table_filt_dada2 = '%s/table_filt.qza' % m2b_dada2_out
    m4_out_deblur = '%s/tree_out/rep_seqs_filt_aligned_masked_tree_rooted.qza' % m2b_deblur_out
    m4_out_dada2 = '%s/tree_out/rep_seqs_filt_aligned_masked_tree_rooted.qza' % m2b_dada2_out
    if aligner == "deblur":
        m9_exporting(m3_out_deblur,m2b_deblur_out,m3_table_filt_deblur,script_path, m4_out_deblur)
    elif aligner == "dada2":
        m9_exporting(m3_out_dada2,m2b_dada2_out,m3_table_filt_dada2,script_path, m4_out_dada2)
    else:
        m9_exporting(m3_out_deblur,m2b_deblur_out,m3_table_filt_deblur,script_path, m4_out_deblur)
        m9_exporting(m3_out_dada2,m2b_dada2_out,m3_table_filt_dada2,script_path, m4_out_dada2)

    m3_out_deblur = '%s/rep_seqs_filt.qza' % m2b_deblur_out
    m3_out_dada2 = '%s/rep_seqs_filt.qza' % m2b_dada2_out
    m3_table_filt_deblur = '%s/table_filt.qza' % m2b_deblur_out
    m3_table_filt_dada2 = '%s/table_filt.qza' % m2b_dada2_out
    ##task
    #m10 do picrust2 analysis
    if aligner == "deblur":
        m10_picrust2(m2b_deblur_out,m3_table_filt_deblur,m2b_deblur_out,nnodes, metadata_category,script_path)
    elif aligner == "dada2":
        m10_picrust2(m2b_dada2_out,m3_table_filt_dada2,m3_out_dada2,nnodes, metadata_category,script_path)
    else:
        m10_picrust2(m2b_deblur_out,m3_table_filt_deblur,m2b_deblur_out,nnodes, metadata_category,script_path)
        m10_picrust2(m2b_dada2_out,m3_table_filt_dada2,m3_out_dada2,nnodes, metadata_category,script_path)
    
    #还有很多其他的功能应该还是缺失的，需要逐步补充


