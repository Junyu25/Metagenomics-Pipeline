#coding:utf-8
from __future__ import print_function
import math
import glob,os,re,argparse


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
    os.system('mkdir humann2_out')
    import math
    nthread=str(int(math.ceil(float(int(nnodes))/int(njobs))))
    #this is the most time consuming step
    os.system("parallel -j %s 'humann2 --threads %s --input {} --output humann2_out/{/.}' ::: cat_reads/*fastq" % (njobs,nthread))
    #3.2 Merge individual sample data together
    os.system('mkdir humann2_final_out')
    os.system('humann2_join_tables -s --input humann2_out/ --file_name pathabundance --output humann2_final_out/humann2_pathabundance.tsv')
    os.system('humann2_join_tables -s --input humann2_out/ --file_name pathcoverage --output humann2_final_out/humann2_pathcoverage.tsv')
    os.system('humann2_join_tables -s --input humann2_out/ --file_name genefamilies --output humann2_final_out/humann2_genefamilies.tsv')
    #3.3 Table output normalization
    #Re-normalize gene family and pathway abundances (so that all samples are in units of copies per million).
    os.system('humann2_renorm_table --input humann2_final_out/humann2_pathabundance.tsv --units cpm --output humann2_final_out/humann2_pathabundance_cpm.tsv')
    os.system('humann2_renorm_table --input humann2_final_out/humann2_genefamilies.tsv --units cpm --output humann2_final_out/humann2_genefamilies_cpm.tsv')
    #3.4 Separate out taxonomic contributions
    #Split HUMAnN2 output abundance tables in stratified and unstratified tables (stratified tables include the taxa associated with a functional profile).
    os.system('humann2_split_stratified_table --input humann2_final_out/humann2_pathabundance_cpm.tsv --output humann2_final_out')
    os.system('humann2_split_stratified_table --input humann2_final_out/humann2_genefamilies_cpm.tsv --output humann2_final_out')
    os.system('humann2_split_stratified_table --input humann2_final_out/humann2_pathcoverage.tsv --output humann2_final_out')
    #3.5 Format STAMP function file
    #Convert unstratified HUMAnN2 abundance tables to STAMP format by changing header-line. These commands remove the comment character and the spaces in the name of the first column. Trailing descriptions of the abundance datatype are also removed from each sample's column name.
    os.system("sed 's/_Abundance-RPKs//g' humann2_final_out/humann2_genefamilies_cpm_unstratified.tsv | sed 's/# Gene Family/GeneFamily/' > humann2_final_out/humann2_genefamilies_cpm_unstratified.spf")
    os.system("sed 's/_Abundance//g' humann2_final_out/humann2_pathabundance_cpm_unstratified.tsv | sed 's/# Pathway/Pathway/' > humann2_final_out/humann2_pathabundance_cpm_unstratified.spf")
    #3.6 Extract MetaPhlAn2 taxonomic compositions
    #Since HUMAnN2 also runs MetaPhlAn2 as an initial step, we can use the output tables already created to get the taxonomic composition of our samples. First we need to gather all the output MetaPhlAn2 results per sample into a single directory and then merge them into a single table using MetaPhlAn2's merge_metaphlan_tables.py command. After this file is created we can fix the header so that each column corresponds to a sample name without the trailing "_metaphlan_bugs_list" description. Note that MetaPhlAn2 works best for well-characterized environments, like the human gut, and has low sensitivity in other environments.
    os.system('mkdir metaphlan2_out')
    os.system('cp humann2_out/*/*/*metaphlan_bugs_list.tsv metaphlan2_out/')
    os.system('merge_metaphlan_tables.py metaphlan2_out/*metaphlan_bugs_list.tsv > metaphlan2_merged.txt')
    os.system("sed -i 's/_metaphlan_bugs_list//g' metaphlan2_merged.txt")
    #3.7 Format STAMP taxonomy file
    #Lastly we can convert this MetaPhlAn2 abundance table to STAMP format
    os.system('perl %s/metaphlan_to_stamp.pl metaphlan2_merged.txt > metaphlan2_merged.spf' % (script_path))

def m3_humann2_rmtmp(nnodes, script_path,njobs):
    """
    :nnodes: the number of nodes to parallel
    :njobs: the number of jobs for humann2 to run in parallel
    :nthread: the number of threads to run per job
    :script_path: the path of custome script, default is /home/yxtan/HUMANN2_SOP_scripts/
    :return: None
    """
    #3.1 Run HUMAnN2
    #Run humann2 with parallel to calculate abundance of UniRef90 gene families and MetaCyc pathways. If you are processing environment data (e.g. soil samples) the vast majority of the reads may not map using this approach. Instead, you can try mapping against the UniRef50 database (which you can point to with the --protein-database option).
    os.system('mkdir humann2_out')
    import math
    nthread=str(int(math.ceil(float(int(nnodes))/int(njobs))))
    #this is the most time consuming step
    os.system("parallel -j %s 'humann2 --threads %s --input {} --output humann2_out/{/.}' ::: cat_reads/*fastq" % (njobs,nthread))
    #3.2 Merge individual sample data together
    os.system('mkdir humann2_final_out')
    os.system('humann2_join_tables -s --input humann2_out/ --file_name pathabundance --output humann2_final_out/humann2_pathabundance.tsv')
    os.system('humann2_join_tables -s --input humann2_out/ --file_name pathcoverage --output humann2_final_out/humann2_pathcoverage.tsv')
    os.system('humann2_join_tables -s --input humann2_out/ --file_name genefamilies --output humann2_final_out/humann2_genefamilies.tsv')
    #3.3 Table output normalization
    #Re-normalize gene family and pathway abundances (so that all samples are in units of copies per million).
    os.system('humann2_renorm_table --input humann2_final_out/humann2_pathabundance.tsv --units cpm --output humann2_final_out/humann2_pathabundance_cpm.tsv')
    os.system('humann2_renorm_table --input humann2_final_out/humann2_genefamilies.tsv --units cpm --output humann2_final_out/humann2_genefamilies_cpm.tsv')
    #3.4 Separate out taxonomic contributions
    #Split HUMAnN2 output abundance tables in stratified and unstratified tables (stratified tables include the taxa associated with a functional profile).
    os.system('humann2_split_stratified_table --input humann2_final_out/humann2_pathabundance_cpm.tsv --output humann2_final_out')
    os.system('humann2_split_stratified_table --input humann2_final_out/humann2_genefamilies_cpm.tsv --output humann2_final_out')
    os.system('humann2_split_stratified_table --input humann2_final_out/humann2_pathcoverage.tsv --output humann2_final_out')
    #3.5 Format STAMP function file
    #Convert unstratified HUMAnN2 abundance tables to STAMP format by changing header-line. These commands remove the comment character and the spaces in the name of the first column. Trailing descriptions of the abundance datatype are also removed from each sample's column name.
    os.system("sed 's/_Abundance-RPKs//g' humann2_final_out/humann2_genefamilies_cpm_unstratified.tsv | sed 's/# Gene Family/GeneFamily/' > humann2_final_out/humann2_genefamilies_cpm_unstratified.spf")
    os.system("sed 's/_Abundance//g' humann2_final_out/humann2_pathabundance_cpm_unstratified.tsv | sed 's/# Pathway/Pathway/' > humann2_final_out/humann2_pathabundance_cpm_unstratified.spf")
    #3.6 Extract MetaPhlAn2 taxonomic compositions
    #Since HUMAnN2 also runs MetaPhlAn2 as an initial step, we can use the output tables already created to get the taxonomic composition of our samples. First we need to gather all the output MetaPhlAn2 results per sample into a single directory and then merge them into a single table using MetaPhlAn2's merge_metaphlan_tables.py command. After this file is created we can fix the header so that each column corresponds to a sample name without the trailing "_metaphlan_bugs_list" description. Note that MetaPhlAn2 works best for well-characterized environments, like the human gut, and has low sensitivity in other environments.
    os.system('mkdir metaphlan2_out')
    os.system('cp humann2_out/*/*/*metaphlan_bugs_list.tsv metaphlan2_out/')
    os.system('merge_metaphlan_tables.py metaphlan2_out/*metaphlan_bugs_list.tsv > metaphlan2_merged.txt')
    os.system("sed -i 's/_metaphlan_bugs_list//g' metaphlan2_merged.txt")
    #3.7 Format STAMP taxonomy file
    #Lastly we can convert this MetaPhlAn2 abundance table to STAMP format
    os.system('perl %s/metaphlan_to_stamp.pl metaphlan2_merged.txt > metaphlan2_merged.spf' % (script_path))

    #3.8 Remove temp file
    #os.system('rm -rf humann2_out/*/*temp*/')
    #os.system('rm -rf kneaddata_out/')
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='rd', type=str, required=True,
                        help="the tabular-table contains paths of the raw data")
    parser.add_argument('-n', '--node', dest='node', type=str, required=False, default='20',
                        help="the number of nodes to request")
    parser.add_argument('-sp', '--spath',dest='sp', type=str, required=False, default='/home/yxtan/HUMANN2_SOP_scripts/',
                        help="path of the custom scripts")
    parser.add_argument('-e', '--pair', dest='pair', type=str, required=False, default='True',
                        help="Is it pair-end seq data? Default is 'True'; Any other strings will be considered False")
    parser.add_argument('-j', '--jobs', dest='jobs', type=str, required=False, default='8',
                        help="The number of jobs run parallell in humann2 step. Default is '8'; It is bounded by the total number of memory available. Each job should have 16GB memory")
    parser.add_argument('-t', '--tpath', dest='tp', type=str, required=False, default='/home/LDlab/anaconda3/envs/humann2/share/trimmomatic-0.39-1/',
                        help="The path of the trimmomatic ref files. Default is /home/LDlab/anaconda3/envs/humann2/share/trimmomatic-0.39-1/")
    parser.add_argument('-d', '--dpath', dest='dp', type=str, required=False, default='/home/yxtan/ref_databases/GRCh38_PhiX_bowtie2_index/GRCh38_PhiX',
                        help="The path of the ref database. Default is /home/yxtan/ref_databases/GRCh38_PhiX_bowtie2_index/GRCh38_PhiX")
    parser.add_argument('-r', '--rmtmp', dest='rt', type=str, required=False, default='FALSE',
                        help="Whether to remove the temp files after running Humann2 to save space. Default is FALSE; Any other strings will be considered True")
    
    args = parser.parse_args()
    print('Usage example with minimum parameters: python /home/yxtan/HUMANN2_SOP_scripts/Metagenomics_HUMANN2.py -i sample_table.txt -n 4')
    rd_dir = os.path.abspath(args.rd)
    script_path = os.path.abspath(args.sp)
    trimmomatic_path = os.path.abspath(args.tp)
    refdb_path = os.path.abspath(args.dp)
    pair_end = args.pair
    nnodes = args.node
    njobs = args.jobs
    rmtmp = args.rt
    
    #需要检查输入的参数是否正确，主要是路径是否存在    
    if not os.path.isfile(rd_dir):
        print('Input sample table is not exist; Exit now.')
        exit(0)
    '''
    if not os.path.isfile("%s.1.bt2" % (refdb_path)):
        print('The ref db bt2 file is not exist; Exit now.')
        exit(0)
    '''
    if not os.path.isdir(script_path):
        print('The folder of custom scripts is not exist; Exit now.')
        exit(0)
    if not os.path.isdir(trimmomatic_path):
        print('The folder of trimmomatic files is not exist; Exit now.')
        exit(0)
    

    #1.2 Inspect read quality
    #if pair_end == 'True':
     #   m1_fastqc(ALL_list,nnodes)
    #else:
     #   m1_fastqc(F_list,nnodes)
    
    #2. Read Quality-Control and Contaminant Screens and connect to a long read
    #if pair_end == 'True':
     #   m2_run_kneaddata(nnodes, F_list,R_list,script_path, trimmomatic_path, njobs, refdb_path)
    #else:
     #   m2_run_kneaddata_single(nnodes, F_list,script_path, trimmomatic_path, njobs, refdb_path)
    
    #3. Determine Functions with HUMAnN2
    if pair_end == 'True':
        m3_humann2(nnodes, script_path, njobs)
    else:
        m3_humann2_rmtmp(nnodes, script_path,njobs)
    

