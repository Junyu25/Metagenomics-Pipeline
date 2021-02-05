#!/bin/bash
#PBS -N humann2
#PBS -l nodes=3:ppn=36 
#PBS -l walltime=2400:00:00
#PBS -j oe
#PBS -q batch

##### "select" define your physical cores ; 
##### "ncpus" define your logical cores ;
##### "mem" define your used memory every "select";

echo Running on host `hostname`
echo Starting Time is `date`
echo Directory is `pwd`
starttime=$(date +"%s")

#actiavte env
source /home/junyuchen/Biosoft/anaconda3/bin/activate /home/junyuchen/Biosoft/anaconda3/envs/humann3
#set database
humann_config --update database_folders utility_mapping /home/hongbinliu/humann_dbs/utility_mapping
humann_config --update database_folders protein /home/hongbinliu/humann_dbs/uniref
humann_config --update database_folders nucleotide /home/hongbinliu/humann_dbs/chocophlan

python /home/junyuchen/Lab/Metagenomics-Pipeline/Scripts/humann3.py \
    -i /home/junyuchen/Lab/Metagenomics-Pipeline/data \
    -o /home/junyuchen/Lab/Metagenomics-Pipeline/result \
    -t 36 \
    -j 6 \
    > log.txt #logfile
# -t  number of threads you want to run per job/sample
# -j  number of jobs/sample you want to run in parallel
endtime=$(date +"%s")
diff=$(($endtime - $starttime))
echo Elapsed time is $(($diff/60)) minutes and $(($diff%60)) seconds.
