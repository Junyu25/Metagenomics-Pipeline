#!/bin/bash
#PBS -N Ass
#PBS -l nodes=3:ppn=36 
#PBS -l walltime=2400:00:00
#PBS -j oe
#PBS -q batch
#PBS Shell_Path_List=/bin/bash 

##### "select" define your physical cores ; 
##### "ncpus" define your logical cores ;
##### "mem" define your used memory every "select";
mpirun -np 108 ./a.out > acpimpilog.out
#mpirun -np 72 ./a.out > acpimpilog.out
mpirun -H node1,node4,node5 ./a.out > acpimpilog.out

echo Running on host `hostname`
echo Starting Time is `date`
echo Directory is `pwd`
starttime=$(date +"%s")

cd /home/junyuchen/Lab/Liuhongbin/Result/Assemble-all

source /home/junyuchen/Biosoft/anaconda3/bin/activate /home/junyuchen/Biosoft/anaconda3/envs/metabgc-r

python /home/junyuchen/Lab/Liuhongbin/ass.py -i /home/junyuchen/Lab/Liuhongbin/Result/kneaddata_out -o /home/junyuchen/Lab/Liuhongbin/Result/Assemble-all 

endtime=$(date +"%s")
diff=$(($endtime - $starttime))
echo Elapsed time is $(($diff/60)) minutes and $(($diff%60)) seconds.
