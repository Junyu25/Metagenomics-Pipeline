import os
import subprocess
from itertools import repeat
from multiprocessing import Pool, freeze_support

def gunzip(file, OutDir):
    cmd = "gunzip -k "+ filePath + " -d " + OutDir
    #print(cmd)
    subprocess.call(cmd, shell=True)

def unzipParallel(fileList, OutDir):
    #numOfprocess = int(20)
    pool = Pool(processes=100)
    pool.starmap(gunzip, zip(fileList, repeat(OutDir)))
    pool.close()
    pool.join()
    pool.terminate()


parser = argparse.ArgumentParser(description='RunKneadData')
parser.add_argument('-i', '--input', dest='fileDir', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OutDir', type=str, required=True,
                    help="the output path of reads")
args = parser.parse_args()

InFile = str(args.fileDir)
OutDir = os.path.abspath(args.OutDir)


#fileList = []
#print(inputfile)
df = pd.read_table(inputfile)
R1 = df["forward-absolute-filepath"].tolist() 
R2 = df["reverse-absolute-filepath"].tolist() 
#print(F_list)
RunCpFParallel(F_list, ouputDir)
