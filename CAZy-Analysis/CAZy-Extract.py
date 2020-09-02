import os
import sys
import re 
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

Dir = "/home/hongbinliu/dbcan2_report"
AssDir = "/home/junyuchen/Lab/Liuhongbin/Result/Assemble-all-1"
OutDir = "/home/junyuchen/Lab/Liuhongbin/Result/GeneCatalogExtraction-2"

def ExtractOverview(overview):
    df = pd.read_table(overview)
    return df["Gene ID"].tolist()

def ExtractUniInput(file, GeneID):
    idList = []
    startList = []
    endList = []
    strandList = []
    for seq in SeqIO.parse(file, "fasta"):
        if seq.id in GeneID:
            desc = seq.description.split(" # ")
            idList.append(desc[0])
            startList.append(desc[1])
            endList.append(desc[2])
            strandList.append(desc[3])
    return idList, startList, endList, strandList

#idList plus Prodigal List
#remove star,end,strand
def ExtractGeneCatalog(fasta, ID, OutDir):
    AllSeq = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    geneCatalog = []
    for i in range(len(idList)):
        ID = re.sub(r'\_[0-9]{0,4}$', '', idList[i]) 
        if ID in AllSeq.keys():
            #print(strandList[i])
            seq = AllSeq[ID]
            #print(strandList[i])
            if strandList[i] == "1":
                Seq = SeqRecord(seq.seq[int(startList[i]) : int(endList[i])])
                #Seq.id = seq.id
                Seq.id = idList[i]
                #print(seq.id)
                #print("1")
                Seq.description = idList[i]
                #Seq.description = seq.id + " " + startList[i] + " " + endList[i] + " " + strandList[i]
            elif strandList[i] == "-1":
                Seq = SeqRecord(seq.seq[int(startList[i]) : int(endList[i])].reverse_complement())
                Seq.id = seq.id
                #print("-1")
                Seq.description = seq.id + " " + startList[i] + " " + endList[i] + " " + strandList[i]        
            geneCatalog.append(Seq)
    SeqIO.write(geneCatalog, os.path.join(OutDir, ID + ".fasta"), "fasta")
    print(ID)
    
for run in os.listdir(Dir):
    runDir = os.path.join(Dir, run)
    overviewDir = os.path.join(runDir, "overview.txt")
    uniInputDir = os.path.join(runDir, "uniInput")
    #print(run)
    #print(runDir)
    GeneID = ExtractOverview(overviewDir)
    idList, startList, endList, strandList = ExtractUniInput(uniInputDir, GeneID)
    fastaDir = os.path.join(AssDir, run, "scaffolds.fasta")
    ExtractGeneCatalog(fastaDir, run, OutDir)