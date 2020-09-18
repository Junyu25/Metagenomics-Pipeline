import os

outFileList = []
for subdir, dirs, files in os.walk(Dir):
    blastp = ""
    for file in files:
        if file.endswith(".fastq_blastp.tsv"):
            blastp = os.path.join(subdir, file)
            df = pd.read_table(blastp, header=None)
            fileName = os.path.split(blastp)[1].replace(".fastq_blastp.tsv", "")
            df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
            df1 = df.loc[(df["length"] >= 25) & (df["pident"] >= 35)]


