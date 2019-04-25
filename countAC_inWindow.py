## 2019-1-5
## write this code to Count A,C for 6mA methylation and 5mC methylation

# import rpy2.robjects as robjects
r_script = """
    library(Rsamtools)
    library(GenomicFeatures)
    library(data.table)
    m4C_gff <- fread("/home/wanghm/whm/pacbio_data/Pacbio/PH-1/6mA/official_methylationResult/m4C.gff", skip = 1)
    m6A_gff <- fread("/home/wanghm/whm/pacbio_data/Pacbio/PH-1/6mA/official_methylationResult/m6A.gff", skip = 1)
    m4CGrange = GRanges(
  seqnames = m4C_gff$V1,
  ranges = IRanges(start = m4C_gff$V4, end = m4C_gff$V5),
  strand = m4C_gff$V7
)
    m6AGrange = GRanges(
  seqnames = m6A_gff$V1,
  ranges = IRanges(start = m6A_gff$V4, end = m6A_gff$V5),
  strand = m6A_gff$V7
)
    txdb <- makeTxDbFromGFF(file = "/home/wanghm/whm/1FG/35genome/Fusarium_graminearum.RR1.35.chr.gff3")
    ge <- genes(txdb)
    start(ge) <- start(ge)-500; end(ge) <- end(ge) + 500
    geneMethylationCount <- countOverlaps(ge, m6AGrange)      # remove no methelation genes
    ge <- ge[names(ge) %in% names(geneMethylationCount[geneMethylationCount!=0])]
    ss <- getSeq(FaFile("/home/wanghm/whm/1FG/35genome/Fusarium_graminearum.RR1_with_mito.fa"), ge)
    writeXStringSet(ss, "/home/wanghm/Desktop/ge.fasta")
"""

"""
python call Rscript
"""
# robjects.r(r_script)

from Bio import SeqIO
from Bio import Seq
import re
import time
import sys
import numpy as np
import multiprocessing
from itertools import islice

def chunk(it, size):
    it = iter(it)
    return iter(lambda: tuple(islice(it, size)), ())

def getAcountMatrix(p,type):
    arr = np.empty(shape=(0, 70))  # make a empty with colnum = 70
    chunk_sequences = record_list[p]
    for seq in chunk_sequences:
        seq_length = len(seq)
        geneBodyLength = seq_length - 1000
        binSize = int(geneBodyLength/50)
        geneBody_lastBinLength = binSize + geneBodyLength%50
        row_list = []
        pos = 0
        for i in range(1,71):
            if i <= 10 or i >= 61:
                binSeq = seq.seq[pos:pos+50]
                pos = pos + 50
            elif i == 60:
                binSeq = seq.seq[pos:pos+geneBody_lastBinLength]
                pos = pos + geneBody_lastBinLength

            else:
                binSeq = seq.seq[pos:pos+binSize]
                pos = pos + binSize

            row_list.append((binSeq.count("A")+binSeq.count("T")) if type == "m6A" else (binSeq.count("C")+binSeq.count("G")))
            # row_list.append(binSeq.count("T") if type == "m6A" else binSeq.count("G"))

        arr = np.row_stack((arr, row_list))
    arr1 = arr.sum(axis=0)    # sum each cols of numpy array
    return arr1

if __name__ == "__main__":
    fh = "/home/wanghm/Desktop/gebody.fasta"
    seqs = SeqIO.parse(fh, format="fasta")
    records = SeqIO.to_dict(seqs)
    record_list = list(chunk(records.values(), int(len(records)/4)))
    res = []
    processor = 4
    pool = multiprocessing.Pool(processes=processor)
    ss = range(processor+1 if len(records)%4 != 0 else processor)
    for i in ss:
        res.append(pool.apply_async(func=getAcountMatrix, args=(i, sys.argv[1]))) # asynchronization
    aa = zip(*(res[i].get().tolist() for i in ss))
    with open("/home/wanghm/wanghm_R/methylation/{0}_genebody".format(sys.argv[1]), "w+") as f:
        for i in aa:
            # print(i)
            total = sum(i)
            f.writelines(str(total)+"\n")


