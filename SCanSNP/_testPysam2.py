import pandas as pd
import subprocess
import pysam


bamFile = "/data/possorted_genome_bam.filtered.bam"
lociDF = pd.read_csv("/data/temp.loci.tsv", sep="\t" )
BarcodeSet = set(pd.read_csv("/data/barcodes.tsv.gz", sep="\t" ).values.flatten().tolist())

readList = []
print(len(readList))

bam=pysam.AlignmentFile(bamFile, "rb")
for locus in lociDF.values.tolist():
    for read in bam.fetch(locus[0], locus[1]-1, locus[1]):
    #Check for position coverage
        CB=read.get_tag("CB")
        try:
            position = read.positions.index(int(locus[1]) - 1)
        except:
            continue
        if read.query_alignment_qualities[position] < 20:
            continue
        if read.mapq < 2:
            continue
        if read.rlen < 30:
            continue
        redBase = read.query_alignment_sequence[position]
        if ( redBase == locus[3] or redBase == locus[4]):
            readList.append([locus[0]+ "_"+str(locus[1]), CB, redBase, locus[3], locus[4] ])

bam.close()
print(len(set([i[1] for  i in readList])))
print(len(readList))


