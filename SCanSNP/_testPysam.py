import pandas as pd
import subprocess
import pysam


bamFile = "/data/possorted_genome_bam.bam"
lociDF = pd.read_csv("/tmp/temp.loci.tsv", sep="\t" )
BarcodeSet = set(pd.read_csv("/data/barcodes.tsv.gz", sep="\t" ).values.flatten().tolist())

readList = []
print("ciao")


bam=pysam.AlignmentFile(bamFile, "rb")
for locus in lociDF.values.tolist():
    for read in bam.fetch(locus[0], locus[1]-1, locus[1]):
    #Check for position coverage
        try:
            CB=read.get_tag(barcodetag)
        except:
            continue
        if not set([CB]).intersection(BarcodeSet):
            continue
        try:
            position = read.positions.index(int(locus[1]) - 1)
        except:
            continue
        if read.query_alignment_qualities[position] < 20:
            continue
        if read.is_secondary:
            continue
        if read.mapq < 2:
            continue
        if read.rlen < 30:
            continue
        if read.flag == 3844:
            continue
        redBase = read.query_alignment_sequence[position]
        if ( redBase == locus[3] or redBase == locus[4]):
            readList.append([locus[0]+ "_"+str(locus[1]), CB, redBase, locus[3], locus[4] ])
        print(len(readList))

bam.close()
print(len(readList))

