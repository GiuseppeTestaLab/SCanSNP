BAM: "/hpcnfs/scratch/temporary/Dav_vc/0_CellrangerAlign/Sample_S20272_157/outs/possorted_genome_bam.bam"
goodBarcodes: "/hpcnfs/scratch/temporary/Dav_vc/0_CellrangerAlign/Sample_S20272_157/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
VCF: "/hpcnfs/scratch/temporary/Dav_vc/2_Combining/Sample_S20273_157/output/combinedVCF/jointGenotype.g.vcf"
field: "PL"
mixedGenotypes: "4"
barcodesMap: NULL
FilteredFeaturesPath: "/hpcnfs/scratch/temporary/Dav_vc/0_CellrangerAlign/Sample_S20272_157/outs/filtered_feature_bc_matrix/"
nThreads=20


#!/usr/bin/env python

import argparse
import os




writePath="/hpcnfs/home/ieo4777/temp"
bamFile="/hpcnfs/scratch/temporary/Dav_vc/0_CellrangerAlign/Sample_S20272_157/outs/possorted_genome_bam.bam"
vcf="/hpcnfs/scratch/temporary/Dav_vc/2_Combining/Sample_S20273_157/output/combinedVCF/jointGenotype.g.vcf"
barcodesFILE="/hpcnfs/home/ieo4777/temp/barcodes.tsv"
nThreads=15
# outdir=args.outdir
# countpath=args.countpath





import pysam
import pandas as pd
import time
import itertools
import io
import numpy as np
from multiprocessing import Pool
import sys
from collections import defaultdict
import re
import math
from collections import Counter
from sklearn.linear_model import LogisticRegression
import scipy






print('Loading Data...')
barcodes=pd.read_csv(barcodesFILE, header=None, names=["b"])
barcodeList=barcodes["b"].tolist()

#Creation of differen loci subsets
cleanLoci = LociPreClean(vcf)
MildcleanLoci = LociPreClean_milds(vcf)
CleanSingularLoci,SingularLoci_Alt,SingularLoci_Ref = SingularLociSCan(vcf,cleanLoci)

#Genotypes map creation
GenotypesDF = creategenotypeDF(vcf)

#Creation of list with all loci to be pileup on
OmniIndex = list(set(CleanSingularLoci + DiffOnlyIndexMaker(vcf, cleanLoci) + DoubletSpecificSingularLociScan(vcf,GenotypesDF, MildcleanLoci)[0]))

print('Splitting Variants into chunks...')
GenotypeChunkDict=ChunkMaker(GenotypesDF, nThreads, OmniIndex)
#Converting genotypeschunk into list
GenotypeChunkList,GenotypeChunkIndexesList = FlattenDict(GenotypeChunkDict)
del GenotypeChunkDict

print('Starting the pileup...')
start = time.time()

#FireUp main Pileupper
results = []
pool=Pool(nThreads)

for chunkIdx in GenotypeChunkIndexesList:
	result = pool.apply_async(ReadCounter, (chunkIdx, bamFile, barcodeList))
	results.append(result)


pool.close()
pool.join()

print('Pileup took', time.time()-start, 'seconds.')

#Mergion pooled process results
Reads=list(itertools.chain.from_iterable(([result.get() for result in results])))
##readList = Reads
##ReadsBU = Reads
##Reads = ReadsBU

start = time.time()
SparseD = {}
SparseD["sparse_Ref"],SparseD["sparse_Alt"],SparseD["Locus"],SparseD["Barcode"] = DFMaker(Reads, barcodeList)

# #Saving countsMatrix
# SparseCounts=scipy.sparse.csc_matrix(countDF.values)
# countBarcodes=countDF.columns
# countLoci=countDF.index
#
# scipy.sparse.save_npz(writePath + '/Counts', SparseCounts)
# pd.Series(countBarcodes).to_csv(writePath + "/countBarcodes.tsv", sep = "\t", header = False, index = False)
# pd.Series(countLoci).to_csv(writePath + "/countLoci.tsv", sep = "\t", header = False, index = False)


##Part 2
##Part 2
##Part 2

#Pre-cleaning loci CHECK!!
#Pre-cleaning loci CHECK!!
cleanLoci = list(set(LociPreClean(vcf)).intersection(list(SparseD["Locus"])))
MildcleanLoci = list(set(LociPreClean_milds(vcf)).intersection(list(SparseD["Locus"])))

#Setting Loci categories
DiffOnlyIndex=DiffOnlyIndexMaker(vcf, cleanLoci)
CleanSingularLoci,SingularLoci_Alt,SingularLoci_Ref = SingularLociSCan(vcf,cleanLoci)
DBLSpecificSingularLociHomoz=DoubletSpecificSingularLociScan(vcf,GenotypesDF, cleanLoci)[1]
DBLSpecificSingularLociHeteroAlt=DoubletSpecificSingularLociScan(vcf,GenotypesDF, cleanLoci)[2]
DBLSpecificSingularLociHeteroRef=DoubletSpecificSingularLociScan(vcf,GenotypesDF, cleanLoci)[3]


LikeliHoodsDF=ComputeLikelihood(vcf,GenotypesDF,SparseD,DiffOnlyIndex)


BestBarcodeID=LikeliHoodsDF.idxmax(axis = 1)
BestBarcodeID.to_csv(writePath + "/Barcode-ID.tsv", sep = "\t", header = False, index = True)

BestInDropDict = defaultdict(list)
for key, value in BestBarcodeID.to_dict().items():
	BestInDropDict[value].append(key)



SingularLociScoreDF = SingularLociCNTR( SingularLoci_Alt, SingularLoci_Ref, SparseD, barcodeList, GenotypesDF,vcf)
DBLsDF=NoiseRregression(BestInDropDict, barcodeList, vcf, SingularLociScoreDF, BestBarcodeID)

#Creation of putative DBL-to-barcode dict
DropToDBLDict = defaultdict(list)
Barcodes_values = list(DBLsDF.index)
DBLtuples_keys=list(zip(DBLsDF["FirstID"],DBLsDF["SecondID"]))
for key, value in dict(zip(Barcodes_values, DBLtuples_keys)).items():
	DropToDBLDict[value].append(key)
