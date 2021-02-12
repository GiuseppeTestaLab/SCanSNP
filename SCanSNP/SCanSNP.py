#!/usr/bin/env python

import argparse
import os

parser = argparse.ArgumentParser(description='Process some integers.')

parser.add_argument('--mode', dest='mode', default="matrixgen",
	required=True,
	choices=["matrixgen","deconvolution"],
	type=str,
	help='creating pileup matrix with ref and alt supporting reads')
parser.add_argument('--vcf', dest='vcf', help='joint vcf with all mixed genotypes SNPs',required=True,type=str)


#conditionally optionals
parser.add_argument('--BAM', dest='bam', help='for matrixgen mode: BAM file of mixed genotypes',type=str)
parser.add_argument('--barcodes', dest='barcodes', help='for matrixgen mode: FILTERED barcodes file in standard 10x .tsv format',type=str)
parser.add_argument('--counts', dest='countpath',type=str,help='for deconvolution mode: path to folder containing Counts.npz countBarcodes.tsv and countLoci.tsv')

#Optionals
parser.add_argument('--threads', dest='nthreads', help='threads to be used',default=10,type=str)
parser.add_argument('--outdir', dest='outdir', help='use',default="currentwd",type=str)


args = parser.parse_args()

if args.mode == "matrixgen":
	if args.bam is None or args.barcodes is None:
		parser.error('please provide --bam and --barcodes for matrixgen mode')

if args.mode == "deconvolution":
	if args.countpath is None:
		parser.error('please provide --counts for deconvolution mode')





mode=args.mode
bamFile=args.bam
vcf=args.vcf
barcodesFILE=args.barcodes
nThreads=int(args.nthreads)
outdir=args.outdir
countpath=args.countpath

if outdir == "currentwd":
	writePath=wd = os.getcwd()
else:
	writePath=args.outdir

if str(writePath).endswith('/'):
	writePath=str(writePath)[:-1]

if str(countpath).endswith('/'):
	countpath=str(countpath)[:-1]



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




from Pileup import *
from VCFUtils import *
from DBLsutils import *
from ComputeLLK  import *








if __name__ == "__main__":
	if mode == "matrixgen":
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
		GenotypeChunkDict = ChunkMaker(GenotypesDF, nThreads, OmniIndex)

		print('Starting the pileup...')
		start = time.time()

		#FireUp main Pileupper
		results = []
		pool=Pool(nThreads)

		for chunkId in GenotypeChunkDict.keys():

			result = pool.apply_async(ReadCounter, (chunkId, bamFile, barcodeList, GenotypeChunkDict))
			results.append(result)


		pool.close()
		pool.join()

		print('Pileup took', time.time()-start, 'seconds.')

		#Mergion pooled process results
		Reads=list(itertools.chain.from_iterable(([result.get() for result in results])))
		start = time.time()
		countDF=DFMaker(Reads, barcodeList)
		print('It took', time.time()-start, 'seconds. To create the Count matrix')

		#Saving countsMatrix
		SparseCounts=scipy.sparse.csc_matrix(countDF.values)
		countBarcodes=countDF.columns
		countLoci=countDF.index

		scipy.sparse.save_npz(writePath + '/Counts', SparseCounts)
		pd.Series(countBarcodes).to_csv(writePath + "/countBarcodes.tsv", sep = "\t", header = False, index = False)
		pd.Series(countLoci).to_csv(writePath + "/countLoci.tsv", sep = "\t", header = False, index = False)
	elif mode == "deconvolution":
		#CountDF reconstruction
		SparseCounts=scipy.sparse.load_npz(str(countpath) + '/Counts.npz')
		Barcodes=pd.read_csv(str(countpath) + '/countBarcodes.tsv', header=None, names=["barcodes"])
		Loci=pd.read_csv(str(countpath) + '/countLoci.tsv', header=None, names=["loci"])


		countDF = pd.DataFrame(SparseCounts.todense())
		countDF.set_index(Loci["loci"], inplace = True)
		countDF.columns = Barcodes["barcodes"]

		barcodeList= list(set([barcode.split("_")[0] for barcode in Barcodes["barcodes"]]))

		#Genotypes map creation
		GenotypesDF = creategenotypeDF(vcf)

		#Pre-cleaning loci
		cleanLoci = list(set(LociPreClean(vcf)).intersection(list(countDF.index)))
		MildcleanLoci = list(set(LociPreClean_milds(vcf)).intersection(list(countDF.index)))

		#Creation of loci subsets
		DiffOnlyIndex=DiffOnlyIndexMaker(vcf, cleanLoci)
		CleanSingularLoci,SingularLoci_Alt,SingularLoci_Ref = SingularLociSCan(vcf,cleanLoci)
		DBLSpecificSingularLociHomoz=DoubletSpecificSingularLociScan(vcf,GenotypesDF, cleanLoci)[1]
		DBLSpecificSingularLociHeteroAlt=DoubletSpecificSingularLociScan(vcf,GenotypesDF, cleanLoci)[2]
		DBLSpecificSingularLociHeteroRef=DoubletSpecificSingularLociScan(vcf,GenotypesDF, cleanLoci)[3]



		LikeliHoodsDF=ComputeLikelihood(vcf,GenotypesDF,countDF,barcodeList,DiffOnlyIndex)


		BestBarcodeID=LikeliHoodsDF.idxmax(axis = 1)
		BestBarcodeID.to_csv(writePath + "/Barcode-ID.tsv", sep = "\t", header = False, index = True)


		BestInDropDict = defaultdict(list)
		for key, value in BestBarcodeID.to_dict().items():
			BestInDropDict[value].append(key)

		#Estimation and regression of genotyped mix per droplet
		SingularLociScoreDF = SingularLociScore( SingularLoci_Alt, SingularLoci_Ref, countDF, barcodeList, GenotypesDF,vcf)
		DBLsDF=NoiseRregression(BestInDropDict, barcodeList, vcf, SingularLociScoreDF, BestBarcodeID)

		#Creation of putative DBL-to-barcode dict
		DropToDBLDict = defaultdict(list)
		Barcodes_values = list(DBLsDF.index)
		DBLtuples_keys=list(zip(DBLsDF["FirstID"],DBLsDF["SecondID"]))
		for key, value in dict(zip(Barcodes_values, DBLtuples_keys)).items():
			DropToDBLDict[value].append(key)

		DBLmetricsDF = DBLScore(DropToDBLDict,countDF,DBLSpecificSingularLociHomoz,DBLSpecificSingularLociHeteroRef,DBLSpecificSingularLociHeteroAlt,vcf, GenotypesDF,SingularLociScoreDF)

		#Editing DFs before concat
		Contributions = SingularLociScoreDF[ExtractSamples(vcf)]
		DBLsContributions = DBLmetricsDF
		BestIDs = DBLsDF.replace("ID_","", regex = True)

		pd.concat([Contributions,DBLsContributions,BestIDs], axis = 1).to_csv(writePath + "/DBLmetricsDF.tsv", sep = "\t", header = True, index = True)
