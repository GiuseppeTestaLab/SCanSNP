
import argparse
import os

parser = argparse.ArgumentParser(description='Process some integers.')


parser.add_argument('--vcf', dest='vcf', help='joint vcf with all mixed genotypes SNPs',required=True,type=str)


#conditionally optionals
parser.add_argument('--bam', dest='bam', help='for matrixgen mode: BAM file of mixed genotypes',type=str)
parser.add_argument('--barcodes', dest='barcodes', help='for matrixgen mode: FILTERED barcodes file in standard 10x .tsv format',type=str)
parser.add_argument('--pickle', dest='countpath',type=str,help='for deconvolution mode: path to folder containing Counts.npz countBarcodes.tsv and countLoci.tsv')

#Optionals
parser.add_argument('--threads', dest='nthreads', help='threads to be used',default=10,type=str)
parser.add_argument('--outdir', dest='outdir', help='use',default="currentwd",type=str)
parser.add_argument('--flagLowQual', dest='LowQual', help='Force-fit GMM to flag lowquality droplets',default=False,type=bool)
parser.add_argument('--mode', dest='mode', default="deconvolution",
	required=False,
	choices=["matrixgen","deconvolution","skipcount"],
	type=str,
	help='creating pileup matrix with ref and alt supporting reads')

args = parser.parse_args()

if args.mode == "matrixgen":
	if args.bam is None or args.barcodes is None:
		parser.error('please provide --bam and --barcodes for matrixgen mode')

if args.mode == "skipcount":
	if args.countpath is None:
		parser.error('please provide path to SparseD.pkl using --pickle for skipcount mode')

if args.mode is None or args.mode == "deconvolution":
	if args.bam is None or args.barcodes is None:
		parser.error('please provide --bam and --barcodes')



mode=args.mode
bamFile=args.bam
vcf=args.vcf
barcodesFILE=args.barcodes
nThreads=int(args.nthreads)
outdir=args.outdir
countpath=args.countpath
LowQual=args.LowQual

if outdir == "currentwd":
	outdir = os.getcwd()
else:
	outdir=args.outdir
	
if str(outdir).endswith(' '):
	outdir=str(outdir)[:-1]

if str(outdir).endswith('/'):
	outdir=str(outdir)[:-1]

if str(countpath).endswith('/'):
	countpath=str(countpath)[:-1]




import pysam
import pandas as pd
import time
import itertools
import io
import numpy as np
import sys
from collections import defaultdict
import re
import math
from collections import Counter
from sklearn.linear_model import LogisticRegression
import scipy
import pickle

#sys.path.append('/hpcnfs/home/ieo4777/gitted/SCanSNP_v1/SCanSNP/')

from VCFUtils import *
from Wrappers import *


if __name__ == "__main__":
	#Creation of differen loci subsets
	cleanLoci = LociPreClean(vcf)
	MildcleanLoci = LociPreClean_milds(vcf)
	CleanSingularLoci,SingularLoci_Alt,SingularLoci_Ref = SingularLociSCan(vcf,cleanLoci)
	
	#barcodeList=pd.read_csv(barcodesFILE, header=None, names=["b"])["b"].tolist()
	
	#Genotypes map creation
	GenotypesDF = creategenotypeDF(vcf)
	
	
	if mode == "matrixgen":
		'''
		Performing only count
		'''
		barcodeList=pd.read_csv(barcodesFILE, header=None, names=["b"])["b"].tolist()
		SparseD = CountsMatrices(CleanSingularLoci, cleanLoci,
			MildcleanLoci, GenotypesDF,
			barcodeList, vcf,
			nThreads, bamFile)
		#Save ReadCounts in pickle
		with open(outdir+ '/SparseD.pkl', 'wb') as f:
			pickle.dump(SparseD, f, pickle.HIGHEST_PROTOCOL)
			
		
		
	elif mode == "skipcount":
		'''
		Performing only deconvolution based on provided counts
		'''
		#Load existing counts
		with open(str(countpath), 'rb') as f:
			SparseD = pickle.load(f)
		#Deconvolution
		Cell_IDs = deconvolution(SparseD, vcf, GenotypesDF,outdir, LowQual)
		Cell_IDs.to_csv(outdir + "/Cell_IDs.tsv", sep = "\t", header = True, index = True)
		
		
	elif mode == "refineDBLs":
		print("ComingSoon")
		
		
	else:
		'''
		Perform both count and ceconvolution
		'''
		barcodeList=pd.read_csv(barcodesFILE, header=None, names=["b"])["b"].tolist()
		#counting
		SparseD = CountsMatrices(CleanSingularLoci, cleanLoci,
			MildcleanLoci, GenotypesDF,
			barcodeList, vcf,
			nThreads, bamFile)
		#Deconvolution
		Cell_IDs = deconvolution(SparseD, vcf, GenotypesDF,outdir,LowQual)
		Cell_IDs.to_csv(outdir + "/Cell_IDs.tsv", sep = "\t", header = True, index = True)
