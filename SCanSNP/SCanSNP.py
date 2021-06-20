
import argparse
import os


parser = argparse.ArgumentParser(description='Process some integers.')


parser.add_argument('--vcf', dest='vcf', help='joint vcf with all mixed genotypes SNPs',required=True,type=str)
parser.add_argument('--barcodes', dest='barcodes', help='for matrixgen mode: FILTERED barcodes file in standard 10x .tsv format',required=True,type=str)


#conditionally optionals
parser.add_argument('--bam', dest='bam', help='for matrixgen mode: BAM file of mixed genotypes',type=str)
parser.add_argument('--pickle', dest='countpath',type=str,help='for deconvolution mode: path to folder containing Counts.npz countBarcodes.tsv and countLoci.tsv')
parser.add_argument('--filtered_matrix_path', dest='filtered_matrix', help='path/to/cellranger/filtered/matrix/',default=None,type=str)


#Optionals
parser.add_argument('--threads', dest='nthreads', help='threads to be used',default=10,type=str)
parser.add_argument('--raw_matrix_path', dest='raw_matrix', help='path/to/cellranger/unfiltered/matrix/',default=None,type=str)
parser.add_argument('--outdir', dest='outdir', help='use',default="currentwd",type=str)
#parser.add_argument('--flagLowQual', dest='LowQual', help='Force-fit GMM to flag lowquality droplets',default=False,type=bool)
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
		parser.error('please provide path to Counts.pkl using --pickle for skipcount mode')

if args.mode is None or args.mode == "deconvolution":
	if args.bam is None or args.barcodes is None:
		parser.error('please provide --bam and --barcodes')

if args.raw_matrix is not None:
	if args.filtered_matrix is None:
		parser.error('If --raw_matrix_path is specified also --filtered_matrix_path must be provided')


mode=args.mode
bamFile=args.bam
vcf=args.vcf
barcodesFILE=args.barcodes
nThreads=int(args.nthreads)
outdir=args.outdir
countpath=args.countpath
rawPath=args.raw_matrix
filteredPath=args.filtered_matrix


#mode="deconvolution"
#bamFile="/hpcnfs/scratch/temporary/Dav_vc/0.2_SynthII_creation/SYNTH5/Synth_Balanced_doublets/SilicoMultiplexed.dblAdded.bam"
#vcf="/hpcnfs/scratch/temporary/Dav_vc/0.1_Synth_creation_Dav_temp/Dav/1_Subsampling/CombineVCF/output/combinedVCF/jointGenotype.g.vcf"
#barcodesFILE="/hpcnfs/scratch/temporary/Dav_vc/0.2_SynthII_creation/SYNTH5/Synth_Balanced_doublets/Barcodes.new.tsv"
#nThreads=10
#outdir=None
#countpath="/hpcnfs/scratch/temporary/Dav_vc/20.4_DemultiMatteo/scRNAseq/demultiplexing/organoidMultiplexing-ensembl_human_93/davide/SYNTH5/scansnp/Counts.pkl"
#rawPath=None
#filteredPath=None



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

if not os.path.exists(outdir):
	os.mkdir(outdir)



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

#sys.path.append('/hpcnfs/scratch/temporary/Dav_vc/20.4_DemultiMatteo/SCanSNP_V1_ssh_EmptyFit/SCanSNP/SCanSNP/')

from VCFUtils import *
from Wrappers import *
from RawBCMatrix_Utils import *
from GenUtils import *

if __name__ == "__main__":
	#Creation of differen loci subsets
	cleanLoci = LociPreClean(vcf)
	MildcleanLoci = LociPreClean_milds(vcf)
	CleanSingularLoci,SingularLoci_Alt,SingularLoci_Ref = SingularLociSCan(vcf,cleanLoci)
	
	FullDrops = pd.read_csv(barcodesFILE, header=None, names=["b"])["b"].tolist()
	
	#If raw counts cellranger matrix is provided launch RawBCMatrix module
	if rawPath is not None:
		FullDropsKNNseries,EmptyBarcodeList = UnfilteredAdataAdata(filteredPath, rawPath, outdir, nHKgenes=50, raw_to_filt_ratio=1)
		barcodeList=FullDrops+EmptyBarcodeList
	else:
		barcodeList = FullDrops
		FullDropsKNNseries = None
	
	
	
	#Genotypes map creation
	GenotypesDF = creategenotypeDF(vcf)
	
	
	if mode == "matrixgen":
		'''
		Performing only count
		'''
		Counts = CountsMatrices(CleanSingularLoci, cleanLoci,
			MildcleanLoci, GenotypesDF,
			barcodeList, vcf,
			nThreads, bamFile)
			
		if rawPath is not None:
			#Save ReadCounts with emptyDrops in pickle
			with open(outdir+ '/Counts.withEmpty.pkl', 'wb') as f:
				pickle.dump(Counts, f, pickle.HIGHEST_PROTOCOL)
			#Save ReadCounts without emptyDrops in pickle
			with open(outdir+ '/Counts.pkl', 'wb') as f:
				pickle.dump(Counts.slice(barcodeList = FullDrops)[0], f, pickle.HIGHEST_PROTOCOL)
		else:
			#Save ReadCounts without emptyDrops in pickle
			with open(outdir+ '/Counts.pkl', 'wb') as f:
				pickle.dump(Counts, f, pickle.HIGHEST_PROTOCOL)
	
	elif mode == "skipcount":
		'''
		Performing only deconvolution based on provided counts
		'''
		#Load existing counts
		with open(str(countpath), 'rb') as f:
			Counts = pickle.load(f)
		#Deconvolution
		Cell_IDs = deconvolution(Counts, vcf, GenotypesDF,outdir,FullDrops, FullDropsKNNseries)
		Cell_IDs.to_csv(outdir + "/Cell_IDs.tsv", sep = "\t", header = True, index = True, index_label = "barcode")
		
		
	else:
		'''
		Perform both count and ceconvolution
		'''
		#counting
		Counts = CountsMatrices(CleanSingularLoci, cleanLoci,
			MildcleanLoci, GenotypesDF,
			barcodeList, vcf,
			nThreads, bamFile)
		
		#Save ReadCounts in pickle
		if rawPath is not None:
			#Save ReadCounts with emptyDrops in pickle
			with open(outdir+ '/Counts.withEmpty.pkl', 'wb') as f:
				pickle.dump(Counts, f, pickle.HIGHEST_PROTOCOL)
			#Save ReadCounts without emptyDrops in pickle
			with open(outdir+ '/Counts.pkl', 'wb') as f:
				pickle.dump(Counts.slice(barcodeList = FullDrops)[0], f, pickle.HIGHEST_PROTOCOL)
		else:
			#Save ReadCounts without emptyDrops in pickle
			with open(outdir+ '/Counts.pkl', 'wb') as f:
				pickle.dump(Counts, f, pickle.HIGHEST_PROTOCOL)
		
		#Deconvolution
		Cell_IDs = deconvolution(Counts, vcf, GenotypesDF,outdir,FullDrops, FullDropsKNNseries)
		Cell_IDs.to_csv(outdir + "/Cell_IDs.tsv", sep = "\t", header = True, index = True, index_label = "barcode")
