
import argparse
import os
#os.environ['R_HOME'] = '/home/davide.castaldi/miniconda3/envs/SCanSNP_devel/lib/R'

parser = argparse.ArgumentParser(description='Process some integers.')


parser.add_argument('--vcf', dest='vcf', help='joint vcf with all mixed genotypes SNPs',required=True,type=str)
parser.add_argument('--barcodes', dest='barcodes', help='for matrixgen mode: FILTERED barcodes file in standard 10x .tsv format',required=True,type=str)


#conditionally optionals
parser.add_argument('--bam', dest='bam', help='for matrixgen mode: BAM file of mixed genotypes',type=str)
parser.add_argument('--counts', dest='countpath',type=str,help='pre-calculated anndata.h5ad with var counts:  adata.layers["RefReads"] and adata.layers["AltReads"]')
parser.add_argument('--filtered_matrix_path', dest='filtered_matrix', help='path/to/cellranger/filtered/matrix/',default=None,type=str)


#Optionals
parser.add_argument('--threads', dest='nthreads', help='threads to be used',default=10,type=str)
parser.add_argument('--raw_matrix_path', dest='raw_matrix', help='path/to/cellranger/unfiltered/matrix/',default=None,type=str)
parser.add_argument('--umitag', dest='umitag', help='tag in the bam file for UMI',default="UB",type=str)
parser.add_argument('--celltag', dest='barcodetag', help='tag in the bam file for cells',default="CB",type=str)
parser.add_argument('--segmentation', dest='segmentation', help='tsv with barcode-nuclei information',default=None,type=str)
parser.add_argument('--outdir', dest='outdir', help='use',default="currentwd",type=str)
#parser.add_argument('--flagLowQual', dest='LowQual', help='Force-fit GMM to flag lowquality droplets',default=False,type=bool)
parser.add_argument('--mode', dest='mode', default="deconvolution",
	required=False,
	choices=["matrixgen","deconvolution","skipcount","pileup"],
	type=str,
	help='creating pileup matrix with ref and alt supporting reads')

parser.add_argument('--platform', dest='platform', default="chromium",
	required=False,
	choices=["chromium","visium"],
	type=str,
	help='select the library platform')

args = parser.parse_args()


if args.mode == "matrixgen":
	if args.bam is None or args.barcodes is None:
		parser.error('please provide --bam and --barcodes for matrixgen mode')

if args.mode == "skipcount":
	if args.countpath is None:
		parser.error('please provide path to anndata.h5ad using --counts for skipcount mode')

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
barcodetag=args.barcodetag
umitag=args.umitag
platform=args.platform
segmentation=args.segmentation


#mode="deconvolution"
#platform="visium"
#bamFile="/group/testa/Users/davide.castaldi/lymph_node_lymphoma_14k_atac_possorted_bam.rmDup.bam"
#vcf="/home/davide.castaldi/projects/scMultiomeData/atac_vc/parsed.VCF"
#barcodesFILE="/home/davide.castaldi/projects/scMultiomeData/filtered_feature_bc_matrix/barcodes.tsv.gz"
#nThreads=20
#outdir="/home/davide.castaldi/projects/scMultiomeData/atac_vc/"
#countpath="/hpcnfs/scratch/temporary/Dav_vc/20.4_DemultiMatteo/scRNAseq/demultiplexing/organoidMultiplexing-ensembl_human_93/davide/SYNTH5/scansnp/Counts.pkl"
#rawPath=None
#filteredPath=None

if platform == "visium" and segmentation is None:
	print("No segmentation provided. SCanSNP will only output information about First and Second IDs per barcode")

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
import anndata as ad
import scanpy as sc
import os


#sys.path.append(os.getcwd())
#sys.path.append('/home/davide.castaldi/git/SCanSNP/SCanSNP')


from VCFUtils_copy import *
from Wrappers import *
from RawBCMatrix_Utils import *
from GenUtils import *
from VCFUtils_copy import *




#Creation of differen loci subsets
if mode != "pileup":
	VariantsFile = VCFdata(vcf)
	#cleanLoci = LociPreClean(vcf)
	#CleanSingularLoci,SingularLoci_Alt,SingularLoci_Ref = SingularLociSCan(vcf,cleanLoci)

FullDrops = pd.read_csv(barcodesFILE, header=None, names=["b"])["b"].astype("string").tolist()

#If raw counts cellranger matrix is provided launch RawBCMatrix module
if rawPath is not None:
	FullDropsKNNseries,EmptyBarcodeList = UnfilteredAdataAdata(filteredPath, rawPath, outdir, nHKgenes=50, raw_to_filt_ratio=1)
	barcodeList=FullDrops+EmptyBarcodeList
else:
	barcodeList = FullDrops
	FullDropsKNNseries = None


GenotypesDF = VariantsFile.creategenotypeDF()

GenotypeChunkList,GenotypeChunkIndexesList = FlattenDict(ChunkMaker(GenotypesDF, nThreads))



locilistChunkX = [GenotypeChunkList[i] for i in GenotypeChunkIndexesList[0]]
locilistChunkX = pd.DataFrame(locilistChunkX, columns = ["chr","start","ref","alt"])
locilistChunkX['end'] = locilistChunkX['start']
locilistChunkX = locilistChunkX[["chr","start","end", "ref","alt"]]

locilistChunkX.to_csv('/tmp/temp.loci.tsv', sep='\t', header=False, index=False )

print("hellp")



