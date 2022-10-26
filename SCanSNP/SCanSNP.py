
import argparse
import os


parser = argparse.ArgumentParser(description='Process some integers.')


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
parser.add_argument('--mito_contig', dest='mitoContig', help='use',type=str)
parser.add_argument('--vcf', dest='vcf', help='joint vcf with all mixed genotypes SNPs',required=False,type=str)


#parser.add_argument('--flagLowQual', dest='LowQual', help='Force-fit GMM to flag lowquality droplets',default=False,type=bool)
parser.add_argument('--mode', dest='mode', default="deconvolution",
	required=False,
	choices=["matrixgen","deconvolution","skipcount","pileup","mito"],
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

if args.mode != "mito":
	if args.vcf is None:
		parser.error('please provide --vcf reference for deconvolution')

if args.mode == "mito":
	if args.mitoContig is None:
		parser.error('please provide --mito_contig')

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
mitoContig=args.mitoContig


#mode="mito"
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

#sys.path.append('/home/davide.castaldi/git/SCanSNP/SCanSNP')

from SCanSNP.VCFUtils import *
from SCanSNP.Wrappers import *
from SCanSNP.RawBCMatrix_Utils import *
from SCanSNP.GenUtils import *
from SCanSNP.mitoUtils import *


def main():
	#Creation of differen loci subsets
	if mode not in ["pileup","mito"]:
		cleanLoci = LociPreClean(vcf)
		MildcleanLoci = LociPreClean_milds(vcf)
		CleanSingularLoci,SingularLoci_Alt,SingularLoci_Ref = SingularLociSCan(vcf,cleanLoci)
		GenotypesDF = creategenotypeDF(vcf)
	
	FullDrops = pd.read_csv(barcodesFILE, header=None, names=["b"])["b"].astype("string").tolist()
	
	#If raw counts cellranger matrix is provided launch RawBCMatrix module
	if rawPath is not None:
		FullDropsKNNseries,EmptyBarcodeList = UnfilteredAdataAdata(filteredPath, rawPath, outdir, nHKgenes=50, raw_to_filt_ratio=1)
		barcodeList=FullDrops+EmptyBarcodeList
	else:
		barcodeList = FullDrops
		FullDropsKNNseries = None
	
	
	
	
	if mode == "matrixgen":
		'''
		Performing only count
		'''
		Counts = CountsMatrices(CleanSingularLoci, cleanLoci,
			MildcleanLoci, GenotypesDF,
			barcodeList, vcf,
			nThreads, bamFile , barcodetag, umitag, Ref=True)
			
		if rawPath is not None:
			#Save ReadCounts with emptyDrops in Anndata
			Counts.write_h5ad(outdir+'/varAdata.h5ad')	
			#Save ReadCounts without emptyDrops in Anndata
			Counts.slice(barcodeList = FullDrops)[0].write_h5ad(outdir+'/varAdata.h5ad')
		else:
			#Save ReadCounts without emptyDrops in Anndata
			Counts.write_h5ad(outdir+'/varAdata.h5ad')
	
	elif mode == "skipcount":
		'''
		Performing only deconvolution based on provided counts
		'''
		#Load existing counts
		varAdata = sc.read_h5ad(countpath)
		Counts = CountData(varAdata.layers["RefReads"], varAdata.layers["AltReads"], varAdata.var_names, varAdata.obs_names, Ref=True)
		del varAdata
		#Deconvolution
		Cell_IDs = deconvolution(Counts, vcf, GenotypesDF,outdir,FullDrops, FullDropsKNNseries, platform, segmentation)
		Cell_IDs.to_csv(outdir + "/Cell_IDs.tsv", sep = "\t", header = True, index = True, index_label = "barcode")
		
		
	elif mode == "deconvolution":
		'''
		Perform both count and ceconvolution
		'''
		#counting
		Counts = CountsMatrices(CleanSingularLoci, cleanLoci,
			MildcleanLoci, GenotypesDF,
			barcodeList, vcf,
			nThreads, bamFile , barcodetag, umitag)
		
		#Save ReadCounts in Anndata
		if rawPath is not None:
			#Save ReadCounts with emptyDrops in Anndata
			Counts.write_h5ad(outdir+'/varAdata.h5ad')	
			#Save ReadCounts without emptyDrops in Anndata
			Counts.slice(barcodeList = FullDrops)[0].write_h5ad(outdir+'/varAdata.h5ad')
		else:
			#Save ReadCounts without emptyDrops in Anndata
			Counts.write_h5ad(outdir+'/varAdata.h5ad')
		
		#Deconvolution
		Cell_IDs = deconvolution(Counts, vcf, GenotypesDF,outdir,FullDrops, FullDropsKNNseries, platform, segmentation)
		Cell_IDs.to_csv(outdir + "/Cell_IDs.tsv", sep = "\t", header = True, index = True, index_label = "barcode")

	elif mode == "pileup":
		'''
		Performing only count on provided loci list, this mode assumes single-sample VCF file.
		'''
		
		GenotypesDF =  pd.read_csv(vcf, sep ="\t", header=None, names=["CHROM","POS","REF","ALT"])
		GenotypesDF["CHROM"] = GenotypesDF["CHROM"].astype(str)
		GenotypesDF.index = GenotypesDF["CHROM"]+"_"+GenotypesDF["POS"].astype(str)
		MildcleanLoci = GenotypesDF[(GenotypesDF["REF"].str.len() == 1) & (GenotypesDF["ALT"].str.len() == 1)].index.tolist()
		
		
		Counts = CountsPileup(MildcleanLoci, GenotypesDF,barcodeList, vcf,nThreads, bamFile, barcodetag, umitag, Ref=True)
			
		if rawPath is not None:
			#Save ReadCounts with emptyDrops in Anndata
			Counts.write_h5ad(outdir+'/varAdata.h5ad')	
			#Save ReadCounts without emptyDrops in Anndata
			Counts.slice(barcodeList = FullDrops)[0].write_h5ad(outdir+'/varAdata.h5ad')
		else:
			#Save ReadCounts without emptyDrops in Anndata
			Counts.write_h5ad(outdir+'/varAdata.h5ad')

	elif mode == "mito":
		'''
		Pileup on all positions of provided (mito) chromosome only
		'''
		
		GenotypesDF = ExtractMitoPositions(bamFile, mitoContig)
		MildcleanLoci = GenotypesDF.index.tolist()
	
	
		Counts = CountsPileup(MildcleanLoci, GenotypesDF,barcodeList,nThreads, bamFile, barcodetag, umitag, Ref=False)
			
		Counts.to_csv(outdir+"/{}.Counts.tsv".format(mitoContig),sep="\t")

if __name__ == "__main__":
	main()
