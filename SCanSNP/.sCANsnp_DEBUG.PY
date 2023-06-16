
import argparse
import os


mode="skipcount"
platform="chromium"
bamFile="/group/testa/Common/multiome/align/AssembloidsSeqRun2201multiome/amaia.tintori/R11/outs/gex_possorted_bam.bam"
vcf="/group/testa/Project/AssembloidsDV/scRNASeq/Data/SeqRun2022Multiome/UpStream/WGS_aggregatedVariantCall/test/final.vcf"
barcodesFILE="/group/testa/Common/multiome/align/AssembloidsSeqRun2201multiome/amaia.tintori/R11/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
nThreads=4
umitag="UB"
outdir="/home/davide.castaldi/amaia/amaiaTest"
barcodetag="CB"
countpath="/home/davide.castaldi/amaia/amaiaTest/varAdata.h5ad"
rawPath=None
filteredPath=None
segmentation=None

if platform == "visium" and segmentation is None:
	print("No segmentation provided. SCanSNP will only output information about First and Second IDs per barcode")




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

def main():
	SampleParseDict = ExtractSamples(vcf, return_dict=True)
	#Creation of differen loci subsets
	if mode != "pileup":
		cleanLoci = LociPreClean(vcf)
		MildcleanLoci = LociPreClean_milds(vcf)
		CleanSingularLoci,SingularLoci_Alt,SingularLoci_Ref = SingularLociSCan(vcf,cleanLoci)
	
	FullDrops = pd.read_csv(barcodesFILE, header=None, names=["b"])["b"].astype("string").tolist()
	
	#If raw counts cellranger matrix is provided launch RawBCMatrix module
	if rawPath is not None:
		FullDropsKNNseries,EmptyBarcodeList = UnfilteredAdataAdata(filteredPath, rawPath, outdir, nHKgenes=50, raw_to_filt_ratio=1)
		barcodeList=FullDrops+EmptyBarcodeList
	else:
		barcodeList = FullDrops
		FullDropsKNNseries = None
	
	
	
	#Genotypes map creation
	if mode != "pileup":
		GenotypesDF = creategenotypeDF(vcf)
	
	
	
	if mode == "matrixgen":
		'''
		Performing only count
		'''
		Counts = CountsMatrices(CleanSingularLoci, cleanLoci,
			MildcleanLoci, GenotypesDF,
			barcodeList, vcf,
			nThreads, bamFile , barcodetag, umitag)
			
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
		Counts = CountData(varAdata.layers["RefReads"].T, varAdata.layers["AltReads"].T, varAdata.var_names, varAdata.obs_names)
		del varAdata
		#Deconvolution
		Cell_IDs = deconvolution(Counts, vcf, GenotypesDF,outdir,FullDrops, FullDropsKNNseries, platform, segmentation)
		Cell_IDs.replace(SampleParseDict).to_csv(outdir + "/Cell_IDs.tsv", sep = "\t", header = True, index = True, index_label = "barcode")
		
		
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
		Cell_IDs.replace(SampleParseDict).to_csv(outdir + "/Cell_IDs.tsv", sep = "\t", header = True, index = True, index_label = "barcode")

	elif mode == "pileup":
		'''
		Performing only count on provided loci list, this mode assumes single-sample VCF file.
		'''
		
		GenotypesDF =  pd.read_csv(vcf, sep ="\t", header=None, names=["CHROM","POS","REF","ALT"])
		GenotypesDF["CHROM"] = GenotypesDF["CHROM"].astype(str)
		GenotypesDF.index = GenotypesDF["CHROM"]+"_"+GenotypesDF["POS"].astype(str)
		MildcleanLoci = GenotypesDF[(GenotypesDF["REF"].str.len() == 1) & (GenotypesDF["ALT"].str.len() == 1)].index.tolist()
		
		
		Counts = CountsPileup(MildcleanLoci, GenotypesDF,barcodeList, vcf,nThreads, bamFile, barcodetag, umitag)
			
		if rawPath is not None:
			#Save ReadCounts with emptyDrops in Anndata
			Counts.write_h5ad(outdir+'/varAdata.h5ad')	
			#Save ReadCounts without emptyDrops in Anndata
			Counts.slice(barcodeList = FullDrops)[0].write_h5ad(outdir+'/varAdata.h5ad')
		else:
			#Save ReadCounts without emptyDrops in Anndata
			Counts.write_h5ad(outdir+'/varAdata.h5ad')



if __name__ == "__main__":
	main()