#!/usr/bin/env python

from SCanSNP.VCFUtils import *
from SCanSNP.DBLsutils import *
#from SCanSNP.Pileup import *
from scipy.sparse import csr_matrix
import copy
import anndata as ad
import numpy as np

#from Wrappers import *


# def MultiSlice(Counts, GenotypesDF, LocusList, BarcodeList=None):
# 	'''
# 	Simultaneus slicing of Sparse counts and genotypes DF according do provided locusList and BarcodesList (if given)
# 	'''
# 	#Locus-based Mask
# 	Counts_Locusmask = Counts.loci.isin(LocusList)
# 	#Barcodes-based Mask
# 	Counts_BarcodesMask = Counts["Barcode"].isin(BarcodeList if BarcodeList != None  else Counts["Barcode"])
# 	#Slicing RefCounts
# 	#Slicing RefCounts
# 	SRef = Counts["sparse_Ref"][Counts_Locusmask]
# 	SRef = SRef[:,Counts_BarcodesMask]
# 	#Slicing AltCounts
# 	#Slicing AltCounts
# 	SAlt = Counts["sparse_Alt"][Counts_Locusmask]
# 	SAlt = SAlt[:,Counts_BarcodesMask]
# 	#Slicing GenotypesDF
# 	GenotypesDF_sliced = GenotypesDF.loc[Counts.loci[Counts_Locusmask]]
# 	return SRef,SAlt,GenotypesDF_sliced


def SingularLociCNTR( SingularLoci_Alt, SingularLoci_Ref, Counts, barcodeList, GenotypesDF,vcf):
	'''
	We compute Genotype contributions for i loci across n cells
	-Giving same weight at Het and Homoz singular loci to account for General Id contribution
	'''
	#RefQuery
	SingularLoci_Ref_Counts,RefSingularGenotype  = Counts.slice(lociList=SingularLoci_Ref, Genotypes=GenotypesDF)
	RefSingularGenotype =  RefSingularGenotype[[ID+'_RefAl' for ID in list(ExtractSamples(vcf))]].replace(1,2)
	
	#Alt Query
	SingularLoci_Alt_Counts,AltSingularGenotype  = Counts.slice(lociList=SingularLoci_Alt, Genotypes=GenotypesDF)
	AltSingularGenotype = AltSingularGenotype[[ID+'_AltAl' for ID in list(ExtractSamples(vcf))]].replace(1,2)
	
	
	#Init data strucutres
	CNTRPerBArcodeAlt = {}
	CNTRPerBArcodeRef = {}
	TotalCNTR = pd.DataFrame( columns = list(ExtractSamples(vcf)), index = Counts.barcodes).fillna(0.0)
	IDsToQueryAlt = [ID+'_AltAl' for ID in list(ExtractSamples(vcf))]
	IDsToQueryRef = [ID+'_RefAl' for ID in list(ExtractSamples(vcf))]
	#PileUp for Alt informative loci
	for ID in IDsToQueryAlt:
		AltGenotypeID = AltSingularGenotype[ID].to_numpy()
		CNTRPerBArcodeAlt["_".join(ID.split("_")[:-1])] = np.array(SingularLoci_Alt_Counts.sparseAlt.multiply(AltGenotypeID).sum(axis = 1))
	#PileUp for Ref informative loci
	for ID in IDsToQueryRef:
		RefGenotypeID = RefSingularGenotype[ID]
		CNTRPerBArcodeRef["_".join(ID.split("_")[:-1])] = np.array(SingularLoci_Ref_Counts.sparseRef.multiply(RefGenotypeID).sum(axis = 1))
	#Sum both pileups into DF
	for ID in CNTRPerBArcodeRef.keys():
		TotalCNTR[ID] =  CNTRPerBArcodeRef[ID]+CNTRPerBArcodeAlt[ID]
	NormCNTR=TotalCNTR.divide(TotalCNTR.sum(axis = 1), axis =0).round(2)
	NormCNTR.columns = [Samp + "_Norm" for Samp in TotalCNTR.columns]
	TotalCNTR = pd.concat([TotalCNTR,NormCNTR], axis = 1)
	return TotalCNTR



class CountData:
	
	def __init__(self, sparseRef, sparseAlt, loci, barcodes):
		self.sparseRef = sparseRef.T
		self.sparseAlt = sparseAlt.T
		self.loci = loci
		self.barcodes = barcodes
	
	
	def copy(self):
		return copy.copy(self)
	
	def slice(self, lociList=None, barcodeList=None, Genotypes=None):
		'''
		Simultaneus slicing of Sparse counts and genotypes DF according do provided locusList and BarcodesList (if given)
		'''
		slicedCounts = self.copy()
		
		#Locus-based Mask
		Counts_Locusmask = self.loci.isin(lociList if lociList is not None else self.loci)
		#Barcodes-based Mask
		Counts_BarcodesMask = self.barcodes.isin(barcodeList if barcodeList is not None else self.barcodes)
		#Slicing RefCounts
		#Slicing RefCounts
		SRef = self.sparseRef[:,Counts_Locusmask]
		slicedCounts.sparseRef = SRef[Counts_BarcodesMask]
		#Slicing AltCounts
		#Slicing AltCounts
		SAlt = self.sparseAlt[:,Counts_Locusmask]
		slicedCounts.sparseAlt = SAlt[Counts_BarcodesMask]
		
		#Slice Barcodes
		slicedCounts.barcodes = self.barcodes[Counts_BarcodesMask]
		
		#Slice loci
		slicedCounts.loci = self.loci[Counts_Locusmask]
		
		if Genotypes is not None:
			SGenotypes = Genotypes.loc[slicedCounts.loci]
		else:
			SGenotypes =None
		
		return slicedCounts, SGenotypes
	
	def write_h5ad(self, writepath=None):
		varAdata = ad.AnnData(X=self.sparseRef, dtype=self.sparseRef.dtype)
		varAdata.obs_names = self.barcodes.tolist()
		varAdata.var_names = self.loci.tolist()
		varAdata.layers["RefReads"] = self.sparseRef
		varAdata.layers["AltReads"] = self.sparseAlt
		varAdata.write_h5ad(writepath)
		
		return None

