#!/usr/bin/env python

from VCFUtils import *
from DBLsutils import *
from Pileup import *
from scipy.sparse import csr_matrix
#from Wrappers import *


def MultiSlice(SparseD, GenotypesDF, LocusList, BarcodeList=None):
	'''
	Simultaneus slicing of Sparse counts and genotypes DF according do provided locusList and BarcodesList (if given)
	'''
	#Locus-based Mask
	SparseD_Locusmask = SparseD["Locus"].isin(LocusList)
	#Barcodes-based Mask
	SparseD_BarcodesMask = SparseD["Barcode"].isin(BarcodeList if BarcodeList != None  else SparseD["Barcode"])
	#Slicing RefCounts
	#Slicing RefCounts
	SRef = SparseD["sparse_Ref"][SparseD_Locusmask]
	SRef = SRef[:,SparseD_BarcodesMask]
	#Slicing AltCounts
	#Slicing AltCounts
	SAlt = SparseD["sparse_Alt"][SparseD_Locusmask]
	SAlt = SAlt[:,SparseD_BarcodesMask]
	#Slicing GenotypesDF
	GenotypesDF_sliced = GenotypesDF.loc[SparseD["Locus"][SparseD_Locusmask]]
	return SRef,SAlt,GenotypesDF_sliced


def SingularLociCNTR( SingularLoci_Alt, SingularLoci_Ref, SparseD, barcodeList, GenotypesDF,vcf):
	'''
	We compute Genotype contributions for i loci across n cells
	-Giving same weight at Het and Homoz singular loci to account for General Id contribution
	'''
	#RefQuery
	SparseD_BolMask_Ref = SparseD["Locus"].isin(SingularLoci_Ref)
	RefSingularGenotype =  GenotypesDF.loc[SparseD["Locus"][SparseD_BolMask_Ref],[ID+'_RefAl' for ID in list(ExtractSamples(vcf))]].replace(1,2)
	SRef = SparseD["sparse_Ref"][SparseD_BolMask_Ref]
	#Alt Query
	SparseD_BolMask_Alt = SparseD["Locus"].isin(SingularLoci_Alt)
	AltSingularGenotype = GenotypesDF.loc[SparseD["Locus"][SparseD_BolMask_Alt],[ID+'_AltAl' for ID in list(ExtractSamples(vcf))]].replace(1,2)
	SAlt = SparseD["sparse_Alt"][SparseD_BolMask_Alt]
	#Init data strucutres
	CNTRPerBArcodeAlt = {}
	CNTRPerBArcodeRef = {}
	TotalCNTR = pd.DataFrame( columns = list(ExtractSamples(vcf)), index = SparseD["Barcode"]).fillna(0.0)
	IDsToQueryAlt = [ID+'_AltAl' for ID in list(ExtractSamples(vcf))]
	IDsToQueryRef = [ID+'_RefAl' for ID in list(ExtractSamples(vcf))]
	#PileUp for Alt informative loci
	for ID in IDsToQueryAlt:
		AltGenotypeID = AltSingularGenotype[ID].to_numpy()
		CNTRPerBArcodeAlt["_".join(ID.split("_")[:-1])] = np.array(SAlt.transpose().multiply(AltGenotypeID).sum(axis = 1))
	#PileUp for Ref informative loci
	for ID in IDsToQueryRef:
		RefGenotypeID = RefSingularGenotype[ID]
		CNTRPerBArcodeRef["_".join(ID.split("_")[:-1])] = np.array(SRef.transpose().multiply(RefGenotypeID).sum(axis = 1))
	#Sum both pileups into DF
	for ID in CNTRPerBArcodeRef.keys():
		TotalCNTR[ID] =  CNTRPerBArcodeRef[ID]+CNTRPerBArcodeAlt[ID]
	NormCNTR=TotalCNTR.divide(TotalCNTR.sum(axis = 1), axis =0).round(2)
	NormCNTR.columns = [Samp + "_Norm" for Samp in TotalCNTR.columns]
	TotalCNTR = pd.concat([TotalCNTR,NormCNTR], axis = 1)
	return TotalCNTR
