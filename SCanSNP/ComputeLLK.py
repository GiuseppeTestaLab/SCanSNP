#!/usr/bin/env python

#Likelihood calculaiton model

import pysam
from VCFUtils import *


def ComputeLikelihood(vcf,GenotypesDF,countDf,barcodeList,DiffOnlyIndex):
	'''
	The ComputeLikelihood function calculates and sum up Altscores and Refscores of the possible mixed genotypes for every locus barcode-wise
	Altscore: (#Alt reads)*(#Alt alleles in the genotype)/(#Total Alt alleles across all mixed genotypes)
	Refscore: (#Ref reads)*(#Ref alleles in the genotype)/(#Total Ref alleles across all mixed genotypes)
	PerGenotypeLLK = Altscore + Refscore
	IGNORING DOUBLETS!!
	'''
	GenotypesDF_DiffOnly = GenotypesDF.loc[DiffOnlyIndex]
	countDf_DiffOnly = countDf.loc[DiffOnlyIndex]
	totRefAlleles = GenotypesDF_DiffOnly[[ sample+"_RefAl" for sample in ExtractSamples(vcf)]].sum(axis =1)
	totAltAlleles = GenotypesDF_DiffOnly[[ sample+"_AltAl" for sample in ExtractSamples(vcf)]].sum(axis =1)
	ScorePerBArcode = pd.DataFrame( columns = list(ExtractSamples(vcf)), index = [ cell for cell in barcodeList])
	for geno in list(ExtractSamples(vcf)):
		RefPerGeno=(countDf_DiffOnly[[cell+"_RefReads" for cell  in barcodeList]].multiply(GenotypesDF_DiffOnly[geno+"_RefAl"]*0.5, axis = 0)).divide(totRefAlleles, axis = 0)
		AltPerGeno=(countDf_DiffOnly[[cell+"_AltReads" for cell  in barcodeList]].multiply(GenotypesDF_DiffOnly[geno+"_AltAl"]*0.5, axis = 0)).divide(totAltAlleles, axis = 0)
		RefPerGeno.columns = barcodeList
		AltPerGeno.columns = barcodeList
		ScorePerBArcode[geno]=RefPerGeno.add(AltPerGeno, axis = 0 ,fill_value = 0.0).sum(axis = 0)
	return ScorePerBArcode


def ComputeSecondIDLikelihood( SingularLoci_Alt, SingularLoci_Ref, countDf, barcodeList):
	AltSingularGenotype = GenotypesDF.loc[SingularLoci_Alt,[ID+'_AltAl' for ID in list(ExtractSamples(vcf))]].replace(2,1)
	RefSingularGenotype =  GenotypesDF.loc[SingularLoci_Ref,[ID+'_RefAl' for ID in list(ExtractSamples(vcf))]].replace(2,1)
	AltCuntsSS= countDf.loc[SingularLoci_Alt,[barcode+'_AltReads' for barcode in barcodeList]]
	RefCuntsSS= countDf.loc[SingularLoci_Ref,[barcode+'_RefReads' for barcode in barcodeList]]
	ScorePerBArcodeAlt = pd.DataFrame( columns = list(ExtractSamples(vcf)), index = barcodeList).fillna(0.0)
	ScorePerBArcodeRef = pd.DataFrame( columns = list(ExtractSamples(vcf)), index = barcodeList).fillna(0.0)
	IDsToQueryAlt = [ID+'_AltAl' for ID in list(ExtractSamples(vcf))]
	IDsToQueryRef = [ID+'_RefAl' for ID in list(ExtractSamples(vcf))]
	#PileUp for Alt informative loci
	for ID in IDsToQueryAlt:
		AltGenotypeID = AltSingularGenotype[ID]
		AltCuntsSS.columns = barcodeList
		ScorePerBArcodeAlt["_".join(ID.split("_")[:-1])] = AltCuntsSS.multiply(AltGenotypeID, axis = 0).sum()
	#PileUp for Ref informative loci
	for ID in IDsToQueryRef:
		RefGenotypeID = RefSingularGenotype[ID]
		RefCuntsSS.columns = barcodeList
		ScorePerBArcodeRef["_".join(ID.split("_")[:-1])] = RefCuntsSS.multiply(RefGenotypeID, axis = 0).sum()
	SingularLociScores = ScorePerBArcodeRef.add(ScorePerBArcodeAlt,fill_value = 0.0)
	return SingularLociScores
