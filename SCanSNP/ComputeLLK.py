#!/usr/bin/env python

#Likelihood calculaiton model

import pysam
from SCanSNP.VCFUtils import *
from SCanSNP.GenUtils import *

def ComputeLikelihood(vcf,GenotypesDF,Counts,DiffOnlyIndex):
	'''
	The ComputeLikelihood function calculates and sum up Altscores and Refscores of the possible mixed genotypes for every locus barcode-wise
	Altscore: (#Alt reads)*(#Alt alleles in the genotype)/(#Total Alt alleles across all mixed genotypes)
	Refscore: (#Ref reads)*(#Ref alleles in the genotype)/(#Total Ref alleles across all mixed genotypes)
	PerGenotypeLLK = Altscore + Refscore
	IGNORING DOUBLETS!!
	'''
	LLKCounts,LLKGenotypesDF =  Counts.slice(lociList = DiffOnlyIndex, Genotypes = GenotypesDF)
	totRefAlleles = LLKGenotypesDF[[ sample+"_RefAl" for sample in ExtractSamples(vcf)]].sum(axis =1)
	totAltAlleles = LLKGenotypesDF[[ sample+"_AltAl" for sample in ExtractSamples(vcf)]].sum(axis =1)
	ScorePerBArcode = pd.DataFrame( columns = list(ExtractSamples(vcf)), index = Counts.barcodes)
	for geno in list(ExtractSamples(vcf)):
		# in order to account for ploidy
		npRefGenotypes = (LLKGenotypesDF[geno+"_RefAl"]*0.5).to_numpy()
		RefPerGeno=LLKCounts.sparseRef.transpose().multiply(npRefGenotypes).multiply(1/totRefAlleles.to_numpy()).transpose()
		npAltGenotypes = (LLKGenotypesDF[geno+"_AltAl"]*0.5).to_numpy()
		AltPerGeno=LLKCounts.sparseAlt.transpose().multiply(npAltGenotypes).multiply(1/totAltAlleles.to_numpy()).transpose()
		ScorePerBArcode[geno]=np.array((AltPerGeno+RefPerGeno).transpose().sum(axis = 1))
	return ScorePerBArcode


def ComputeSecondIDLikelihood(SingularLoci_Alt, SingularLoci_Ref, countDf, barcodeList):
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
