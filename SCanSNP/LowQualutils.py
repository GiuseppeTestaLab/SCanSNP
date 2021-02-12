#!/usr/bin/env python

import itertools
from VCFUtils import *
from sklearn.linear_model import LogisticRegression
from DBLsutils import *
from GenUtils import *

def LocusSpecNoisCalc(DropToDBLDict,Doublet,SingularLociScoreDF,GenotypesDF,vcf,HomozlociList,RefHetLoci,AltHetLoci):
	comp1 = Doublet[0]
	comp2 = Doublet[1]
	#Step 1 work on HomozLoci
	#Step 1 work on HomozLoci
	## - Calculate normalized Noise
	NormalizedNoise=SingularLociScoreDF.loc[DropToDBLDict[Doublet], ["".join([component,"_Norm"]) for component in ExtractSamples(vcf) if component not in [comp1,comp2]]]
	TotalNoise=NormalizedNoise.sum(axis = 1)
	## 1.1 Calculate normalized Ref Allelic contribution of noise
	RefAl=GenotypesDF.loc[HomozlociList, ["".join([component,"_RefAl"]) for component in ExtractSamples(vcf)]]
	RelRefAlAboundance=RefAl.divide(RefAl.sum(axis = 1), axis = 0)
	Noise_RefAl_Ratio=RelRefAlAboundance[["".join([component,"_RefAl"]) for component in ExtractSamples(vcf) if component not in [comp1,comp2]]].sum(axis = 1)
	LocuS_Ref_Noise=Noise_RefAl_Ratio.to_frame().dot(TotalNoise.to_frame().T)
	CleanRefSignal=1-LocuS_Ref_Noise
	## 1.2 Calculate normalized Alt Allelic contribution of noise
	AltAl=GenotypesDF.loc[HomozlociList, ["".join([component,"_AltAl"]) for component in ExtractSamples(vcf)]]
	RelAltAlAboundance=AltAl.divide(AltAl.sum(axis = 1), axis = 0)
	Noise_AltAl_Ratio=RelAltAlAboundance[["".join([component,"_AltAl"]) for component in ExtractSamples(vcf) if component not in [comp1,comp2]]].sum(axis = 1)
	LocuS_Alt_Noise=Noise_AltAl_Ratio.to_frame().dot(TotalNoise.to_frame().T)
	CleanAltSignal=1-LocuS_Alt_Noise
	#Step 2 work on HetLoci
	#Step 2 work on HetLoci
	## 2.1 Calculate normalized Ref Allelic contribution of noise
	Ref_HET_Al=GenotypesDF.loc[RefHetLoci, ["".join([component,"_RefAl"]) for component in ExtractSamples(vcf)]]
	RelRef_HET_AlAboundance=Ref_HET_Al.divide(Ref_HET_Al.sum(axis = 1), axis = 0)
	Noise_Ref_HET_Al_Ratio=RelRef_HET_AlAboundance[["".join([component,"_RefAl"]) for component in ExtractSamples(vcf) if component not in [comp1,comp2]]].sum(axis = 1)
	LocuS_Ref_HET_Noise=Noise_Ref_HET_Al_Ratio.to_frame().dot(TotalNoise.to_frame().T)
	CleanRef_HET_Signal=1-LocuS_Ref_HET_Noise
	## 2.2 Calculate normalized Alt Allelic contribution of noise
	Alt_HET_Al=GenotypesDF.loc[AltHetLoci, ["".join([component,"_AltAl"]) for component in ExtractSamples(vcf)]]
	RelAlt_HET_AlAboundance=Alt_HET_Al.divide(Alt_HET_Al.sum(axis = 1), axis = 0)
	Noise_Alt_HET_Al_Ratio=RelAlt_HET_AlAboundance[["".join([component,"_AltAl"]) for component in ExtractSamples(vcf) if component not in [comp1,comp2]]].sum(axis = 1)
	LocuS_Alt_HET_Noise=Noise_Alt_HET_Al_Ratio.to_frame().dot(TotalNoise.to_frame().T)
	CleanAlt_HET_Signal=1-LocuS_Alt_HET_Noise
	return CleanRefSignal,CleanAltSignal,CleanRef_HET_Signal,CleanAlt_HET_Signal



def LowQualScore(DropToDBLDict,SparseD,DBLSpecificSingularLociHomoz,DBLSpecificSingularLociHeteroRef,DBLSpecificSingularLociHeteroAlt,vcf, GenotypesDF,SingularLociScoreDF):
	DBLInfo = pd.DataFrame( columns = ["DBL_FirstID_Score","DBL_SecondID_Score"])
	for Doublet in ExtractOrderedDoublets(vcf):
		
		comp1 = Doublet[0]
		comp2 = Doublet[1]
		
		#Slice for relevant Loci/Barcodes
		#Homoz Loci
		SRef_HOMOZ,SAlt_HOMOZ,GenotypesDF_sliced_HOMOZ=MultiSlice(SparseD, GenotypesDF, DBLSpecificSingularLociHomoz[Doublet], DropToDBLDict[Doublet])
		#Heteroz Loci Ref
		SRef_RefHet,SAlt_RefHet,GenotypesDF_sliced_RefHet=MultiSlice(SparseD, GenotypesDF, DBLSpecificSingularLociHeteroRef[Doublet], DropToDBLDict[Doublet])
		#Heteroz Loci Alt
		SRef_AltHet,SAlt_AltHet,GenotypesDF_sliced_AltHet=MultiSlice(SparseD, GenotypesDF, DBLSpecificSingularLociHeteroAlt[Doublet], DropToDBLDict[Doublet])
		# Locus-specific Noise calculation
		CleanRefSignal,CleanAltSignal,CleanRef_HET_Signal,CleanAlt_HET_Signal=LocusSpecNoisCalc(DropToDBLDict,
			Doublet,SingularLociScoreDF,GenotypesDF,vcf,
			DBLSpecificSingularLociHomoz[Doublet],
			DBLSpecificSingularLociHeteroRef[Doublet],
			DBLSpecificSingularLociHeteroAlt[Doublet])
			
		GenotypesDF_sliced_HOMOZ=GenotypesDF_sliced_HOMOZ.replace(1,2)
		GenotypesDF_sliced_RefHet=GenotypesDF_sliced_RefHet.replace(1,2)
		GenotypesDF_sliced_AltHet=GenotypesDF_sliced_AltHet.replace(1,2)
		##First thing use HomoLoci information
		##First thing use HomoLoci information
		##First thing use HomoLoci information
		# Cleaning Signals from Locus Nloise
		RefReads_HOMOZ=SRef_HOMOZ.multiply(CleanRefSignal)
		AltReads_HOMOZ=SAlt_HOMOZ.multiply(CleanAltSignal)
		#For each component of Putative DBl we calc the score
		HomozRefScoresComp1 = RefReads_HOMOZ.T.multiply(GenotypesDF_sliced_HOMOZ[comp1+"_RefAl"].to_numpy()).T
		HomozRefScoresComp2 = RefReads_HOMOZ.T.multiply(GenotypesDF_sliced_HOMOZ[comp2+"_RefAl"].to_numpy()).T
		HomozAltScoresComp1 = AltReads_HOMOZ.T.multiply(GenotypesDF_sliced_HOMOZ[comp1+"_AltAl"].to_numpy()).T
		HomozAltScoresComp2 = AltReads_HOMOZ.T.multiply(GenotypesDF_sliced_HOMOZ[comp2+"_AltAl"].to_numpy()).T
		HomozComp1Total = HomozRefScoresComp1+HomozAltScoresComp1
		HomozComp2Total = HomozRefScoresComp2+HomozAltScoresComp2
		
		#Now info about Singular Ref loci
		#Now info about Singular Ref loci
		#Now info about Singular Ref loci
		# Cleaning Signals from Locus Nloise
		RefReads_HETEROZ=SRef_RefHet.multiply(CleanRef_HET_Signal)
		HeterozRefScoreComp1 = RefReads_HETEROZ.T.multiply(GenotypesDF_sliced_RefHet[comp1+"_RefAl"].to_numpy()).T
		HeterozRefScoreComp2 = RefReads_HETEROZ.T.multiply(GenotypesDF_sliced_RefHet[comp2+"_RefAl"].to_numpy()).T
		
		#Now those of Singular Alt loci
		#Now those of Singular Alt loci
		#Now those of Singular Alt loci
		AltReads_HETEROZ=SAlt_AltHet.multiply(CleanAlt_HET_Signal)
		HeterozAltScoreComp1 = AltReads_HETEROZ.T.multiply(GenotypesDF_sliced_AltHet[comp1+"_AltAl"].to_numpy()).T
		HeterozAltScoreComp2 = AltReads_HETEROZ.T.multiply(GenotypesDF_sliced_AltHet[comp2+"_AltAl"].to_numpy()).T
		
		### COMPUTE TOTAL SCORES
		### COMPUTE TOTAL SCORES
		### COMPUTE TOTAL SCORES
		Comp1CumulativeScore = HomozComp1Total.sum(axis =0)+HeterozRefScoreComp1.sum(axis =0)+HeterozAltScoreComp1.sum(axis =0)
		Comp2CumulativeScore = HomozComp2Total.sum(axis =0)+HeterozRefScoreComp2.sum(axis =0)+HeterozAltScoreComp2.sum(axis =0)
		
		DBLstatus=pd.concat([pd.DataFrame(Comp1CumulativeScore.T, index = DropToDBLDict[Doublet], columns = ["DBL_FirstID_Score"]),
			pd.DataFrame(Comp2CumulativeScore.T, index = DropToDBLDict[Doublet], columns = ["DBL_SecondID_Score"])], axis = 1)
			
		DBLInfo = pd.concat([DBLInfo,DBLstatus.round(2)], axis = 0)
		
	return DBLInfo
