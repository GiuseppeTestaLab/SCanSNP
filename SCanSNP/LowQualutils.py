import itertools
from VCFUtils import *
from sklearn.linear_model import LogisticRegression


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


# def DBLScore(DropToDBLDict,countDF,DBLSpecificSingularLociHomoz,DBLSpecificSingularLociHeteroRef,DBLSpecificSingularLociHeteroAlt,vcf, GenotypesDF,SingularLociScoreDF):
# 	DBLInfo = pd.DataFrame( columns = ["DBL_InformativeUMIs","DBL_FirstID_Score","DBL_SecondID_Score"])
# 	for Doublet in ExtractOrderedDoublets(vcf):
# 		HomozlociList=DBLSpecificSingularLociHomoz[Doublet]
# 		RefHetLoci=DBLSpecificSingularLociHeteroRef[Doublet]
# 		AltHetLoci=DBLSpecificSingularLociHeteroAlt[Doublet]
# 		AltReads=countDF[[Drop+'_AltReads' for Drop in DropToDBLDict[Doublet]]]
# 		RefReads=countDF[[Drop+'_RefReads' for Drop in DropToDBLDict[Doublet]]]
# 		AltReads.columns = [col.split("_")[0] for col in  AltReads.columns]
# 		RefReads.columns = [col.split("_")[0] for col in  RefReads.columns]
# 		InformativeReads=AltReads.loc[AltHetLoci].sum()+RefReads.loc[RefHetLoci].sum()+RefReads.loc[HomozlociList].sum()+AltReads.loc[HomozlociList].sum()
# 		comp1 = Doublet[0]
# 		comp2 = Doublet[1]
# 		HomozRefSingularGenotype = GenotypesDF.loc[HomozlociList,[comp1+"_RefAl",comp2+"_RefAl"]].replace(1,2)
# 		HomozAltSingularGenotype = GenotypesDF.loc[HomozlociList,[comp1+"_AltAl",comp2+"_AltAl"]].replace(1,2)
# 		HeterozRefSingularGenotype = GenotypesDF.loc[RefHetLoci,[comp1+"_RefAl",comp2+"_RefAl"]].replace(1,2)
# 		HeterozAltSingularGenotype = GenotypesDF.loc[AltHetLoci,[comp1+"_AltAl",comp2+"_AltAl"]].replace(1,2)
# 		CleanRefSignal,CleanAltSignal,CleanRef_HET_Signal,CleanAlt_HET_Signal=LocusSpecNoisCalc(DropToDBLDict,Doublet,SingularLociScoreDF,GenotypesDF,vcf,HomozlociList,RefHetLoci,AltHetLoci)
# 		#First thing use HomoLoci information
# 		#First thing use HomoLoci information
# 		#First thing use HomoLoci information
# 		RefReads_HOMOZ=RefReads.loc[HomozlociList].multiply(CleanRefSignal)
# 		AltReads_HOMOZ=AltReads.loc[HomozlociList].multiply(CleanAltSignal)
# 		HomozRefScoresComp1 = RefReads_HOMOZ.multiply(HomozRefSingularGenotype[comp1+"_RefAl"], axis = 0)
# 		HomozRefScoresComp2 = RefReads_HOMOZ.multiply(HomozRefSingularGenotype[comp2+"_RefAl"], axis = 0)
# 		HomozAltScoresComp1 = AltReads_HOMOZ.multiply(HomozAltSingularGenotype[comp1+"_AltAl"], axis = 0)
# 		HomozAltScoresComp2 = AltReads_HOMOZ.multiply(HomozAltSingularGenotype[comp2+"_AltAl"], axis = 0)
# 		HomozComp1Total = HomozRefScoresComp1.add(HomozAltScoresComp1,fill_value = 0.0)
# 		HomozComp2Total = HomozRefScoresComp2.add(HomozAltScoresComp2,fill_value = 0.0)
# 		#Now info about Singular Ref loci
# 		#Now info about Singular Ref loci
# 		#Now info about Singular Ref loci
# 		RefReads_HETEROZ=RefReads.loc[RefHetLoci].multiply(CleanRef_HET_Signal)
# 		HeterozRefScoreComp1 = RefReads_HETEROZ.multiply(HeterozRefSingularGenotype[comp1+"_RefAl"], axis = 0)
# 		HeterozRefScoreComp2 = RefReads_HETEROZ.multiply(HeterozRefSingularGenotype[comp2+"_RefAl"], axis = 0)
# 		#Now those of Singular Alt loci
# 		#Now those of Singular Alt loci
# 		#Now those of Singular Alt loci
# 		AltReads_HETEROZ=AltReads.loc[AltHetLoci].multiply(CleanAlt_HET_Signal)
# 		HeterozAltScoreComp1 = AltReads_HETEROZ.multiply(HeterozAltSingularGenotype[comp1+"_AltAl"], axis = 0)
# 		HeterozAltScoreComp2 = AltReads_HETEROZ.multiply(HeterozAltSingularGenotype[comp2+"_AltAl"], axis = 0)
# 		### TOTAL SCORES
# 		Comp1CumulativeScore = HomozComp1Total.sum()+HeterozRefScoreComp1.sum()+HeterozAltScoreComp1.sum()
# 		Comp2CumulativeScore = HomozComp2Total.sum()+HeterozRefScoreComp2.sum()+HeterozAltScoreComp2.sum()
# 		TotalScore = Comp1CumulativeScore.add(Comp2CumulativeScore,fill_value = 0.0)
# 		#SNG score will be the ratio between reads assigned to 1st id and reads assigned to second id
# 		#Wrapping together into DF slice
# 		DBLstatus = pd.concat([InformativeReads.to_frame(name = "DBL_InformativeUMIs"),
# 			Comp1CumulativeScore.to_frame(name = "DBL_FirstID_Score"),
# 			Comp2CumulativeScore.to_frame(name = "DBL_SecondID_Score")], axis = 1)
# 		DBLInfo = pd.concat([DBLInfo,DBLstatus], axis = 0)
# 	return DBLInfo








def QualityMetrics(DropToDBLDict,SparseD,DBLSpecificSingularLociHomoz,DBLSpecificSingularLociHeteroRef,DBLSpecificSingularLociHeteroAlt,vcf, GenotypesDF,SingularLociScoreDF):
	'''

	'''
	DBLInfo = pd.DataFrame( columns = ["DBL_InformativeUMIs","DBL_FirstID_Score","DBL_SecondID_Score"])

#for Doublet in ExtractOrderedDoublets(vcf):

Doublet = ExtractOrderedDoublets(vcf)[0]


#DBLs specific Mask selection
SparseD_BolMask_DBL = SparseD["Barcode"].isin(DropToDBLDict[Doublet])

#Extracting specific loci
HomozlociList=DBLSpecificSingularLociHomoz[Doublet]
RefHetLoci=DBLSpecificSingularLociHeteroRef[Doublet]
AltHetLoci=DBLSpecificSingularLociHeteroAlt[Doublet]


AltReads=SparseD["sparse_Alt"][:,SparseD_BolMask_DBL]
RefReads=SparseD["sparse_Ref"][:,SparseD_BolMask_DBL]



comp1 = Doublet[0]
comp2 = Doublet[1]

HomozRefSingularGenotype = GenotypesDF.loc[HomozlociList,[comp1+"_RefAl",comp2+"_RefAl"]].replace(1,2)
HomozAltSingularGenotype = GenotypesDF.loc[HomozlociList,[comp1+"_AltAl",comp2+"_AltAl"]].replace(1,2)
HeterozRefSingularGenotype = GenotypesDF.loc[RefHetLoci,[comp1+"_RefAl",comp2+"_RefAl"]].replace(1,2)
HeterozAltSingularGenotype = GenotypesDF.loc[AltHetLoci,[comp1+"_AltAl",comp2+"_AltAl"]].replace(1,2)
CleanRefSignal,CleanAltSignal,CleanRef_HET_Signal,CleanAlt_HET_Signal=LocusSpecNoisCalc(DropToDBLDict,Doublet,SingularLociScoreDF,GenotypesDF,vcf,HomozlociList,RefHetLoci,AltHetLoci)




#First thing use HomoLoci information
#First thing use HomoLoci information
#First thing use HomoLoci information


RefReads_HOMOZ=RefReads.loc[HomozlociList].multiply(CleanRefSignal)
AltReads_HOMOZ=AltReads.loc[HomozlociList].multiply(CleanAltSignal)
HomozRefScoresComp1 = RefReads_HOMOZ.multiply(HomozRefSingularGenotype[comp1+"_RefAl"], axis = 0)
HomozRefScoresComp2 = RefReads_HOMOZ.multiply(HomozRefSingularGenotype[comp2+"_RefAl"], axis = 0)
HomozAltScoresComp1 = AltReads_HOMOZ.multiply(HomozAltSingularGenotype[comp1+"_AltAl"], axis = 0)
HomozAltScoresComp2 = AltReads_HOMOZ.multiply(HomozAltSingularGenotype[comp2+"_AltAl"], axis = 0)
HomozComp1Total = HomozRefScoresComp1.add(HomozAltScoresComp1,fill_value = 0.0)
HomozComp2Total = HomozRefScoresComp2.add(HomozAltScoresComp2,fill_value = 0.0)










#Now info about Singular Ref loci
#Now info about Singular Ref loci
#Now info about Singular Ref loci
RefReads_HETEROZ=RefReads.loc[RefHetLoci].multiply(CleanRef_HET_Signal)
HeterozRefScoreComp1 = RefReads_HETEROZ.multiply(HeterozRefSingularGenotype[comp1+"_RefAl"], axis = 0)
HeterozRefScoreComp2 = RefReads_HETEROZ.multiply(HeterozRefSingularGenotype[comp2+"_RefAl"], axis = 0)
#Now those of Singular Alt loci
#Now those of Singular Alt loci
#Now those of Singular Alt loci
AltReads_HETEROZ=AltReads.loc[AltHetLoci].multiply(CleanAlt_HET_Signal)
HeterozAltScoreComp1 = AltReads_HETEROZ.multiply(HeterozAltSingularGenotype[comp1+"_AltAl"], axis = 0)
HeterozAltScoreComp2 = AltReads_HETEROZ.multiply(HeterozAltSingularGenotype[comp2+"_AltAl"], axis = 0)
### TOTAL SCORES
Comp1CumulativeScore = HomozComp1Total.sum()+HeterozRefScoreComp1.sum()+HeterozAltScoreComp1.sum()
Comp2CumulativeScore = HomozComp2Total.sum()+HeterozRefScoreComp2.sum()+HeterozAltScoreComp2.sum()
TotalScore = Comp1CumulativeScore.add(Comp2CumulativeScore,fill_value = 0.0)
#SNG score will be the ratio between reads assigned to 1st id and reads assigned to second id
#Wrapping together into DF slice
DBLstatus = pd.concat([InformativeReads.to_frame(name = "DBL_InformativeUMIs"),
	Comp1CumulativeScore.to_frame(name = "DBL_FirstID_Score"),
	Comp2CumulativeScore.to_frame(name = "DBL_SecondID_Score")], axis = 1)
DBLInfo = pd.concat([DBLInfo,DBLstatus], axis = 0)



	return DBLInfo
