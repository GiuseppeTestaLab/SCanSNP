#!/usr/bin/env python

#Some function for doublets handle

import itertools
from VCFUtils import *
from sklearn.linear_model import LogisticRegression


def ExtractOrderedDoublets(vcf):
	'''
	Extraction of possible doublets belonging to different donors
	'''
	doublets=list(itertools.permutations(ExtractSamples(vcf), 2))
	return doublets



def FindSecondID(SingularLociScore,BestInDropDict):
	SecondBestIDList = pd.Series([])
	for BestID in list(BestInDropDict.keys()):
		SingularScoreIDspec = SingularLociScore.loc[BestInDropDict[BestID],list(set(ExtractSamples(vcf)) - set([BestID]))]
		SingularScoreIDspecSS = SingularScoreIDspec[SingularScoreIDspec.apply(lambda row: row.nlargest(2).values[-1],axis=1) != SingularScoreIDspec.apply(lambda row: row.nlargest(1).values[-1],axis=1)]
		SecondBestIDList = SecondBestIDList.append(SingularScoreIDspecSS.idxmax(axis = 1))
	DBLsDF = pd.DataFrame( columns = ["FirstID","SecondID"], index = [ cell for cell in barcodeList])
	DBLsDF = pd.concat([DBLsDF, LikeliHoodsDF2.idxmax(axis = 1).to_frame(name = "BestScoringLL")], axis=1)
	DBLsDF = pd.concat([DBLsDF, LikeliHoodsDF2.T.apply(lambda x: x.nlargest(2).idxmin()).to_frame(name = "SecondScoringLL")], axis=1)
	DBLsDF = pd.concat([DBLsDF, SecondBestIDList.to_frame(name = "BestOverallSecondID")], axis=1)
	DBLsDF["FirstID"] = DBLsDF["BestScoringLL"]
	DBLsDF["SecondID"] = DBLsDF["BestOverallSecondID"]
	DBLsDF.SecondID.fillna(DBLsDF.SecondScoringLL, inplace=True)
	DBLsDF = DBLsDF[["FirstID","SecondID"]]
	return DBLsDF




def CalculateSpecs(DBLsDF,BestInDropDict,SingularLociScore):
	SpecsTable=pd.DataFrame( index = barcodeList, columns = ["FirstToSecond_Ratio","FirstToSecond_logRatio","BestIDcounts","SecondIDcounts","NoiseCounts","InformativeUMIs"]).fillna(0)
	#First we register the ratio between the Overall SCORE of First ID and Second ID across all loci
	MainScore_FirstID = LikeliHoodsDF2.lookup(list(DBLsDF.index), list(DBLsDF["FirstID"]))
	MainScore_SecondID = LikeliHoodsDF2.lookup(list(DBLsDF.index), list(DBLsDF["SecondID"]))
	SpecsTable["FirstToSecond_Ratio"] = (MainScore_FirstID/MainScore_SecondID).tolist()
	SpecsTable["FirstToSecond_logRatio"] = (np.log(MainScore_FirstID/MainScore_SecondID)).tolist()
	#We then calculate the informative UMIs as the UMIs mapping to loci with singular information in cohort
	SpecsTable["BestIDcounts"] = SingularLociScore.lookup(list(DBLsDF.index), list(DBLsDF["FirstID"])).tolist()
	SpecsTable["SecondIDcounts"] = SingularLociScore.lookup(list(DBLsDF.index), list(DBLsDF["SecondID"])).tolist()
	SpecsTable["InformativeUMIs"] = SingularLociScore.sum(axis = 1)
	#We then calculate the "Confident Noise" as the reads mapping to unique locy not ascribable to first nor second genotype
	SpecsTable["NoiseCounts"] = SpecsTable["InformativeUMIs"]-(SpecsTable["BestIDcounts"]+SpecsTable["SecondIDcounts"])
	return SpecsTable

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


def DBLScore(DropToDBLDict,countDF,DBLSpecificSingularLociHomoz,DBLSpecificSingularLociHeteroRef,DBLSpecificSingularLociHeteroAlt,vcf, GenotypesDF,SingularLociScoreDF):
	DBLInfo = pd.DataFrame( columns = ["DBL_InformativeUMIs","DBL_FirstID_Score","DBL_SecondID_Score"])
	for Doublet in ExtractOrderedDoublets(vcf):
		HomozlociList=DBLSpecificSingularLociHomoz[Doublet]
		RefHetLoci=DBLSpecificSingularLociHeteroRef[Doublet]
		AltHetLoci=DBLSpecificSingularLociHeteroAlt[Doublet]
		AltReads=countDF[[Drop+'_AltReads' for Drop in DropToDBLDict[Doublet]]]
		RefReads=countDF[[Drop+'_RefReads' for Drop in DropToDBLDict[Doublet]]]
		AltReads.columns = [col.split("_")[0] for col in  AltReads.columns]
		RefReads.columns = [col.split("_")[0] for col in  RefReads.columns]
		InformativeReads=AltReads.loc[AltHetLoci].sum()+RefReads.loc[RefHetLoci].sum()+RefReads.loc[HomozlociList].sum()+AltReads.loc[HomozlociList].sum()
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


def DoubletSpecificSingularLociScan(vcf,GenotypesDF, cleanLoci):
	'''
	Here we extract SingularLoci for every possible genotypes couple (putative doublet)
	'''
	#Generating genotypes from Index and Colnames
	GenotypesClean = GenotypesDF.loc[GenotypesDF.index.isin(cleanLoci)]
	DBLSpecificSingularLociHomoz = {}
	DBLSpecificSingularLociHeteroRef = {}
	DBLSpecificSingularLociHeteroAlt = {}
	SNGLociList = []
	for PutativeDBL in ExtractOrderedDoublets(vcf):
		Comp1 = PutativeDBL[0]
		Comp2 = PutativeDBL[1]
		HomoZSingularLoci = list(GenotypesClean[((GenotypesClean[Comp1+'_RefAl'] == 2) & (GenotypesClean[Comp2+'_AltAl'] == 2 )) | ((GenotypesClean[Comp1+'_AltAl'] == 2) & (GenotypesClean[Comp2+'_RefAl'] == 2 )) ].index)
		HeterozSingularAltLoci = list(GenotypesClean[((GenotypesClean[Comp1+'_RefAl'] == 2) & (GenotypesClean[Comp2+'_AltAl'] == 1 )) | ((GenotypesClean[Comp1+'_AltAl'] == 1) & (GenotypesClean[Comp2+'_RefAl'] == 2 )) ].index)
		HeterozSingularRefLoci = list(GenotypesClean[((GenotypesClean[Comp1+'_RefAl'] == 1) & (GenotypesClean[Comp2+'_AltAl'] == 2 )) | ((GenotypesClean[Comp1+'_AltAl'] == 2) & (GenotypesClean[Comp2+'_RefAl'] == 1 )) ].index)
		DBLSpecificSingularLociHomoz[PutativeDBL] =  list(set(HomoZSingularLoci))
		DBLSpecificSingularLociHeteroAlt[PutativeDBL] =  list(set(HeterozSingularAltLoci))
		DBLSpecificSingularLociHeteroRef[PutativeDBL] =  list(set(HeterozSingularRefLoci))
		SNGLociList.extend(list(set(HeterozSingularRefLoci+HeterozSingularAltLoci+HomoZSingularLoci)))
	SNGLociList = list(set(SNGLociList))
	return SNGLociList,DBLSpecificSingularLociHomoz,DBLSpecificSingularLociHeteroAlt,DBLSpecificSingularLociHeteroRef


def NoiseRregression(BestInDropDict, barcodeList, vcf, SingularLociScoreDF, BestID):
	LogReg = LogisticRegression(multi_class = "multinomial", solver =  "newton-cg", max_iter=100000, penalty = "l2" )
	DBLsDF = pd.DataFrame()
	NotNormScore = SingularLociScoreDF[ExtractSamples(vcf)]
	for FirstID in list(BestInDropDict.keys()):
		BestIDSlice = NotNormScore.loc[BestInDropDict[FirstID]]
		PredictSet = NotNormScore.loc[BestInDropDict[FirstID],list(set(ExtractSamples(vcf)) - set([FirstID]))]
		TrainingSet = NotNormScore.loc[list(set(barcodeList) - set(BestInDropDict[FirstID])),list(set(ExtractSamples(vcf)) - set([FirstID]))]
		TrainingSet["BestID"] = BestID.loc[list(set(barcodeList ) - set(BestInDropDict[FirstID]))]
		TrainingSet = TrainingSet[TrainingSet.sum(axis = 1) > 0]
		FittedModel = LogReg.fit(TrainingSet.iloc[:,:-1],TrainingSet.iloc[:,-1])
		ClassPrediction = FittedModel.predict(PredictSet)
		ClassPredictionFrame = pd.Series(ClassPrediction, index = BestIDSlice.index).to_frame(name = "SecondID")
		BestIDSlice["FirstID"] = FirstID
		BestIDSlice = pd.concat([BestIDSlice, ClassPredictionFrame], axis = 1)
		DBLsDF = pd.concat([DBLsDF, BestIDSlice], axis = 0)
	return DBLsDF[["FirstID","SecondID"]]



def SingularLociScore( SingularLoci_Alt, SingularLoci_Ref, countDf, barcodeList, GenotypesDF,vcf):
	AltSingularGenotype = GenotypesDF.loc[SingularLoci_Alt,[ID+'_AltAl' for ID in list(ExtractSamples(vcf))]].replace(1,2)
	RefSingularGenotype =  GenotypesDF.loc[SingularLoci_Ref,[ID+'_RefAl' for ID in list(ExtractSamples(vcf))]].replace(1,2)
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
	NormScores = SingularLociScores.divide(SingularLociScores.sum(axis = 1), axis = 0)
	NormScores.columns = [Samp + "_Norm" for Samp in NormScores.columns]
	SingularLociScores = pd.concat([SingularLociScores,NormScores], axis = 1)
	return SingularLociScores




def SingularLociScoreNorm( SingularLoci_Alt, SingularLoci_Ref, countDF, barcodeList, GenotypesDF,vcf):
	AltSingularGenotype = GenotypesDF.loc[SingularLoci_Alt,[ID+'_AltAl' for ID in list(ExtractSamples(vcf))]].replace([1, 2],[2, 1])
	RefSingularGenotype =  GenotypesDF.loc[SingularLoci_Ref,[ID+'_RefAl' for ID in list(ExtractSamples(vcf))]].replace([1, 2],[2, 1])
	AltCuntsSS= countDF.loc[SingularLoci_Alt,[barcode+'_AltReads' for barcode in barcodeList]]
	#Total Counts per Alt singular locus will be useful for normalization
	TotalAltCountSS_Alt = countDF.loc[SingularLoci_Alt,[barcode+'_AltReads' for barcode in barcodeList]]
	TotalAltCountSS_Ref = countDF.loc[SingularLoci_Alt,[barcode+'_RefReads' for barcode in barcodeList]]
	TotalAltCountSS_Alt.columns = [col.split("_")[0] for col in  TotalAltCountSS_Alt.columns]
	TotalAltCountSS_Ref.columns = [col.split("_")[0] for col in  TotalAltCountSS_Ref.columns]
	TotalAltCountSS = TotalAltCountSS_Alt.add(TotalAltCountSS_Ref,fill_value = 0.0)
	RefCuntsSS= countDF.loc[SingularLoci_Ref,[barcode+'_RefReads' for barcode in barcodeList]]
	#Total Counts per Ref singular locus will be useful for normalization
	TotalRefCountSS_Alt = countDF.loc[SingularLoci_Ref,[barcode+'_AltReads' for barcode in barcodeList]]
	TotalRefCountSS_Ref = countDF.loc[SingularLoci_Ref,[barcode+'_RefReads' for barcode in barcodeList]]
	TotalRefCountSS_Alt.columns = [col.split("_")[0] for col in  TotalRefCountSS_Alt.columns]
	TotalRefCountSS_Ref.columns = [col.split("_")[0] for col in  TotalRefCountSS_Ref.columns]
	TotalRefCountSS = TotalRefCountSS_Alt.add(TotalRefCountSS_Ref,fill_value = 0.0)
	ScorePerBArcodeAlt = pd.DataFrame( columns = list(ExtractSamples(vcf)), index = barcodeList).fillna(0.0)
	ScorePerBArcodeRef = pd.DataFrame( columns = list(ExtractSamples(vcf)), index = barcodeList).fillna(0.0)
	IDsToQueryAlt = [ID+'_AltAl' for ID in list(ExtractSamples(vcf))]
	IDsToQueryRef = [ID+'_RefAl' for ID in list(ExtractSamples(vcf))]
	#PileUp for Alt informative loci
	for ID in IDsToQueryAlt:
		AltGenotypeID = AltSingularGenotype[ID]
		AltCuntsSS.columns = barcodeList
		ScorePerBArcodeAlt["_".join(ID.split("_")[:-1])] = AltCuntsSS.divide(TotalAltCountSS).multiply(AltGenotypeID, axis = 0).sum()
	#PileUp for Ref informative loci
	for ID in IDsToQueryRef:
		RefGenotypeID = RefSingularGenotype[ID]
		RefCuntsSS.columns = barcodeList
		ScorePerBArcodeRef["_".join(ID.split("_")[:-1])] = RefCuntsSS.divide(TotalRefCountSS).multiply(RefGenotypeID, axis = 0).sum()
	SingularLociScoresII = ScorePerBArcodeRef.add(ScorePerBArcodeAlt,fill_value = 0.0)
	return SingularLociScoresII
