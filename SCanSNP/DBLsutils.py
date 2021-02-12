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
