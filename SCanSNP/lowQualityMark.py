#!/usr/bin/env python

import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.mixture import GaussianMixture
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm
import numpy as np
from scipy import stats
from plots import FittedMixturePlot




def AmbientFactorsEstimate(DBLmetricsDF, DBLsList):
	SNGsmetricsDF = DBLmetricsDF.loc[DBLmetricsDF.index.difference(DBLsList)]
	SNGsmetricsDF = SNGsmetricsDF[SNGsmetricsDF.DBL_FirstID_Score != 0]
	IDs = list(SNGsmetricsDF.columns[SNGsmetricsDF.columns.str.startswith("ID")])
	#
	#
	#Doublet combinations found
	DBLscomb = set(zip(SNGsmetricsDF.FirstID, SNGsmetricsDF.SecondID))
	#
	#Model the noise/ambient only using IDs for each droplet not marked as Best nor Second best ID
	AggregatedNoiseDF = pd.DataFrame()
	for comb in list(DBLscomb):
		NotDBLids = [id for id in IDs if id not in ["ID_"+comb[0],"ID_"+comb[1]]]
		SliceNoise = SNGsmetricsDF[NotDBLids]
		SliceNoise = SNGsmetricsDF.loc[(SNGsmetricsDF['FirstID'] == comb[0]) & (SNGsmetricsDF['SecondID'] == comb[1]) ,NotDBLids]
		SliceNoise = SliceNoise[NotDBLids]
		AggregatedNoiseDF=pd.concat([SliceNoise, AggregatedNoiseDF])
	#
	#Noise normalization vs First ID score
	NormalizedNoise=AggregatedNoiseDF.divide(SNGsmetricsDF.loc[AggregatedNoiseDF.index,"DBL_FirstID_Score"], axis = 0)
	#
	#NormalizeFactor per ID by averaging Noises
	AmbientFactors=NormalizedNoise.mean()
	return AmbientFactors



def AddPseudoCounts(DBLmetricsDF, AmbientFactors, DBLsList):
	NoisedIDs = DBLmetricsDF.loc[DBLmetricsDF.index.difference(DBLsList)]
	#Add 2 columns with correct ambientFactors
	NoisedIDs["NoiseFactor_ID1"] = "ID_"+NoisedIDs["FirstID"]
	NoisedIDs["NoiseFactor_ID2"] = "ID_"+NoisedIDs["SecondID"]
	NoisedIDs["NoiseFactor_ID1"]=NoisedIDs["NoiseFactor_ID1"].replace(  NoisedIDs["NoiseFactor_ID1"].unique().tolist(),    AmbientFactors.loc[NoisedIDs["NoiseFactor_ID1"].unique().tolist()].values.tolist()  )
	NoisedIDs["NoiseFactor_ID2"]=NoisedIDs["NoiseFactor_ID2"].replace(  NoisedIDs["NoiseFactor_ID2"].unique().tolist(),    AmbientFactors.loc[NoisedIDs["NoiseFactor_ID2"].unique().tolist()].values.tolist()  )
	'''
	Add Ambient factor to Both ID scores, **ambient rate is calculate on total of First ID (more often != 0 and more representative of droplet counts)**
	and coherent with First ID used for normalization in AmbientFactorsEstimate()
	'''
	NoisedIDs["Noised_FirstID_Score"] = NoisedIDs["DBL_FirstID_Score"]+(NoisedIDs["NoiseFactor_ID1"].multiply(NoisedIDs["DBL_FirstID_Score"]))
	NoisedIDs["Noised_SecondID_Score"] = NoisedIDs["DBL_SecondID_Score"]+(NoisedIDs["NoiseFactor_ID2"].multiply(NoisedIDs["DBL_FirstID_Score"]))
	NoisedIDs = NoisedIDs[["Noised_FirstID_Score","Noised_SecondID_Score"]]
	return NoisedIDs



def mixtureModel(NoisedIDs, outdir, DBLsList, HighOutlierThreshold = .95,LowOutlierThreshold = .05 ):
	NoisedIDs["logFC"] = np.log(NoisedIDs["Noised_FirstID_Score"].divide(NoisedIDs["Noised_SecondID_Score"]))
	#Marking low and high outliers as bonafitr LowQual and GoodQual
	SafeLowQual = NoisedIDs[NoisedIDs["logFC"] < NoisedIDs["logFC"].quantile(LowOutlierThreshold)].index.tolist()
	SafeGoodQual = NoisedIDs[NoisedIDs["logFC"] > NoisedIDs["logFC"].quantile(HighOutlierThreshold)].index.tolist()
	# Fit mixModel on non-outliers
	NoisedIDs = NoisedIDs[(NoisedIDs["logFC"] < NoisedIDs["logFC"].quantile(HighOutlierThreshold)) & (NoisedIDs["logFC"] > NoisedIDs["logFC"].quantile(LowOutlierThreshold))]
	gm = GaussianMixture(n_components=2, random_state=0).fit(NoisedIDs["logFC"].to_numpy().reshape(-1,1))
	QualLabels = gm.predict(NoisedIDs["logFC"].to_numpy().reshape(-1,1))
	
	#Plot Fitting
	FittedMixturePlot(gm, outdir, NoisedIDs)
	
	#Flagging lowquals
	QualFlag = np.where(QualLabels==gm.means_.argmax(), "GoodQuality", "LowQuality")
	QualDF=pd.DataFrame(QualFlag, index = NoisedIDs.index, columns = ["Quality"])
	
	SafeLowQual = pd.DataFrame(["LowQuality" for x in range(len(SafeLowQual))], index = SafeLowQual, columns = ["Quality"])
	SafeGoodQual = pd.DataFrame(["GoodQuality" for x in range(len(SafeGoodQual))], index = SafeGoodQual, columns = ["Quality"])
	Doublets = pd.DataFrame(["Doublet" for x in range(len(DBLsList))], index = DBLsList, columns = ["Quality"])
	
	QualDF = pd.concat([QualDF, SafeLowQual, SafeGoodQual, Doublets])
	
	return(QualDF)



def main_FlagLowQual(DBLmetricsDF, DBLsList, outdir):
	AmbientFactors = AmbientFactorsEstimate(DBLmetricsDF, DBLsList)
	NoisedIDs = AddPseudoCounts(DBLmetricsDF, AmbientFactors, DBLsList)
	QualDF = mixtureModel(NoisedIDs, outdir, DBLsList, HighOutlierThreshold = .95,LowOutlierThreshold = .05 )
	
	return QualDF
