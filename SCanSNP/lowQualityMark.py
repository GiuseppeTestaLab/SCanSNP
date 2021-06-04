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
		NotDBLids = [id for id in IDs if id not in ["ID_"+comb[0],"ID_"+comb[1],"ID"]]
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



# def mixtureModel(NoisedIDs, outdir, DBLsList,DBLmetricsDF, HighOutlierThreshold = .95,LowOutlierThreshold = .05 ):
# 	logFC_DF = NoisedIDs.copy()
	
# 	#Calc logFC
# 	logFC_DF["logFC"] = np.log(logFC_DF["Noised_FirstID_Score"].divide(logFC_DF["Noised_SecondID_Score"]))
	
	
# 	#Retrieve Best and Second Best IDs for next fitting
# 	logFC_DF = pd.concat([logFC_DF,DBLmetricsDF.loc[list(logFC_DF.index), ["FirstID","SecondID"]]], axis = 1)
	
# 	#Doublet combinations 
# 	DBLscomb = list(itertools.combinations(list(set(list(SNGsmetricsDF.FirstID) + list(SNGsmetricsDF.SecondID))), 2))
	
# 	# Fit mixModel on non-outliers in Pairwise manner
# 	QualDF = pd.DataFrame(columns = ["Quality"])
# 	for comb in list(DBLscomb):
# 		logFC_DFDBL =  logFC_DF[(logFC_DF["FirstID"].isin(comb)) & (logFC_DF["SecondID"].isin(comb))]
# 		#
# 		logFC_DFDBL.shape
# 		#Marking low and high outliers as bonafide LowQual and GoodQual
# 		SafeLowQual = logFC_DFDBL[logFC_DFDBL["logFC"] < logFC_DFDBL["logFC"].quantile(LowOutlierThreshold)].index.tolist()
# 		SafeGoodQual = logFC_DFDBL[logFC_DFDBL["logFC"] > logFC_DFDBL["logFC"].quantile(HighOutlierThreshold)].index.tolist()
# 		len(SafeLowQual)
# 		len(SafeGoodQual)
# 		#
		
# 		logFC_DFDBL = logFC_DFDBL[(logFC_DFDBL["logFC"] <= logFC_DFDBL["logFC"].quantile(HighOutlierThreshold)) & (logFC_DFDBL["logFC"] >= logFC_DFDBL["logFC"].quantile(LowOutlierThreshold))]
# 		logFC_DFDBL.shape
# 		gm = GaussianMixture(n_components=2, random_state=0, weights_init = [.05,.95]).fit(logFC_DFDBL["logFC"].to_numpy().reshape(-1,1))
# 		QualLabels = gm.predict(logFC_DFDBL["logFC"].to_numpy().reshape(-1,1))
# 		#
# 		#Plot Fitting
# 		FittedMixturePlot(gm, outdir,comb, logFC_DFDBL)
# 		#
# 		#Flagging lowquals
# 		QualFlag = np.where(QualLabels==gm.means_.argmax(), "GoodQuality", "LowQuality")
# 		QualDF_Pair=pd.DataFrame(QualFlag, index = logFC_DFDBL.index, columns = ["Quality"])
# 		#
# 		SafeLowQual = pd.DataFrame(["LowQuality" for x in range(len(SafeLowQual))], index = SafeLowQual, columns = ["Quality"])
# 		SafeGoodQual = pd.DataFrame(["GoodQuality" for x in range(len(SafeGoodQual))], index = SafeGoodQual, columns = ["Quality"])
		
# 		#
# 		QualDF_Pair = pd.concat([QualDF_Pair, SafeLowQual, SafeGoodQual])
# 		QualDF_Pair.shape
# 		#
# 		QualDF = pd.concat([QualDF, QualDF_Pair], axis = 0)
# 	#
# 	Doublets = pd.DataFrame(["Doublet" for x in range(len(DBLsList))], index = DBLsList, columns = ["Quality"])
# 	QualDF = pd.concat([QualDF, Doublets])
# 	return(QualDF)



def mixtureModel(NoisedIDs, outdir, DBLsList,DBLmetricsDF, HighOutlierThreshold = .95,LowOutlierThreshold = .05 ):
	NoisedIDs_processed=NoisedIDs.copy()
	
	NoisedIDs_processed["logFC"] = np.log(NoisedIDs_processed["Noised_FirstID_Score"].divide(NoisedIDs_processed["Noised_SecondID_Score"]))
	#Marking low and high outliers as bonafitr LowQual and GoodQual
	SafeLowQual = NoisedIDs_processed[NoisedIDs_processed["logFC"] <= NoisedIDs_processed["logFC"].quantile(LowOutlierThreshold)].index.tolist()
	SafeGoodQual = NoisedIDs_processed[NoisedIDs_processed["logFC"] >= NoisedIDs_processed["logFC"].quantile(HighOutlierThreshold)].index.tolist()
	# Fit mixModel on non-outliers
	NoisedIDs_processed = NoisedIDs_processed[(NoisedIDs_processed["logFC"] < NoisedIDs_processed["logFC"].quantile(HighOutlierThreshold)) & (NoisedIDs_processed["logFC"] > NoisedIDs_processed["logFC"].quantile(LowOutlierThreshold))]
	gm = GaussianMixture(n_components=2, random_state=0, weights_init = [.05,.95]).fit(NoisedIDs_processed["logFC"].to_numpy().reshape(-1,1))
	QualLabels = gm.predict(NoisedIDs_processed["logFC"].to_numpy().reshape(-1,1))
	
	#Plot Fitting
	FittedMixturePlot(gm, outdir, NoisedIDs_processed)
	
	#Flagging lowquals
	QualFlag = np.where(QualLabels==gm.means_.argmax(), "GoodQuality", "LowQuality")
	QualDF=pd.DataFrame(QualFlag, index = NoisedIDs_processed.index, columns = ["Quality"])
	
	SafeLowQual = pd.DataFrame(["LowQuality" for x in range(len(SafeLowQual))], index = SafeLowQual, columns = ["Quality"])
	SafeGoodQual = pd.DataFrame(["GoodQuality" for x in range(len(SafeGoodQual))], index = SafeGoodQual, columns = ["Quality"])
	Doublets = pd.DataFrame(["Doublet" for x in range(len(DBLsList))], index = DBLsList, columns = ["Quality"])
	
	QualDF = pd.concat([QualDF, SafeLowQual, SafeGoodQual, Doublets])
	
	return QualDF



def main_FlagLowQual(DBLmetricsDF, DBLsList, outdir):
	AmbientFactors = AmbientFactorsEstimate(DBLmetricsDF, DBLsList)
	NoisedIDs = AddPseudoCounts(DBLmetricsDF, AmbientFactors, DBLsList)
	QualDF = mixtureModel(NoisedIDs, outdir, DBLsList,DBLmetricsDF, HighOutlierThreshold = .95,LowOutlierThreshold = .05)
	
	Cell_IDs = pd.concat([DBLmetricsDF, QualDF], axis = 1)
	
	
	return Cell_IDs
