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




def AmbientFactorsEstimate_wEmpty(DBLmetricsDF, DBLsList, FullDrops, HighOutlierThreshold,LowOutlierThreshold):
	DBLmetricsDF_Empty = DBLmetricsDF.loc[~DBLmetricsDF.index.isin(FullDrops)]
	DBLmetricsDF_Empty = DBLmetricsDF_Empty[DBLmetricsDF_Empty.DBL_FirstID_Score != 0]
	IDs = list(DBLmetricsDF_Empty.columns[DBLmetricsDF_Empty.columns.str.startswith("ID")])
	
	# QuantileMask=pd.DataFrame.from_records([[DBLmetricsDF_Empty[genotype].quantile(HighOutlierThreshold) for genotype in IDs]]*len(DBLmetricsDF_Empty),index =  DBLmetricsDF_Empty.index, columns = IDs )-DBLmetricsDF_Empty[IDs] >= 0
	# DBLmetricsDF_Empty = DBLmetricsDF_Empty.loc[QuantileMask[QuantileMask.all(axis=1)].index]
	
	
	
	
	# NormalizedNoise=(DBLmetricsDF_Empty[IDs]).divide(DBLmetricsDF_Empty.loc[DBLmetricsDF_Empty.index,"DBL_FirstID_Score"], axis = 0)
	# AmbientFactors=NormalizedNoise.mean()
	
	
	#Doublet combinations found
	DBLscomb = set(zip(DBLmetricsDF_Empty.FirstID, DBLmetricsDF_Empty.SecondID))
	#
	#Model the noise/ambient only using IDs for each droplet not marked as Best nor Second best ID
	AggregatedNoiseDF = pd.DataFrame()
	for comb in list(DBLscomb):
		NotDBLids = [id for id in IDs if id not in ["ID_"+comb[0],"ID_"+comb[1],"ID"]]
		SliceNoise = DBLmetricsDF_Empty[NotDBLids]
		SliceNoise = DBLmetricsDF_Empty.loc[(DBLmetricsDF_Empty['FirstID'] == comb[0]) & (DBLmetricsDF_Empty['SecondID'] == comb[1]) ,NotDBLids]
		SliceNoise = SliceNoise[NotDBLids]
		AggregatedNoiseDF=pd.concat([SliceNoise, AggregatedNoiseDF])
	
	#Noise normalization vs First ID score
	NormalizedNoise=AggregatedNoiseDF.divide(DBLmetricsDF_Empty.loc[AggregatedNoiseDF.index,"DBL_FirstID_Score"], axis = 0)
	
	#NormalizeFactor per ID by averaging Noises
	AmbientFactors=NormalizedNoise.mean()
	return AmbientFactors



def AddPseudoCounts_wEmpty(DBLmetricsDF, AmbientFactors, DBLsList, FullDrops):
	DBLmetricsDF_Full = DBLmetricsDF.loc[DBLmetricsDF.index.isin(FullDrops)]
	
	NoisedIDs = DBLmetricsDF_Full.loc[DBLmetricsDF_Full.index.difference(DBLsList)]
	#Add 2 columns with correct ambientFactors
	NoisedIDs["NoiseFactor_ID1"] = "ID_"+NoisedIDs["FirstID"]
	NoisedIDs["NoiseFactor_ID2"] = "ID_"+NoisedIDs["SecondID"]
	NoisedIDs["NoiseFactor_ID1"]=NoisedIDs["NoiseFactor_ID1"].replace(  NoisedIDs["NoiseFactor_ID1"].unique().tolist(), AmbientFactors.loc[NoisedIDs["NoiseFactor_ID1"].unique().tolist()].values.tolist()  )
	NoisedIDs["NoiseFactor_ID2"]=NoisedIDs["NoiseFactor_ID2"].replace(  NoisedIDs["NoiseFactor_ID2"].unique().tolist(), AmbientFactors.loc[NoisedIDs["NoiseFactor_ID2"].unique().tolist()].values.tolist()  )
	'''
	Add Ambient factor to Both ID scores, **ambient rate is calculate on total of First ID (more often != 0 and more representative of droplet counts)**
	and coherent with First ID used for normalization in AmbientFactorsEstimate()
	'''
	NoisedIDs["Noised_FirstID_Score"] = NoisedIDs["DBL_FirstID_Score"]+(NoisedIDs["NoiseFactor_ID1"].multiply(NoisedIDs["DBL_FirstID_Score"]))
	NoisedIDs["Noised_SecondID_Score"] = NoisedIDs["DBL_SecondID_Score"]+(NoisedIDs["NoiseFactor_ID2"].multiply(NoisedIDs["DBL_FirstID_Score"]))
	NoisedIDs = NoisedIDs[["Noised_FirstID_Score","Noised_SecondID_Score"]]
	return NoisedIDs


def mixtureModel_wEmpty(NoisedIDs, outdir, DBLsList,DBLmetricsDF, HighOutlierThreshold,LowOutlierThreshold):
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
	
	return QualDF,SafeLowQual,SafeGoodQual,Doublets


def intersectKNNs(QualDF,SafeLowQual,SafeGoodQual,Doublets,FullDropsKNNseries,DBLmetricsDF, FullDrops):
	DBLmetricsDF_Full = DBLmetricsDF.loc[DBLmetricsDF.index.isin(FullDrops)]
	
	
	QualDF_KNN = pd.concat([QualDF,FullDropsKNNseries.loc[QualDF.index]], axis =1)
	QualDF_KNN.loc[(QualDF_KNN["Quality"] == "LowQuality")  &  (QualDF_KNN["FullDropsKNN"] <= 0.5),"Final_Quality"] = "LowQuality"
	
	QualDF_KNN.Quality = QualDF_KNN.Final_Quality.fillna("GoodQuality")
	
	
	Cell_IDs = pd.concat([DBLmetricsDF_Full, (pd.concat([QualDF_KNN["Quality"].to_frame(), SafeLowQual,SafeGoodQual,Doublets]))], axis = 1)
	
	return Cell_IDs




def main_FlagLowQual_wEmpty(DBLmetricsDF, DBLsList, outdir,FullDrops,FullDropsKNNseries, HighOutlierThreshold = .95,LowOutlierThreshold = .05):
	AmbientFactors = AmbientFactorsEstimate_wEmpty(DBLmetricsDF, DBLsList, FullDrops, HighOutlierThreshold,LowOutlierThreshold)
	NoisedIDs = AddPseudoCounts_wEmpty(DBLmetricsDF, AmbientFactors, DBLsList, FullDrops)
	QualDF,SafeLowQual,SafeGoodQual,Doublets = mixtureModel_wEmpty(NoisedIDs, outdir, DBLsList,DBLmetricsDF, HighOutlierThreshold,LowOutlierThreshold)
	Cell_IDs = intersectKNNs(QualDF,SafeLowQual,SafeGoodQual,Doublets,FullDropsKNNseries,DBLmetricsDF, FullDrops)
	
	return Cell_IDs








