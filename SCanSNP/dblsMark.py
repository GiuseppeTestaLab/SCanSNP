
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
import os


fitdistPath=os.path.join(os.path.dirname(__file__), "fitdistrplus/fitdist.R")
quantilePath=os.path.join(os.path.dirname(__file__), "fitdistrplus/quantile.R")


# def iterateFitting(DBLmetricsDF):
# 	'''
# 	Iterating over ranger of Second ID / Total droplet Private reads contributions
# 	To identify peaks in Doublets detection and best thresholding
# 	'''
# 	r=robjects.r
# 	r.source(fitdistPath)
# 	r.source(quantilePath)
# 	print("Iterative fitting negbin for different ID1/ID2 ratios")
	
# 	FittedIterations = pd.DataFrame(index = np.arange(0.1,0.61, 0.01).round(2), dtype="float32")
# 	for SecondIDfraction in np.arange(0.1,0.61, 0.01).round(2):
# 		Outliers="Trimmed"
# 		ReadSignals = pd.DataFrame(index = DBLmetricsDF.index, dtype="int")
# 		ID_specificThresholds = {}
# 		#Defining "negative" droplets for each ID
# 		for id in DBLmetricsDF.FirstID.unique():
# 			IDNoise = DBLmetricsDF.loc[(DBLmetricsDF["FirstID"] != id) & ~((DBLmetricsDF["SecondID"] == id) & (DBLmetricsDF["".join(["ID_",id])] >= DBLmetricsDF[["".join(["ID_", ident]) for ident in DBLmetricsDF.FirstID.unique() ]].sum(axis = 1)*SecondIDfraction)),"".join(["ID_", id])]
# 			#We try to trim top 5% of droplet signal before fitting negbin in order to remove outlier
# 			CroppedNoise = IDNoise[IDNoise < IDNoise.quantile(.95)]
# 			if CroppedNoise.sum() == 0:
# 				#In case trimming fails (result in all 0 values) we fall back on untrimmed values
# 				#storing the result in $Outliers variable for subsequent assessment
# 				print(".5 outliers trimming was not possible")
# 				Outliers="NotTrimmed"
# 				CroppedNoise = IDNoise
# 			elif len(CroppedNoise) < 2:
# 				continue
# 			#Wrapping fitdistrplus::fitdist() and stats::quantile() functions from R scripts
# 			distFitt=r.fitdist(robjects.vectors.FloatVector(CroppedNoise), distr = "nbinom")
# 			Threshold=r.quantile(distFitt, probs = .99)[0][0][0]
# 			ID_specificThresholds["".join([id,"_R.Threshold"])] = Threshold
# 			Positiveid = DBLmetricsDF[DBLmetricsDF["".join(["ID_",id])] > Threshold].index
# 			ReadSignals.loc[Positiveid,id] = 1
# 			ReadSignals=ReadSignals.fillna(0)
# 		#Subsetting for double (or more) positive droplets
# 		MultiSIg = ReadSignals[ReadSignals.sum(axis = 1) > 1]
# 		FittedIterations.loc[SecondIDfraction,"DBLs"] = MultiSIg.shape[0]
# 		FittedIterations.loc[SecondIDfraction,"Outliers"] = Outliers
# 		print("".join(["Completed ", str(round(SecondIDfraction,2)), " ID2 fraction iteration"]))
	
# 	return FittedIterations



# def ID2_RatioSelect(FittedIterations):
# 	'''
# 	For each iteration we store the number of lost doublets with respect of previous one to spot the highest "Loss"
# 	'''
# 	FittedIterationsLoss = FittedIterations.copy()
# 	FittedIterationsLoss["DBLsLoss"] = FittedIterationsLoss["DBLs"].diff().fillna(0.0).abs()
# 	#If 5% outliers trimming failed we ignore the first peak because it will likely be due to non trimmed outliers retention
# 	if FittedIterationsLoss["Outliers"].unique()[0] == "NotTrimmed":
# 		FittedIterationsLoss =  FittedIterationsLoss.iloc[2:]
# 		MaxLoss = FittedIterationsLoss["DBLsLoss"].idxmax()
# 		ID2_Ratio = FittedIterationsLoss.index[FittedIterationsLoss.index.get_loc(MaxLoss)-1]
# 	else:
# 		MaxLoss = FittedIterationsLoss["DBLsLoss"].idxmax()
# 		ID2_Ratio = FittedIterationsLoss.index[FittedIterationsLoss.index.get_loc(MaxLoss)-1]
# 	return ID2_Ratio


# def DBLsMark(DBLmetricsDF, ID2_Ratio):
# 	'''
# 	Perform actual DBLs marking based on selected ID2 ratio negbin fitting
# 	'''
# 	r=robjects.r
# 	r.source(fitdistPath)
# 	r.source(quantilePath)
# 	ReadSignals = pd.DataFrame(index = DBLmetricsDF.index, dtype="int")
# 	for id in DBLmetricsDF.FirstID.unique():
# 		IDNoise = DBLmetricsDF.loc[(DBLmetricsDF["FirstID"] != id) & ~((DBLmetricsDF["SecondID"] == id) & (DBLmetricsDF["".join(["ID_",id])] >= DBLmetricsDF[["".join(["ID_", ident]) for ident in DBLmetricsDF.FirstID.unique() ]].sum(axis = 1)*ID2_Ratio)),"".join(["ID_", id])]
# 		#We try to trim top 5% of droplet signal before fitting negbin in order to remove outlier
# 		CroppedNoise = IDNoise[IDNoise < IDNoise.quantile(.95)]
# 		if CroppedNoise.sum() == 0:
# 			#In case trimming fails (result in all 0 values) we fall back on untrimmed values
# 			#storing the result in $Outliers variable for subsequent assessment
# 			print(".5 outliers trimming was not possible")
# 			Outliers="NotTrimmed"
# 			CroppedNoise = IDNoise
# 		elif len(CroppedNoise) < 2:
# 			continue
# 		#Wrapping fitdistrplus::fitdist() and stats::quantile() functions from R scripts
# 		distFitt=r.fitdist(robjects.vectors.FloatVector(CroppedNoise), distr = "nbinom")
# 		Threshold=r.quantile(distFitt, probs = .99)[0][0][0]
# 		Positiveid = DBLmetricsDF[DBLmetricsDF["".join(["ID_",id])] > Threshold].index
# 		ReadSignals.loc[Positiveid,id] = 1
# 		ReadSignals=ReadSignals.fillna(0)
# 	#Subsetting for double (or more) positive droplets
# 	DBLs = list(ReadSignals[ReadSignals.sum(axis = 1) > 1].index)
# 	return DBLs


# def main__DBLsMark(DBLmetricsDF):
# 	'''
# 	Wrapping together previous 3 in the main dbls mark function
# 	'''
	
# 	#DBLs detection Module
# 	FittedIterations=iterateFitting(DBLmetricsDF)
# 	ID2_Ratio=ID2_RatioSelect(FittedIterations)
# 	DBLsList=DBLsMark(DBLmetricsDF, ID2_Ratio)
	
# 	return DBLsList


def main__DBLsMark(DBLmetricsDF):
	
	r=robjects.r
	r.source(fitdistPath)
	r.source(quantilePath)
	print("Iterative fitting negbin for different ID1/ID2 ratios")
	
	ReadSignals = pd.DataFrame(index = DBLmetricsDF.index, dtype="int")
	#Defining "negative" droplets for each ID
	SafeSNGs = round(len(DBLmetricsDF)-(len(DBLmetricsDF)/130000*len(DBLmetricsDF)))
	SafeSNGs = DBLmetricsDF.loc[((DBLmetricsDF["DBL_FirstID_Score"]+0.01)/(DBLmetricsDF["DBL_SecondID_Score"]+0.01)).sort_values().tail(SafeSNGs).index]
	for id in DBLmetricsDF.FirstID.unique():
		#FCthresholded = (DBLmetricsDF["DBL_FirstID_Score"]+0.01)/(DBLmetricsDF["DBL_SecondID_Score"]+0.01) > 10
		IDNoise = SafeSNGs.loc[SafeSNGs["FirstID"] != id ,"".join(["ID_", id])]
		#We try to trim top 5% of droplet signal before fitting negbin in order to remove outlier
		CroppedNoise = IDNoise
		if CroppedNoise.sum() == 0:
			print("Impossible predict positive threshold for Genotype "+ id + ", it won't be included in doublets detection")
			continue
		elif len(CroppedNoise) < 2:
			print("Impossible predict positive threshold for Genotype "+ id + ", it won't be included in doublets detection")
			continue
		else:
			#Wrapping fitdistrplus::fitdist() and stats::quantile() functions from R scripts
			try:
				distFitt=r.fitdist(robjects.vectors.FloatVector(CroppedNoise), distr = "nbinom")
				#Using 3 as fallBack threhsold in case of 0 thresholding
				Threshold=r.quantile(distFitt, probs = .99)[0][0][0] if r.quantile(distFitt, probs = .99)[0][0][0] > 0 else 2
				print("For ID {} Threshold is {}".format(id, Threshold))
			except Exception as e:
				print("mle failed to fit nbinom for genotype Genotype {} , with following exception,  it won't be included in doublets detection {}".format(id, e))
				continue
		Positiveid = DBLmetricsDF[DBLmetricsDF["".join(["ID_",id])] > Threshold].index
		ReadSignals.loc[Positiveid,id] = 1
		ReadSignals=ReadSignals.fillna(0)
	DBLs = list(ReadSignals[ReadSignals.sum(axis = 1) > 1].index)
	
	return DBLs
