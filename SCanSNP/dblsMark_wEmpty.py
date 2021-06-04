
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
import os


fitdistPath=os.path.join(os.path.dirname(__file__), "fitdistrplus/fitdist.R")
quantilePath=os.path.join(os.path.dirname(__file__), "fitdistrplus/quantile.R")


def main__DBLsMark_wEmpty(DBLmetricsDF, FullDrops):
	'''
	Iterating over ranger of Second ID / Total droplet Private reads contributions
	To identify peaks in Doublets detection and best thresholding
	'''
	r=robjects.r
	r.source(fitdistPath)
	r.source(quantilePath)
	print("Fitting negbin on EmptyDrops ")
	
	DBLmetricsDF_Empty = DBLmetricsDF.loc[~DBLmetricsDF.index.isin(FullDrops)]
	DBLmetricsDF_Full = DBLmetricsDF.loc[DBLmetricsDF.index.isin(FullDrops)]
	
	ReadSignals = pd.DataFrame(index = DBLmetricsDF_Full.loc[FullDrops].index, dtype="int")
	
	for id in DBLmetricsDF_Empty.FirstID.unique():
		IDNoise = DBLmetricsDF_Empty["".join(["ID_", id])]
		#We try to trim top 5% of droplet signal before fitting negbin in order to remove outlier
		CroppedNoise = IDNoise[IDNoise < IDNoise.quantile(.95)]
		if sum(CroppedNoise) == 0:
			print("Impossible predict positive threshold for Genotype "+ id + ", it won't be included in doublets detection")
			continue
		distFitt=r.fitdist(robjects.vectors.FloatVector(CroppedNoise), distr = "nbinom")
		Threshold=r.quantile(distFitt, probs = .99)[0][0][0]
		Positiveid = DBLmetricsDF_Full[DBLmetricsDF_Full["".join(["ID_",id])] > Threshold].index
		ReadSignals.loc[Positiveid,id] = 1
		ReadSignals=ReadSignals.fillna(0)
	
	
	DBLs = list(ReadSignals[ReadSignals.sum(axis = 1) > 1].index)
	
	return DBLs


