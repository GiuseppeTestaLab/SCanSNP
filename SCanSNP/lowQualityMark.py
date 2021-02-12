#!/usr/bin/env python

import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.mixture import GaussianMixture
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm


DBLmetricsDF=pd.read_csv("/hpcnfs/scratch/temporary/Dav_vc/0.3_SATIJA_CellHashing_data/SynthDemuxing2/SynthDemux/CELLHASHING/SCanSNP/DBLmetricsDF.tsv", sep ="\t", index_col=0)



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

AmbientFactors=AmbientFactorsEstimate(DBLmetricsDF, DBLsList)


def AddPseudoCounts(DBLmetricsDF, AmbientFactors):
	SNGsnoisedDF = DBLmetricsDF.loc[DBLmetricsDF.index.difference(DBLsList)]
	#Add 2 columns with correct ambientFactors
	SNGsnoisedDF["NoiseFactor_ID1"] = "ID_"+SNGsnoisedDF["FirstID"]
	SNGsnoisedDF["NoiseFactor_ID2"] = "ID_"+SNGsnoisedDF["SecondID"]
	SNGsnoisedDF["NoiseFactor_ID1"]=SNGsnoisedDF["NoiseFactor_ID1"].replace(  SNGsnoisedDF["NoiseFactor_ID1"].unique().tolist(),    AmbientFactors.loc[SNGsnoisedDF["NoiseFactor_ID1"].unique().tolist()].values.tolist()  )
	SNGsnoisedDF["NoiseFactor_ID2"]=SNGsnoisedDF["NoiseFactor_ID2"].replace(  SNGsnoisedDF["NoiseFactor_ID2"].unique().tolist(),    AmbientFactors.loc[SNGsnoisedDF["NoiseFactor_ID2"].unique().tolist()].values.tolist()  )
	'''
	Add Ambient factor to Both ID scores, **ambient rate is calculate on total of First ID (more often != 0 and more representative of droplet counts)**
	and coherent with First ID used for normalization in AmbientFactorsEstimate()
	'''
	SNGsnoisedDF["Noised_FirstID_Score"] = SNGsnoisedDF["DBL_FirstID_Score"]+(SNGsnoisedDF["NoiseFactor_ID1"].multiply(SNGsnoisedDF["DBL_FirstID_Score"]))
	SNGsnoisedDF["Noised_SecondID_Score"] = SNGsnoisedDF["DBL_SecondID_Score"]+(SNGsnoisedDF["NoiseFactor_ID2"].multiply(SNGsnoisedDF["DBL_FirstID_Score"]))
	SNGsnoisedDF = SNGsnoisedDF[["Noised_FirstID_Score","Noised_SecondID_Score"]]
	return SNGsnoisedDF


SNGsnoisedDF=AddPseudoCounts(DBLmetricsDF, AmbientFactors)


def mixtureModel(SNGsnoisedDF, HighOutlierThresh = .95,LowOutlierThresh = .05 ):

HighOutlierThresh = .95
LowOutlierThresh = .05
SNGsnoisedDF["logFC"] = np.log(SNGsnoisedDF["Noised_FirstID_Score"].divide(SNGsnoisedDF["Noised_SecondID_Score"]))
#Marking low and high outliers as bonafitr LowQual and GoodQual
SafeLowQual = SNGsnoisedDF[SNGsnoisedDF["logFC"] < SNGsnoisedDF["logFC"].quantile(LowOutlierThresh)].index.tolist()
SafeGoodQual = SNGsnoisedDF[SNGsnoisedDF["logFC"] > SNGsnoisedDF["logFC"].quantile(HighOutlierThresh)].index.tolist()
# Fit mixModel on non-outliers
SNGsnoisedDF = SNGsnoisedDF[(SNGsnoisedDF["logFC"] < SNGsnoisedDF["logFC"].quantile(HighOutlierThresh)) & (SNGsnoisedDF["logFC"] > SNGsnoisedDF["logFC"].quantile(LowOutlierThresh))]
gm = GaussianMixture(n_components=2, random_state=0).fit(SNGsnoisedDF["logFC"].to_numpy().reshape(-1,1))
QualLabels = gm.predict(SNGsnoisedDF["logFC"].to_numpy().reshape(-1,1))

#Simulate fitted components
C1 = np.random.normal(gm.means_[0][0], np.sqrt(gm.covariances_[0][0][0]), 10000)
C2 = np.random.normal(gm.means_[1][0], np.sqrt(gm.covariances_[1][0][0]), 10000)


###################ààà

plt.clf()


import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


plt.clf()
figure(num=None, figsize=(20, 10), dpi=80, facecolor='w', edgecolor='k')

C1 = stats.norm(gm.means_[0][0], np.sqrt(gm.covariances_[0][0][0]))
C2 = stats.norm(gm.means_[1][0], np.sqrt(gm.covariances_[1][0][0]))


mc = gm.weights_

x = np.linspace(min(SNGsnoisedDF["logFC"]), max(SNGsnoisedDF["logFC"]), 501)


C1 = C1.pdf(x) * mc[0]
C2 = C2.pdf(x) * mc[1]


gm.means_[0][0] > gm.means_[1][0]

sns.distplot(SNGsnoisedDF["logFC"].to_numpy(), hist=False, label='Mixture')
plt.plot(x, C1,'--', label='Fitted LowQuality' if gm.means_[0][0] < gm.means_[1][0] else 'Fitted_HighQuality' )
plt.plot(x, C2,'--',  label='Fitted LowQuality' if gm.means_[0][0] > gm.means_[1][0] else 'Fitted_HighQuality' )

plt.legend(prop={'size': 16}, title = 'Variants', fontsize=25)
plt.title('logFC Density Plot '+ Var, fontsize=25)
plt.xlabel("ID1-ID2_logFC", fontsize=25)
plt.ylabel('Density', fontsize=25)

plt.savefig("/hpcnfs/scratch/temporary/Dav_vc/output.png")
