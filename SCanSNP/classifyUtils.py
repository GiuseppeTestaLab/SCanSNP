import pandas as pd
import numpy as np
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.misc import derivative
import matplotlib.pyplot as plt






def normCenter(DBLmetricsDF,mixedIDs):
	
	mixedIDs = list(set(DBLmetricsDF["FirstID"]).union(set(DBLmetricsDF["SecondID"])))
	
	
	logVal = np.log(DBLmetricsDF[["ID_"+ID for ID in mixedIDs]])
	logVal[logVal < 0 ] = 0
	normCent = logVal-np.mean(logVal.values)
	
	normCent = pd.concat([normCent, DBLmetricsDF[["FirstID","SecondID"]]], axis = 1)
	
	
	return normCent



def singletscontribution(segmentationDF, normCounts, ID):
	
	normCounts_seg = pd.concat([normCounts, segmentationDF], axis = 1)  
	
	#ID_specificNegatives = normCounts_seg[(normCounts_seg["segmentation"] == 1) & (normCounts_seg["FirstID"] != ID)].index
	#ID_specificPositives = normCounts_seg[(normCounts_seg["segmentation"] == 1) & (normCounts_seg["FirstID"] == ID)].index
	
	ID_specificNegatives = normCounts_seg[(    (normCounts_seg["segmentation"] == 1) & (normCounts_seg["FirstID"] != ID)    ) |
	(    (normCounts_seg["segmentation"] == 2) & ((normCounts_seg["FirstID"] != ID) &  (normCounts_seg["SecondID"] != ID))    )
	].index
	
	ID_specificPositives = normCounts_seg[(    (normCounts_seg["segmentation"] == 1) & (normCounts_seg["FirstID"] == ID)    ) |
	(    (normCounts_seg["segmentation"] == 2) & ((normCounts_seg["FirstID"] == ID) | (normCounts_seg["SecondID"] == ID))    )
	].index 
	
	ID_specificNegatives_contr = pd.DataFrame({"normCounts":normCounts_seg.loc[ID_specificNegatives,"ID_"+ID].values,"label":"negative"}, index=ID_specificNegatives)
	ID_specificPositives_contr = pd.DataFrame({"normCounts":normCounts_seg.loc[ID_specificPositives,"ID_"+ID].values,"label":"positive"}, index=ID_specificPositives)
	
	countsDF = pd.concat([ID_specificNegatives_contr,ID_specificPositives_contr], axis = 0)
	
	return countsDF





def definePositiveRates(countsDF):
	
	cuts = np.linspace(np.quantile(countsDF["normCounts"], 0.01), 
	np.quantile(countsDF["normCounts"], 0.99),
	num=10)
	
	bin_posRatesDF = pd.DataFrame({"x":cuts[:-1],"PositiveRate":np.nan})
	
	
	countsDF["bin"] = pd.cut(countsDF["normCounts"], cuts, include_lowest=True, retbins=False, labels=cuts[:-1].copy())
	
	for binval in countsDF["bin"].cat.categories:
		
		positive_bin_Cells = countsDF[(countsDF["bin"] == binval) & (countsDF["label"] == "positive")].shape[0]
		negative_bin_Cells = countsDF[(countsDF["bin"] == binval) & (countsDF["label"] == "negative")].shape[0]
		
		total_bin_Cells = positive_bin_Cells + negative_bin_Cells
		
		PosRate_bin = positive_bin_Cells/total_bin_Cells
		
		bin_posRatesDF.loc[ bin_posRatesDF["x"] == binval, "PositiveRate"] = PosRate_bin
	
	return bin_posRatesDF



def logFun(x, L ,x0, k, b):
	y = L / (1 + np.exp(-k*(x-x0)))+b
	return y



def findThreshold(posRatesDF, ID):
	
	fig, axis = plt.subplots(3,1,figsize=(20,20)) #
	
	#--------------------plot binned ratios
	sns.scatterplot(posRatesDF["x"], posRatesDF["PositiveRate"], hue=posRatesDF["PositiveRate"], ax=axis[0])
	axis[0].set_title("Raw binned ratios "+ID)
	#---------------------------------------
	
	X = posRatesDF["x"].values
	y = posRatesDF["PositiveRate"].values
	
	
	
	bounds = [max(y), np.median(X),1,min(y)] 
	
	fittedParams, pcov = curve_fit(logFun, X, y,bounds, method='dogbox')
	
	#Defining "continous" x values
	continous_cuntSpace = np.linspace(min(posRatesDF["x"]), max(posRatesDF["x"]), 100)
	derivativeDF = pd.DataFrame({"acceleration":np.nan,"countPoint":continous_cuntSpace}, index = continous_cuntSpace )
	
	#--------------------plot fitting
	sns.lineplot(continous_cuntSpace, logFun(continous_cuntSpace, *fittedParams), ax=axis[1])
	axis[1].set_title("Fitted logistic function "+ID)
	#---------------------------------
	
	
	for countPoint in derivativeDF.index.tolist():
		RateDerivate = derivative(lambda x:  fittedParams[0] / (1 + np.exp(-fittedParams[2]*(x-fittedParams[1])))+fittedParams[3],
		x0=countPoint, 
		dx=1.0, 
		n=2, 
		order=3)
		
		derivativeDF.loc[countPoint, "acceleration"] = RateDerivate
	
	
	threshold = derivativeDF["acceleration"].idxmax()
	#--------------------plot second derivative
	sns.lineplot(derivativeDF["countPoint"], derivativeDF["acceleration"], color='red', ax=axis[2])
	axis[2].set_title("Rates_acceleration "+ID)
	sns.scatterplot(x=[derivativeDF["acceleration"].idxmax()], y=[derivativeDF["acceleration"].max()] , s = 150, ax=axis[2], color = "black")
	#-------------------------------------------
	
	return derivativeDF["acceleration"].idxmax()








def barcodeClassifier_main(DBLmetricsDF, segmentationDF):
	
	mixedIDs = list(set(DBLmetricsDF["FirstID"]).union(set(DBLmetricsDF["SecondID"])))
	
	normCounts = normCenter(DBLmetricsDF, mixedIDs)
	
	positiveCellsID = {}
	
	for ID in mixedIDs:
		
		countsDF = singletscontribution(segmentationDF, normCounts, ID)
		
		posRatesDF = definePositiveRates(countsDF)
		
		threshold = findThreshold(posRatesDF, ID)
	
		positiveCellsID[ID] = pd.DataFrame({"IDs_in_barcode":""},index= normCounts.index)
		positiveCellsID[ID].loc[normCounts["ID_"+ID] >= threshold, "IDs_in_barcode"] = ID
	
	
	positiveCellsID = pd.concat(positiveCellsID.values(), axis = 1).apply(lambda x: ','.join(x[x.notnull()]), axis = 1)
	
	positiveCellsID=positiveCellsID.replace(',$','', regex=True)
	
	positiveCellsID=positiveCellsID.replace(',,',',', regex=True)
	
	positiveCellsID=positiveCellsID.replace('^,','', regex=True)
	
	return positiveCellsID