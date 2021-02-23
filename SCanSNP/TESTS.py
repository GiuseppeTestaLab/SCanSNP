def iterateFitting_Test2(DBLmetricsDF):
	'''
	Iterating over ranger of Second ID / Total droplet Private reads contributions
	To identify peaks in Doublets detection and best thresholding
	Providing DBLs count for each ID separately
	'''
	r=robjects.r
	r.source(fitdistPath)
	r.source(quantilePath)
	print("Iterative fitting negbin for different ID1/ID2 ratios")
	FittedIterations = pd.DataFrame(index = np.arange(0.1,0.61, 0.01).round(2), dtype="float32")
	for SecondIDfraction in np.arange(0.1,0.61, 0.01).round(2):
		Outliers="Trimmed"
		ReadSignals = pd.DataFrame(index = DBLmetricsDF.index, dtype="int")
		ID_specificThresholds = {}
		#Defining "negative" droplets for each ID
		for id in DBLmetricsDF.FirstID.unique():
			IDNoise = DBLmetricsDF.loc[(DBLmetricsDF["FirstID"] != id) & ~((DBLmetricsDF["SecondID"] == id) & (DBLmetricsDF["".join(["ID_",id])] >= DBLmetricsDF[["".join(["ID_", ident]) for ident in DBLmetricsDF.FirstID.unique() ]].sum(axis = 1)*SecondIDfraction)),"".join(["ID_", id])]
			#We try to trim top 5% of droplet signal before fitting negbin in order to remove outlier
			CroppedNoise = IDNoise[IDNoise < IDNoise.quantile(.95)]
			if CroppedNoise.sum() == 0:
				#In case trimming fails (result in all 0 values) we fall back on untrimmed values
				#storing the result in $Outliers variable for subsequent assessment
				print(".5 outliers trimming was not possible")
				Outliers="NotTrimmed"
				CroppedNoise = IDNoise
			#Wrapping fitdistrplus::fitdist() and stats::quantile() functions from R scripts
			distFitt=r.fitdist(robjects.vectors.FloatVector(CroppedNoise), distr = "nbinom")
			Threshold=r.quantile(distFitt, probs = .99)[0][0][0]
			ID_specificThresholds["".join([id,"_R.Threshold"])] = Threshold
			Positiveid = DBLmetricsDF[DBLmetricsDF["".join(["ID_",id])] > Threshold].index
			ReadSignals.loc[Positiveid,id] = 1
			ReadSignals=ReadSignals.fillna(0)
		#Subsetting for double (or more) positive droplets
		MultiSIg = ReadSignals[ReadSignals.sum(axis = 1) > 1]
		
		#Storing Total DBLs at the given ratio and ID spec  presence in DBLs
		FittedIterations.loc[SecondIDfraction,"DBLs_total"] = MultiSIg.shape[0]
		
		for ID in MultiSIg.columns:
			FittedIterations.loc[SecondIDfraction, ID+"_inDBLs"] = MultiSIg[ID].sum()
		
		FittedIterations.loc[SecondIDfraction,"Outliers"] = Outliers
		print("".join(["Completed ", str(round(SecondIDfraction,2)), " ID2 fraction iteration"]))
	return FittedIterations


def ID2_RatioSelect_test2(FittedIterations):
	'''
	For each iteration we store the number of lost doublets with respect of previous one to spot the highest "Loss"
	'''
	FittedIterationsLoss = FittedIterations.copy()
	for ID in [id for id in FittedIterations.columns if "_inDBLs" in id]:
		FittedIterationsLoss[ID+"_DBLsLoss"] = FittedIterationsLoss[ID].diff().fillna(0.0).abs()
	
	FittedIterationsLoss["DBLsLoss"] = FittedIterationsLoss["DBLs"].diff().fillna(0.0).abs()
	FittedIterationsLoss["DBLsLoss"] = FittedIterationsLoss["DBLs"].diff().fillna(0.0).abs()
	#If 5% outliers trimming failed we ignore the first peak because it will likely be due to non trimmed outliers retention
	if FittedIterationsLoss["Outliers"].unique()[0] == "NotTrimmed":
		FittedIterationsLoss <-  FittedIterationsLoss.iloc[2:]
		MaxLoss = FittedIterationsLoss["DBLsLoss"].idxmax()
		ID2_Ratio = FittedIterationsLoss.index[FittedIterationsLoss.index.get_loc(MaxLoss)-1]
	else:
		MaxLoss = FittedIterationsLoss["DBLsLoss"].idxmax()
		ID2_Ratio = FittedIterationsLoss.index[FittedIterationsLoss.index.get_loc(MaxLoss)-1]
	return ID2_Ratio



#### RMSD tests
#### RMSD tests
#### RMSD tests
#### RMSD tests

import itertools
import pandas as pd
from DBLsutils import *


chrLenDict = {"1":"248956422",
"2":"242193529",
"3":"198295559",
"4":"190214555",
"5":"181538259",
"6":"170805979",
"7":"159345973",
"8":"145138636",
"9":"138394717",
"10":"133797422",
"11":"135086622",
"12":"133275309",
"13":"114364328",
"14":"107043718",
"15":"101991189",
"16":"90338345",
"17":"83257441",
"18":"80373285",
"19":"58617616",
"20":"64444167",
"21":"46709983",
"22":"50818468",
"Y":"57227415",
"X":"156040895",
"MT":"16500"}




def chrScan_II(column):
	if column.count() > 1:
		combinations = list(itertools.combinations(column.dropna(), 2))
		combDiff = np.array([Pos[0] - Pos[1] for Pos in combinations])
		sqCombDiff = ( np.power(combDiff,2).sum())/(len(combDiff))
		return sqCombDiff
	else:
		return int(0)


def compute_RMSD(TotalDist, chrLenDict):
	baseLineDistsDF = pd.DataFrame.from_dict(chrLenDict, orient='index').loc[TotalDist.index].astype("int")
	baseLineDistsDF = pd.concat([baseLineDistsDF]*len(TotalDist.columns), axis = 1)
	baseLineDistsDF.columns = TotalDist.columns
	RMSD = TotalDist.subtract(baseLineDistsDF, axis = 0)
	RMSD = RMSD.pow(2)
	RMSD = (RMSD.sum()).divide(len(baseLineDistsDF))
	RMSD = np.sqrt(RMSD)
	return RMSD


def RMSDWRAPPER(HomozAltScore,HeterozAltScore, Doublet,DBLSpecificSingularLociDict, DropToDBLDict):
	pd_HomozAltScore = pd.DataFrame(HomozAltScore.todense(), index = DBLSpecificSingularLociDict["Homoz"][Doublet], columns=DropToDBLDict[Doublet])
	pd_HeterozAltScore = pd.DataFrame(HeterozAltScore.todense(), index = DBLSpecificSingularLociDict["HeteroAlt"][Doublet], columns=DropToDBLDict[Doublet])
	CoveredLoci=pd.concat([pd_HomozAltScore,pd_HeterozAltScore])
	CoveredLoci=CoveredLoci[CoveredLoci > 0]
	CoveredLoci["pos"] = CoveredLoci.index.to_series().str.split("_").str[1].astype(int)
	CoveredLoci = CoveredLoci.where(CoveredLoci.isnull(), CoveredLoci.pos, axis = 0)
	CoveredLoci["chr"] = CoveredLoci.index.to_series().str.split("_").str[0]
	CoveredLoci.drop("pos", axis = 1, inplace = True)
	grouped = CoveredLoci.groupby("chr")
	#Calculations on chromosome with more than 1 covered locus
	SNPdist=grouped[CoveredLoci.columns[:-1]].agg(chrScan_II)
	#Calculations on chromosome with less than 1 covered locus (will account for whole chr size)
	ChrDists = SNPdist.copy()
	ChrDists=ChrDists.replace([0,0.0],[np.nan,np.nan])
	ChrDists=ChrDists.apply(lambda x: np.where(x > 0,0.0,x))
	ChrDists["chrLen"] = [int(chrLenDict[str(x)]) for x in ChrDists.index]
	ChrDists = ChrDists.where(ChrDists.notnull(), ChrDists.chrLen, axis = 0)
	#Sum covered ant empty chromosomes and dividing by sqrt(n chromosomes)
	ChrDistsSum=ChrDists[ChrDists.columns[:-1]]
	TotalDist=(ChrDistsSum.add(SNPdist))
	RMSD = compute_RMSD(TotalDist, chrLenDict)
	return RMSD



DBLInfo = pd.DataFrame( columns = ["DBL_FirstID_Score","DBL_SecondID_Score"])
for Doublet in ExtractOrderedDoublets(vcf):
	
	comp1 = Doublet[0]
	comp2 = Doublet[1]
	
	#Slice for relevant Loci/Barcodes
	#Homoz Loci
	SRef_HOMOZ,SAlt_HOMOZ,GenotypesDF_sliced_HOMOZ=MultiSlice(SparseD, GenotypesDF, DBLSpecificSingularLociDict["Homoz"][Doublet], DropToDBLDict[Doublet])
	#Heteroz Loci Ref
	SRef_RefHet,SAlt_RefHet,GenotypesDF_sliced_RefHet=MultiSlice(SparseD, GenotypesDF, DBLSpecificSingularLociDict["HeteroRef"][Doublet], DropToDBLDict[Doublet])
	#Heteroz Loci Alt
	SRef_AltHet,SAlt_AltHet,GenotypesDF_sliced_AltHet=MultiSlice(SparseD, GenotypesDF, DBLSpecificSingularLociDict["HeteroAlt"][Doublet], DropToDBLDict[Doublet])
	# Locus-specific Noise calculation
	CleanRefSignal,CleanAltSignal,CleanRef_HET_Signal,CleanAlt_HET_Signal=LocusSpecNoisCalc(DropToDBLDict,
		Doublet,SingularLociScoreDF,GenotypesDF,vcf,
		DBLSpecificSingularLociDict["Homoz"][Doublet],
		DBLSpecificSingularLociDict["HeteroRef"][Doublet],
		DBLSpecificSingularLociDict["HeteroAlt"][Doublet])
		
	GenotypesDF_sliced_HOMOZ=GenotypesDF_sliced_HOMOZ.replace(1,2)
	GenotypesDF_sliced_RefHet=GenotypesDF_sliced_RefHet.replace(1,2)
	GenotypesDF_sliced_AltHet=GenotypesDF_sliced_AltHet.replace(1,2)
	##First thing use HomoLoci information
	##First thing use HomoLoci information
	##First thing use HomoLoci information
	# Cleaning Signals from Locus Nloise
	RefReads_HOMOZ=SRef_HOMOZ.multiply(CleanRefSignal)
	AltReads_HOMOZ=SAlt_HOMOZ.multiply(CleanAltSignal)
	#For each component of Putative DBl we calc the score
	HomozRefScoresComp1 = RefReads_HOMOZ.T.multiply(GenotypesDF_sliced_HOMOZ[comp1+"_RefAl"].to_numpy()).T
	HomozRefScoresComp2 = RefReads_HOMOZ.T.multiply(GenotypesDF_sliced_HOMOZ[comp2+"_RefAl"].to_numpy()).T
	HomozAltScoresComp1 = AltReads_HOMOZ.T.multiply(GenotypesDF_sliced_HOMOZ[comp1+"_AltAl"].to_numpy()).T
	HomozAltScoresComp2 = AltReads_HOMOZ.T.multiply(GenotypesDF_sliced_HOMOZ[comp2+"_AltAl"].to_numpy()).T
	HomozComp1Total = HomozRefScoresComp1+HomozAltScoresComp1
	HomozComp2Total = HomozRefScoresComp2+HomozAltScoresComp2
	
	#Now info about Singular Ref loci
	#Now info about Singular Ref loci
	#Now info about Singular Ref loci
	# Cleaning Signals from Locus Nloise
	RefReads_HETEROZ=SRef_RefHet.multiply(CleanRef_HET_Signal)
	HeterozRefScoreComp1 = RefReads_HETEROZ.T.multiply(GenotypesDF_sliced_RefHet[comp1+"_RefAl"].to_numpy()).T
	HeterozRefScoreComp2 = RefReads_HETEROZ.T.multiply(GenotypesDF_sliced_RefHet[comp2+"_RefAl"].to_numpy()).T
	
	#Now those of Singular Alt loci
	#Now those of Singular Alt loci
	#Now those of Singular Alt loci
	AltReads_HETEROZ=SAlt_AltHet.multiply(CleanAlt_HET_Signal)
	HeterozAltScoreComp1 = AltReads_HETEROZ.T.multiply(GenotypesDF_sliced_AltHet[comp1+"_AltAl"].to_numpy()).T
	HeterozAltScoreComp2 = AltReads_HETEROZ.T.multiply(GenotypesDF_sliced_AltHet[comp2+"_AltAl"].to_numpy()).T
	
	### COMPUTE TOTAL SCORES
	### COMPUTE TOTAL SCORES
	### COMPUTE TOTAL SCORES
	Comp1CumulativeScore = HomozComp1Total.sum(axis =0)+HeterozRefScoreComp1.sum(axis =0)+HeterozAltScoreComp1.sum(axis =0)
	Comp2CumulativeScore = HomozComp2Total.sum(axis =0)+HeterozRefScoreComp2.sum(axis =0)+HeterozAltScoreComp2.sum(axis =0)
	
	### COMPUTE RMSD
	### COMPUTE RMSD
	### COMPUTE RMSD
	TotalDistComp1=RMSDWRAPPER(HomozAltScoresComp1,HeterozAltScoreComp1, Doublet,DBLSpecificSingularLociDict, DropToDBLDict)
	TotalDistComp2=RMSDWRAPPER(HomozAltScoresComp2,HeterozAltScoreComp2, Doublet,DBLSpecificSingularLociDict, DropToDBLDict)
	
	DBLstatus=pd.concat([pd.DataFrame(Comp1CumulativeScore.T, index = DropToDBLDict[Doublet], columns = ["DBL_FirstID_Score"]),
		pd.DataFrame(Comp2CumulativeScore.T, index = DropToDBLDict[Doublet], columns = ["DBL_SecondID_Score"]),
		pd.DataFrame(TotalDistComp1, columns = ["TotalDistComp1"]),
		pd.DataFrame(TotalDistComp2, columns = ["TotalDistComp2"])], axis = 1)
		
	DBLInfo = pd.concat([DBLInfo,DBLstatus.round(2)], axis = 0)
	



##MODDED LOWQUALS


def AmbientFactorsEstimate_test2(DBLmetricsDF):
	'''
	Not removing DBLs here
	'''
	SNGsmetricsDF = DBLmetricsDF[DBLmetricsDF.DBL_FirstID_Score != 0]
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


def AddPseudoCounts_test2(DBLmetricsDF, AmbientFactors):
	NoisedIDs = DBLmetricsDF.copy()
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


AmbientFactors = AmbientFactorsEstimate_test2(DBLmetricsDF)

NoisedIDs = AddPseudoCounts_test2(DBLmetricsDF, AmbientFactors)

NoisedIDs["logFC"] = np.log(NoisedIDs["Noised_FirstID_Score"].divide(NoisedIDs["Noised_SecondID_Score"]))




test = pd.concat([DBLmetricsDF,DBLInfo[["TotalDistComp1","TotalDistComp2"]], NoisedIDs], axis = 1)
test.to_csv(outdir + "/test.tsv", sep = "\t", header = True, index = True)
