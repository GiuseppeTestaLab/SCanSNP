#LinkageDisequilibrium test

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


def LowQualScore_RMSD(DropToDBLDict,SparseD,DBLSpecificSingularLociDict,vcf, GenotypesDF,SingularLociScoreDF):
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
	
