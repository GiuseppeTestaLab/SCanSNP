#!/usr/bin/env python

from VCFUtils import *
from DBLsutils import *
from Pileup import *
from scipy.sparse import csr_matrix
import time
from multiprocessing import Pool
import itertools
from ComputeLLK import *
from GenUtils import *
from LowQualutils import *
from lowQualityMark import *
from lowQualityMark_wEmpty import *
from dblsMark import *
from dblsMark_wEmpty import *
from itertools import chain
import scipy.sparse




def deconvolution(Counts, vcf, GenotypesDF, outdir, FullDrops, FullDropsKNNseries, platform):
	#CountDF reconstruction
	# SparseCounts=scipy.sparse.load_npz(str(countpath) + '/Counts.npz')
	# Barcodes=pd.read_csv(str(countpath) + '/countBarcodes.tsv', header=None, names=["barcodes"])
	# Loci=pd.read_csv(str(countpath) + '/countLoci.tsv', header=None, names=["loci"])
	
	emptyTraining = True if len(FullDrops) < len(Counts.barcodes) else False
	
	barcodeList = Counts.barcodes
	
	#Pre-cleaning loci CHECK!!
	#Pre-cleaning loci CHECK!!
	cleanLoci = list(set(LociPreClean(vcf)).intersection(list(Counts.loci)))
	MildcleanLoci = list(set(LociPreClean_milds(vcf)).intersection(list(Counts.loci)))
	
	#Setting Loci categories
	DiffOnlyIndex=DiffOnlyIndexMaker(vcf, cleanLoci)
	CleanSingularLoci,SingularLoci_Alt,SingularLoci_Ref = SingularLociSCan(vcf,cleanLoci)
	
	DBLSpecificSingularLociDict = {}
	DBLSpecificSingularLociDict["Homoz"]=DoubletSpecificSingularLociScan(vcf,GenotypesDF, cleanLoci)[1]
	DBLSpecificSingularLociDict["HeteroAlt"]=DoubletSpecificSingularLociScan(vcf,GenotypesDF, cleanLoci)[2]
	DBLSpecificSingularLociDict["HeteroRef"]=DoubletSpecificSingularLociScan(vcf,GenotypesDF, cleanLoci)[3]
	
	
	LikeliHoodsDF=ComputeLikelihood(vcf,GenotypesDF,Counts,DiffOnlyIndex)
	
	
	BestBarcodeID=LikeliHoodsDF.idxmax(axis = 1)
	
	BestInDropDict = defaultdict(list)
	for key, value in BestBarcodeID.to_dict().items():
		BestInDropDict[value].append(key)
	
	
	SingularLociScoreDF = SingularLociCNTR( SingularLoci_Alt, SingularLoci_Ref, Counts, barcodeList, GenotypesDF,vcf)
	
	if len(ExtractSamples(vcf)) > 2:
		DBLsDF=NoiseRregression(BestInDropDict, barcodeList, vcf, SingularLociScoreDF, BestBarcodeID, LikeliHoodsDF)
	else:
		ID1 = list(BestInDropDict.keys())[0]
		ID2 = list(BestInDropDict.keys())[1]
		DBLsDF = pd.DataFrame(columns = ["FirstID","SecondID"], index = list(chain(*BestInDropDict.values()))  )
		DBLsDF.loc[BestInDropDict[ID1], ["FirstID","SecondID"]] = [ID1,ID2]
		DBLsDF.loc[BestInDropDict[ID2], ["FirstID","SecondID"]] = [ID2,ID1]
	
	
	#Creation of putative DBL-to-barcode dict
	DropToDBLDict = defaultdict(list)
	Barcodes_values = list(DBLsDF.index)
	DBLtuples_keys=list(zip(DBLsDF["FirstID"],DBLsDF["SecondID"]))
	for key, value in dict(zip(Barcodes_values, DBLtuples_keys)).items():
		DropToDBLDict[value].append(key)
	
	DBLs_ADJ_Contributions = LowQualScore(DropToDBLDict,Counts,DBLSpecificSingularLociDict,vcf, GenotypesDF,SingularLociScoreDF)
	
	#Editing DFs before concat
	Contributions = SingularLociScoreDF[ExtractSamples(vcf)]
	BestIDs = DBLsDF.replace("ID_","", regex = True)
	
	DBLmetricsDF = pd.concat([Contributions,DBLs_ADJ_Contributions,BestIDs], axis = 1)
	
	if platform == "chromium":
		#DBLs detection Module
		DBLsList = main__DBLsMark_wEmpty(DBLmetricsDF,FullDrops ) if emptyTraining else main__DBLsMark(DBLmetricsDF)
		DBLmetricsDF["DropletType"] = "Singlet"
		DBLmetricsDF.loc[DBLsList,"DropletType"] = "Doublet"
		DBLmetricsDF["ID"] = DBLmetricsDF["FirstID"]
		DBLmetricsDF.loc[DBLmetricsDF["DropletType"] == "Doublet","ID"] = "Doublet"
		if emptyTraining:
			Cell_IDs = main_FlagLowQual_wEmpty(DBLmetricsDF, DBLsList, outdir,FullDrops,FullDropsKNNseries, HighOutlierThreshold = .95,LowOutlierThreshold = .05)
		else:
			Cell_IDs = main_FlagLowQual(DBLmetricsDF, DBLsList, outdir)	
		#pd.concat([Contributions,DBLsContributions,BestIDs], axis = 1).to_csv(writePath + "/DBLmetricsDF.tsv", sep = "\t", header = True, index = True)
	elif platform == "visium":
		Cell_IDs = DBLmetricsDF
		
	return Cell_IDs
