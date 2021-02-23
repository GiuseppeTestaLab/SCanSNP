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
from dblsMark import *


def CountsMatrices(CleanSingularLoci, cleanLoci, MildcleanLoci, GenotypesDF, barcodeList,vcf,nThreads, bamFile):
	'''
	Fire-up the main Pileupper
	'''

	OmniIndex = list(set(CleanSingularLoci + DiffOnlyIndexMaker(vcf, cleanLoci) + DoubletSpecificSingularLociScan(vcf,GenotypesDF, MildcleanLoci)[0]))

	print('Splitting Variants into chunks...')
	GenotypeChunkList,GenotypeChunkIndexesList = FlattenDict(ChunkMaker(GenotypesDF, nThreads, OmniIndex))

	print('Pileup started...')
	start = time.time()

	#FireUp main Pileupper
	results = []
	pool=Pool(nThreads)

	for chunkIdx in GenotypeChunkIndexesList:
		result = pool.apply_async(ReadCounter, (chunkIdx, bamFile, barcodeList,GenotypeChunkList))
		results.append(result)


	pool.close()
	pool.join()

	print('Pileup took', time.time()-start, 'seconds.')

	#Merging pooled process results
	Reads=list(itertools.chain.from_iterable(([result.get() for result in results])))

	#load reads into sparse matrices
	SparseD = {}
	SparseD["sparse_Ref"],SparseD["sparse_Alt"],SparseD["Locus"],SparseD["Barcode"] = DFMaker(Reads, barcodeList)

	return SparseD



def deconvolution(SparseD, vcf, GenotypesDF, outdir, LowQual):
	#CountDF reconstruction
	# SparseCounts=scipy.sparse.load_npz(str(countpath) + '/Counts.npz')
	# Barcodes=pd.read_csv(str(countpath) + '/countBarcodes.tsv', header=None, names=["barcodes"])
	# Loci=pd.read_csv(str(countpath) + '/countLoci.tsv', header=None, names=["loci"])
	
	barcodeList = SparseD["Barcode"]
	
	#Pre-cleaning loci CHECK!!
	#Pre-cleaning loci CHECK!!
	cleanLoci = list(set(LociPreClean(vcf)).intersection(list(SparseD["Locus"])))
	MildcleanLoci = list(set(LociPreClean_milds(vcf)).intersection(list(SparseD["Locus"])))
	
	#Setting Loci categories
	DiffOnlyIndex=DiffOnlyIndexMaker(vcf, cleanLoci)
	CleanSingularLoci,SingularLoci_Alt,SingularLoci_Ref = SingularLociSCan(vcf,cleanLoci)
	
	DBLSpecificSingularLociDict = {}
	DBLSpecificSingularLociDict["Homoz"]=DoubletSpecificSingularLociScan(vcf,GenotypesDF, cleanLoci)[1]
	DBLSpecificSingularLociDict["HeteroAlt"]=DoubletSpecificSingularLociScan(vcf,GenotypesDF, cleanLoci)[2]
	DBLSpecificSingularLociDict["HeteroRef"]=DoubletSpecificSingularLociScan(vcf,GenotypesDF, cleanLoci)[3]
	
	
	LikeliHoodsDF=ComputeLikelihood(vcf,GenotypesDF,SparseD,DiffOnlyIndex)
	
	
	BestBarcodeID=LikeliHoodsDF.idxmax(axis = 1)
	
	BestInDropDict = defaultdict(list)
	for key, value in BestBarcodeID.to_dict().items():
		BestInDropDict[value].append(key)
		
	
	SingularLociScoreDF = SingularLociCNTR( SingularLoci_Alt, SingularLoci_Ref, SparseD, barcodeList, GenotypesDF,vcf)
	DBLsDF=NoiseRregression(BestInDropDict, barcodeList, vcf, SingularLociScoreDF, BestBarcodeID)
	
	#Creation of putative DBL-to-barcode dict
	DropToDBLDict = defaultdict(list)
	Barcodes_values = list(DBLsDF.index)
	DBLtuples_keys=list(zip(DBLsDF["FirstID"],DBLsDF["SecondID"]))
	for key, value in dict(zip(Barcodes_values, DBLtuples_keys)).items():
		DropToDBLDict[value].append(key)
	
	DBLs_ADJ_Contributions = LowQualScore(DropToDBLDict,SparseD,DBLSpecificSingularLociDict,vcf, GenotypesDF,SingularLociScoreDF)

	#Editing DFs before concat
	Contributions = SingularLociScoreDF[ExtractSamples(vcf)]
	BestIDs = DBLsDF.replace("ID_","", regex = True)

	DBLmetricsDF = pd.concat([Contributions,DBLs_ADJ_Contributions,BestIDs], axis = 1)

	#DBLs detection Module
	DBLsList = main__DBLsMark(DBLmetricsDF)
	
	DBLmetricsDF["DropletType"] = "Singlet"
	DBLmetricsDF.loc[DBLsList,"DropletType"] = "Doublet"
	DBLmetricsDF["ID"] = DBLmetricsDF["FirstID"]
	DBLmetricsDF.loc[DBLmetricsDF["DropletType"] == "Doublet","ID"] = "Doublet"

	if LowQual == True:
		QualDF = main_FlagLowQual(DBLmetricsDF, DBLsList, outdir)
		Cell_IDs = pd.concat([DBLmetricsDF, QualDF], axis = 1)
	else:
		Cell_IDs = DBLmetricsDF
		
		
	
	#pd.concat([Contributions,DBLsContributions,BestIDs], axis = 1).to_csv(writePath + "/DBLmetricsDF.tsv", sep = "\t", header = True, index = True)
	return Cell_IDs
