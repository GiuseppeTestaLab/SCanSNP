
#!/usr/bin/env python

#Main function scanning BAM

import pysam
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import itertools
from scipy.sparse import csr_matrix
import scipy.sparse
from SCanSNP.GenUtils import CountData

from SCanSNP.VCFUtils import *
from SCanSNP.DBLsutils import *

import time
from multiprocessing import Pool
from SCanSNP.ComputeLLK import *
#from SCanSNP.GenUtils import *
from SCanSNP.LowQualutils import *
from SCanSNP.lowQualityMark import *
from SCanSNP.lowQualityMark_wEmpty import *
from SCanSNP.dblsMark import *
from SCanSNP.dblsMark_wEmpty import *
from itertools import chain





def ReadCounter(chunkIdx, bamFile, barcodeList, GenotypeChunkList, barcodetag, umitag, filterbam):
	bam=pysam.AlignmentFile(bamFile, "rb")
	BarcodeSet=set(barcodeList)
	readList = []
	Counts = {"sparse_Ref" : csr_matrix((0, len(barcodeList)),dtype=np.int16), "sparse_Alt":csr_matrix((0, len(barcodeList)),dtype=np.int16), "Locus" : pd.Series(dtype="str") }
	lastPos = chunkIdx[-1]
	for indexPos in chunkIdx:
		
		locus = GenotypeChunkList[indexPos]
		
		readList = Pileupper_preFilt(bam, locus, lastPos,BarcodeSet, readList,barcodetag, umitag) if filterbam else Pileupper(bam, locus, lastPos,BarcodeSet, readList,barcodetag, umitag) 
		
		if sys.getsizeof(readList) >= 30000000 or indexPos == lastPos:
			try:
				sparse_Ref_TMP,sparse_Alt_TMP,locus_TMP = DFMaker(readList, barcodeList)
			except:
				continue
			Counts["sparse_Ref"] = scipy.sparse.vstack((Counts["sparse_Ref"],sparse_Ref_TMP))
			Counts["sparse_Alt"] = scipy.sparse.vstack((Counts["sparse_Alt"],sparse_Alt_TMP))
			Counts["Locus"] = pd.concat([Counts["Locus"], locus_TMP])
			readList = []
	bam.close()
	print("Chunk "+ str(chunkIdx) + " completed")
	return Counts



def Pileupper(bam, locus, lastPos, BarcodeSet, readList, barcodetag, umitag, bannedFlag=3844, mapQuality=2, baseQuality=20 , readLength=30):
	for read in bam.fetch(locus[0], locus[1]-1, locus[1]):
	#Check for position coverage
		try:
			CB=read.get_tag(barcodetag)
		except:
			continue
		if not set([CB]).intersection(BarcodeSet):
			continue
		try:
			position = read.positions.index(int(locus[1]) - 1)
		except:
			continue
		if read.query_alignment_qualities[position] < 20:
			continue
		if read.is_secondary:
			continue
		if read.mapq < 2:
			continue
		if read.rlen < 30:
			continue
		if read.flag == 3844:
			continue
		redBase = read.query_alignment_sequence[position]
		if ( redBase == locus[2] or redBase == locus[3]):
			readList.append([locus[0]+ "_"+str(locus[1]), CB, redBase, locus[2], locus[3] ])
	return readList


def Pileupper_preFilt(bam, locus, lastPos, BarcodeSet, readList, barcodetag, umitag, baseQuality=20 , readLength=30):
	for read in bam.fetch(locus[0], locus[1]-1, locus[1]):
	#Check for position coverage
		CB=read.get_tag(barcodetag)
		try:
			position = read.positions.index(int(locus[1]) - 1)
		except:
			continue
		if read.query_alignment_qualities[position] < baseQuality:
			continue
		if read.rlen < readLength:
			continue
		redBase = read.query_alignment_sequence[position]
		if ( redBase == locus[2] or redBase == locus[3]):
			readList.append([locus[0]+ "_"+str(locus[1]), CB, redBase, locus[2], locus[3] ])
	return readList



def DFMaker(readList, barcodeList):
	ReadsList=pd.DataFrame(readList, columns = ["Pos","Barcode","Base","Ref","Alt"])
	uniqueLoci = list(set(ReadsList.Pos))
	RefReads=ReadsList.loc[ReadsList["Base"] == ReadsList["Ref"], ["Pos","Barcode"]]
	AltReads=ReadsList.loc[ReadsList["Base"] == ReadsList["Alt"], ["Pos","Barcode"]]
	del ReadsList
	#Storing Sparse categories
	Pos_c = pd.CategoricalDtype(sorted(uniqueLoci), ordered=True)
	Barcode_c = pd.CategoricalDtype(sorted(list(barcodeList)), ordered=True)
	#Storing Ref Counts
	#Storing Ref Counts
	RefReads=RefReads.groupby(["Pos","Barcode"], as_index = False)
	RefReads=RefReads.size()
	RefReads.rename(columns={'size':'Counts'}, inplace=True)
	RefReads.index = list(RefReads.Pos)
	Pos_Ref = RefReads.Pos.astype(Pos_c).cat.codes
	Barcode_Ref = RefReads.Barcode.astype(Barcode_c).cat.codes
	sparse_Ref = csr_matrix((RefReads["Counts"], (Pos_Ref, Barcode_Ref)), shape=(len(uniqueLoci), len(barcodeList)))
	#Storing Alt Counts
	#Storing Alt Counts
	AltReads=AltReads.groupby(["Pos","Barcode"], as_index = False)
	AltReads=AltReads.size()
	AltReads.rename(columns={'size':'Counts'}, inplace=True)
	AltReads.index = list(AltReads.Pos)
	Pos_Alt = AltReads.Pos.astype(Pos_c).cat.codes
	Barcode_Alt = AltReads.Barcode.astype(Barcode_c).cat.codes
	sparse_Alt = csr_matrix((AltReads["Counts"], (Pos_Alt, Barcode_Alt)), shape=(len(uniqueLoci), len(barcodeList)))
	return sparse_Ref,sparse_Alt,pd.Series(sorted(uniqueLoci))




def CountsMatrices(CleanSingularLoci, cleanLoci, MildcleanLoci, GenotypesDF, barcodeList,vcf,nThreads, bamFile, barcodetag, umitag, filterbam):
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
		result = pool.apply_async(ReadCounter, (chunkIdx, bamFile, barcodeList,GenotypeChunkList,barcodetag, umitag,filterbam ))
		results.append(result)


	pool.close()
	pool.join()

	print('Pileup took', time.time()-start, 'seconds.')
	
	#Gathering rsults
	CountsDict = {"sparse_Ref" : csr_matrix((0, len(barcodeList))), "sparse_Alt":csr_matrix((0, len(barcodeList))), "Locus" : pd.Series(dtype="str"),  "Barcode" : pd.Series(sorted(list(barcodeList))) }
	for result in results:
		CountsDict["sparse_Ref"] = scipy.sparse.vstack((CountsDict["sparse_Ref"],result.get()["sparse_Ref"]))
		CountsDict["sparse_Alt"] = scipy.sparse.vstack((CountsDict["sparse_Alt"],result.get()["sparse_Alt"]))
		CountsDict["Locus"] =  pd.concat([CountsDict["Locus"],result.get()["Locus"]])
	
	Counts = CountData(CountsDict["sparse_Ref"], CountsDict["sparse_Alt"], CountsDict["Locus"], CountsDict["Barcode"])
	
	return Counts
	
	


def CountsPileup(MildcleanLoci, GenotypesDF, barcodeList,vcf,nThreads, bamFile,barcodetag, umitag, filterbam):
	'''
	Fire-up the main Pileupper assuming 1 only ID in VCF i.e. sc pileup over list of loci
	'''

	print('Splitting Variants into chunks...')
	GenotypeChunkList,GenotypeChunkIndexesList = FlattenDict(ChunkMaker(GenotypesDF, nThreads, MildcleanLoci))

	print('Pileup started...')
	start = time.time()

	#FireUp main Pileupper
	results = []
	pool=Pool(nThreads)

	for chunkIdx in GenotypeChunkIndexesList:
		result = pool.apply_async(ReadCounter, (chunkIdx, bamFile, barcodeList,GenotypeChunkList, barcodetag, umitag, filterbam))
		results.append(result)


	pool.close()
	pool.join()

	print('Pileup took', time.time()-start, 'seconds.')
	
	#Gathering rsults
	CountsDict = {"sparse_Ref" : csr_matrix((0, len(barcodeList))), "sparse_Alt":csr_matrix((0, len(barcodeList))), "Locus" : pd.Series(dtype="str"),  "Barcode" : pd.Series(sorted(list(barcodeList))) }
	for result in results:
		CountsDict["sparse_Ref"] = scipy.sparse.vstack((CountsDict["sparse_Ref"],result.get()["sparse_Ref"]))
		CountsDict["sparse_Alt"] = scipy.sparse.vstack((CountsDict["sparse_Alt"],result.get()["sparse_Alt"]))
		CountsDict["Locus"] =  pd.concat([CountsDict["Locus"],result.get()["Locus"]])
	
	Counts = CountData(CountsDict["sparse_Ref"], CountsDict["sparse_Alt"], CountsDict["Locus"], CountsDict["Barcode"])
	
	return Counts