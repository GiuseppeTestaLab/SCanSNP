
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







def ReadCounter(chunkIdx, bamFile, barcodeList, GenotypeChunkList):
	bam=pysam.AlignmentFile(bamFile, "rb")
	BarcodeSet=set(barcodeList)
	readList = []
	SparseD = {"sparse_Ref" : csr_matrix((0, len(barcodeList))), "sparse_Alt":csr_matrix((0, len(barcodeList))), "Locus" : pd.Series()  }
	lastPos = chunkIdx[-1]
	for indexPos in chunkIdx:
		
		locus = GenotypeChunkList[indexPos]
		
		readList=Pileupper(bam, locus, lastPos,BarcodeSet, readList)
		
		if sys.getsizeof(readList) >= 30000000 or indexPos == lastPos:
			try:
				sparse_Ref_TMP,sparse_Alt_TMP,locus_TMP = DFMaker(readList, barcodeList)
			except:
				continue
			SparseD["sparse_Ref"] = scipy.sparse.vstack((SparseD["sparse_Ref"],sparse_Ref_TMP))
			SparseD["sparse_Alt"] = scipy.sparse.vstack((SparseD["sparse_Alt"],sparse_Alt_TMP))
			SparseD["Locus"] = SparseD["Locus"].append(locus_TMP)
			readList = []
	bam.close()
	print("Chunk "+ str(chunkIdx) + " completed")
	return SparseD



def Pileupper(bam, locus, lastPos, BarcodeSet, readList, bannedFlag=3844, mapQuality=2, baseQuality=20 , readLength=30):
	##!!! Accessing readList "fake global" because of diverse GIL spawned for childs<<<<
	for read in bam.fetch(locus[0], locus[1]-1, locus[1]):
	#Check for position coverage
		try:
			CB=read.get_tag("CB")
		except:
			continue
		if not set([CB]).intersection(BarcodeSet):
			continue
		try:
			position = read.positions.index(int(locus[1]) - 1)
		except:
			continue
		try:
			UB=read.get_tag("UB")
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


