#!/usr/bin/env python

#Main function scanning BAM

import pysam
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import itertools
from scipy.sparse import csr_matrix

sys.path.append('./SCanSNPv3_REWORKEDLIST/')







# def ReadCounter2(chunkId, bamFile, barcodeList, GenotypeChunkDict, bannedFlag=3844, mapQuality=2, baseQuality=20 , readLength=30):
# 	bam=pysam.AlignmentFile(bamFile, "rb")
# 	chunk=GenotypeChunkDict[chunkId]
# 	RegionList=chunk.loc[:,["CHROM","POS","REF","ALT"]]
# 	RegionList["POS"] =RegionList["POS"].astype(int)
# 	RegionList=RegionList.values.tolist()
# 	readList=[]
# 	for region in RegionList:
# 		UsedUMI = []
# 		for pileupcolumn in bam.pileup(contig=region[0], start=region[1]-1,stop=region[1], flag_filter = bannedFlag, min_mapping_quality = 2,truncate=True):
# 			for read in pileupcolumn.pileups:
# 				if (region[1]-1 in  read.alignment.get_reference_positions()
# 				and read.alignment.has_tag("CB")
# 				and read.alignment.has_tag("UB")
# 				and read.alignment.rlen > readLength
# 				and read.alignment.get_tag('CB') in barcodeList
# 				and read.alignment.get_tag('CB') + "_" + read.alignment.get_tag('UB') not in UsedUMI
# 				and ( read.alignment.query_alignment_sequence[read.alignment.positions.index(int(region[1]) - 1)] == region[2] or
# 				read.alignment.query_alignment_sequence[read.alignment.positions.index(int(region[1]) - 1)] == region[3])
# 				and read.alignment.query_alignment_qualities[read.alignment.positions.index(int(region[1]) - 1)] >= baseQuality
# 				and not read.alignment.is_secondary and not read.alignment.is_unmapped):
# 					Base = read.alignment.query_alignment_sequence[read.alignment.positions.index(int(region[1]) - 1)]
# 					readList.append([region[0]+ "_"+str(region[1]), read.alignment.get_tag('CB'), Base, region[2], region[3] ])
# 					UsedUMI.append(read.alignment.get_tag('CB') + "_" + read.alignment.get_tag('UB'))
# 	bam.close()
# 	print("Chunk number "+ str(chunkId) + " completed")
# 	return readList
#
#
#
# def ReadCounter3(chunkId, bamFile, barcodeList, GenotypeChunkDict, bannedFlag=3844, mapQuality=2, baseQuality=20 , readLength=30):
# 	bam=pysam.AlignmentFile(bamFile, "rb")
# 	chunk=GenotypeChunkDict[chunkId]
# 	locilist=chunk.loc[:,["CHROM","POS","REF","ALT"]]
# 	locilist["POS"] =locilist["POS"].astype(int)
# 	locilist=locilist.values.tolist()
# 	readList=[]
# 	for locus in locilist:
# 		UsedUMI = []
# 		for read in bam.fetch(locus[0], locus[1]-1, locus[1]):
# 		#Check for position coverage
# 			if locus[1]-1 in read.get_reference_positions():
# 				baseID = read.positions.index(locus[1]-1)
# 				#Reads filtering
# 				if (read.flag == bannedFlag or read.mapq < mapQuality or read.rlen < readLength or read.query_alignment_qualities[baseID] < baseQuality ):
# 					continue
# 				if not read.has_tag("CB") or read.is_secondary or read.is_unmapped or not read.has_tag("UB") :
# 					continue
# 				if read.get_tag('CB') not in barcodeList:
# 					continue
# 				UMI = read.get_tag('CB') + "_" + read.get_tag('UB')
# 				if UMI in UsedUMI:
# 					continue
# 				#AltCounts
# 				CBtag = read.get_tag('CB')
# 				base = read.query_alignment_sequence[baseID].upper()
# 				if (read.query_alignment_sequence[read.positions.index(locus[1]-1)] == locus[2]
# 				or read.query_alignment_sequence[read.positions.index(locus[1]-1)] == locus[3]):
# 					readList.append([locus[0]+ "_"+str(locus[1]), read.get_tag('CB'), base, locus[2], locus[3] ])
# 					UsedUMI.append(UMI)
# 	bam.close()
# 	print("Chunk number "+ str(chunkId) + " completed")
# 	return readList


# def ReadCounter(chunkId, bamFile, barcodeList, GenotypeChunkDict, bannedFlag=3844, mapQuality=2, baseQuality=20 , readLength=30):
# 	bam=pysam.AlignmentFile(bamFile, "rb")
# 	chunk=GenotypeChunkDict[chunkId]
# 	locilist=chunk.loc[:,["CHROM","POS","REF","ALT"]]
# 	locilist["POS"] =locilist["POS"].astype(int)
# 	locilist=locilist.values.tolist()
# 	readList=[]
# 	for locus in locilist:
# 		UsedUMI = []
# 		for read in bam.fetch(locus[0], locus[1]-1, locus[1]):
# 		#Check for position coverage
# 			if (locus[1]-1 in  read.get_reference_positions()
# 			and read.has_tag("CB")
# 			and read.mapq >= mapQuality
# 			and read.flag != bannedFlag
# 			and read.has_tag("UB")
# 			and read.rlen > readLength
# 			and read.get_tag('CB') in barcodeList
# 			and read.get_tag('CB') + "_" + read.get_tag('UB') not in UsedUMI
# 			and read.query_alignment_qualities[read.positions.index(int(locus[1]) - 1)] >= baseQuality
# 			and not read.is_secondary and not read.is_unmapped):
# 				UsedUMI.append(read.get_tag('CB') + "_" + read.get_tag('UB'))
# 				if ( read.query_alignment_sequence[read.positions.index(int(locus[1]) - 1)] == locus[2] or
# 				read.query_alignment_sequence[read.positions.index(int(locus[1]) - 1)] == locus[3]):
# 					Base = read.query_alignment_sequence[read.positions.index(int(locus[1]) - 1)]
# 					readList.append([locus[0]+ "_"+str(locus[1]), read.get_tag('CB'), Base, locus[2], locus[3] ])
# 	bam.close()
# 	print("Chunk number "+ str(chunkId) + " completed")
# 	return readList


def ReadCounter(chunkIdx, bamFile, barcodeList, bannedFlag=3844, mapQuality=2, baseQuality=20 , readLength=30):
	bam=pysam.AlignmentFile(bamFile, "rb")
	BarcodeSet=set(barcodeList)
	readList=[]
	for indexPos in chunkIdx:
		locus = GenotypeChunkList[indexPos]
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
	bam.close()
	print("Chunk "+ str(chunkIdx) + " completed")
	return readList



def DFMaker(readList, barcodeList):
	ReadsList=pd.DataFrame(readList, columns = ["Pos","Barcode","Base","Ref","Alt"])
	uniqueLoci = list(set(ReadsList.Pos))
	RefReads=ReadsList[ReadsList["Base"] == ReadsList["Ref"]][["Pos","Barcode"]]
	AltReads=ReadsList[ReadsList["Base"] == ReadsList["Alt"]][["Pos","Barcode"]]
	del ReadsList
	#Storing Sparse categories
	Pos_c = pd.CategoricalDtype(sorted(uniqueLoci), ordered=True)
	Barcode_c = pd.CategoricalDtype(sorted(list(barcodeList)), ordered=True)
	#Storing Ref Counts
	#Storing Ref Counts
	print("Storing RefReads...")
	RefReads=RefReads.groupby(["Pos","Barcode"], as_index = False)
	RefReads=RefReads.size()
	RefReads.rename(columns={'size':'Counts'}, inplace=True)
	RefReads.index = list(RefReads.Pos)
	Pos_Ref = RefReads.Pos.astype(Pos_c).cat.codes
	Barcode_Ref = RefReads.Barcode.astype(Barcode_c).cat.codes
	sparse_Ref = csr_matrix((RefReads["Counts"], (Pos_Ref, Barcode_Ref)), shape=(len(uniqueLoci), len(barcodeList)))
	print("RefReads stored")
	#Storing Alt Counts
	#Storing Alt Counts
	print("Storing AltReads...")
	AltReads=AltReads.groupby(["Pos","Barcode"], as_index = False)
	AltReads=AltReads.size()
	AltReads.rename(columns={'size':'Counts'}, inplace=True)
	AltReads.index = list(AltReads.Pos)
	Pos_Alt = AltReads.Pos.astype(Pos_c).cat.codes
	Barcode_Alt = AltReads.Barcode.astype(Barcode_c).cat.codes
	sparse_Alt = csr_matrix((AltReads["Counts"], (Pos_Alt, Barcode_Alt)), shape=(len(uniqueLoci), len(barcodeList)))
	print("RefReads and AltReads stored")
	return sparse_Ref,sparse_Alt,pd.Series(sorted(uniqueLoci)),pd.Series(sorted(list(barcodeList)))
