#!/usr/bin/env python

#VCF preprocessing utils for mito


import itertools
import pandas as pd
import io
import re
from collections import Counter
import math
from SCanSNP.VCFUtils import *






def ExtractMitoPositions(bamFile, mitoContig):
	
	bam=pysam.AlignmentFile(bamFile, "rb")
	
	GenotypesDF = pd.DataFrame()
	
	posList = []
	for i in bam.pileup(mitoContig,max_depth =1):
		posList.append(i.reference_pos)
	
	if len(posList) != bam.get_reference_length(mitoContig):
		print("WARNING detected positions do not match chromoseome length!!!")
	
	GenotypesDF["CHROM"] = [str(mitoContig)]*len(posList)
	GenotypesDF["POS"] = posList
	GenotypesDF.index = GenotypesDF["CHROM"]+"_"+GenotypesDF["POS"].astype(str)
	
	
	bam.close()
	return GenotypesDF
	
	



def Pileupper_noRef(bam, locus, BarcodeSet, readList, barcodetag, umitag, bannedFlag=3844, mapQuality=2, baseQuality=20 , readLength=30):
	##!!! Accessing readList "fake global" because of diverse GIL spawned for childs<<<<
	# UMItag not used in this version
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
		readList.append([locus[0]+ "_"+str(locus[1]), CB, redBase])
	return readList
	

def DFMaker_noRef(readList, barcodeList, locus):
	locusIndex = str(locus[0])+"_"+str(locus[1])
	ReadsList=pd.DataFrame(readList, columns = ["Pos","Barcode","Base"])[["Barcode","Base"]]
	Counts = ReadsList.groupby(["Barcode"])["Base"].apply(lambda x: "".join(set(np.unique(x)))).to_frame(name=locusIndex)
	Counts[locusIndex] =  np.where(Counts[locusIndex].str.len() == 2, Counts[locusIndex], Counts[locusIndex]*2)
	return Counts