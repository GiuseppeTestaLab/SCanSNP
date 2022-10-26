#!/usr/bin/env python

#VCF preprocessing utils for mito


import itertools
import pandas as pd
import io
import re
from collections import Counter
import math
from SCanSNP.VCFUtils import *
import pysam
import numpy as np




def ExtractMitoPositions(bamFile, mitoContig):
	
	bam=pysam.AlignmentFile(bamFile, "rb")
	
	GenotypesDF = pd.DataFrame()
	
	posList = []
	for i in bam.pileup(mitoContig,max_depth =1):
		posList.append(i.reference_pos)
	
	if len(posList) != bam.get_reference_length(mitoContig):
		print("WARNING detected positions do not match chromoseome length!!!")
	
	#temp!!
	posList = posList[1:20]
	
	GenotypesDF["CHROM"] = [str(mitoContig)]*len(posList)
	GenotypesDF["POS"] = posList
	GenotypesDF.index = GenotypesDF["CHROM"]+"_"+GenotypesDF["POS"].astype(str)
	
	
	bam.close()
	return GenotypesDF
	
	



def Pileupper_noRef(bam, locus, BarcodeSet, readList, barcodetag, umitag, bannedFlag=3844, mapQuality=2, baseQuality=20 , readLength=30, maxReadsPerBC=10):
	##!!! Accessing readList "fake global" because of diverse GIL spawned for childs<<<<
	# UMItag not used in this version
	i = 0
	cblist = []
	BarcodeSetLocal = BarcodeSet.copy()	
	for read in bam.fetch(locus[0], locus[1]-1, locus[1]):
		# This will try to avoid endless counting on few hugely covered loci, its not optimal.
		if (i > len(BarcodeSetLocal)*maxReadsPerBC) or (len(BarcodeSetLocal) == 0):
			break
	#Check for position coverage
		try:
			CB=read.get_tag(barcodetag)
		except:
			continue
		if not set([CB]).intersection(BarcodeSetLocal):
			continue
		try:
			position = read.positions.index(int(locus[1]) - 1)
		except:
			continue
		if read.query_alignment_qualities[position] < 20:
			continue
		if read.mapq < 2:
			continue
		if read.rlen < 30:
			continue
		if read.flag == 3844:
			continue
		redBase = read.query_alignment_sequence[position]
		readList.append([locus[0]+ "_"+str(locus[1]), CB, redBase])
		cblist.append(CB)
		if cblist.count(CB) == maxReadsPerBC:
			BarcodeSetLocal.remove(CB)
		i = i+1
	return readList
	


def DFMaker_noRef(readList, barcodeList):
	ReadsList=pd.DataFrame(readList, columns = ["Pos","Barcode","Base"])
	Counts = ReadsList.groupby(["Pos","Barcode"])["Base"].apply(lambda x: "".join(set(np.unique(x)))).to_frame()
	Counts["Base"] =  np.where(Counts["Base"].str.len() == 2, Counts["Base"], Counts["Base"]*2)
	Counts = Counts.unstack().T.droplevel(0)
	#Remove columns with no variant loci
	#Counts = Counts[Counts.apply(lambda col: col.unique()).apply(lambda col: len(col))[Counts.apply(lambda col: col.unique()).apply(lambda col: len(col)) > 2].index.tolist()]
	return Counts