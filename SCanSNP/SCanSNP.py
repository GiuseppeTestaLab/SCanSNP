
#!/usr/bin/env python

import argparse
import os




writePath="/hpcnfs/home/ieo4777/temp"
bamFile="/hpcnfs/scratch/temporary/Dav_vc/0_CellrangerAlign/Sample_S20272_157/outs/possorted_genome_bam.bam"
vcf="/hpcnfs/scratch/temporary/Dav_vc/2_Combining/Sample_S20273_157/output/combinedVCF/jointGenotype.g.vcf"
barcodesFILE="/hpcnfs/home/ieo4777/temp/barcodes.tsv"
writePath="/hpcnfs/home/ieo4777/temp"
countpath="/hpcnfs/home/ieo4777/temp"

nThreads=15
# outdir=args.outdir
# countpath=args.countpath





import pysam
import pandas as pd
import time
import itertools
import io
import numpy as np
import sys
from collections import defaultdict
import re
import math
from collections import Counter
from sklearn.linear_model import LogisticRegression
import scipy
import pickle


sys.path.append('/hpcnfs/home/ieo4777/gitted/SCanSNP_v1/SCanSNP')

from VCFUtils import *
from Wrappers import *


if __name__ == "__main__":
	#Creation of differen loci subsets
	cleanLoci = LociPreClean(vcf)
	MildcleanLoci = LociPreClean_milds(vcf)
	CleanSingularLoci,SingularLoci_Alt,SingularLoci_Ref = SingularLociSCan(vcf,cleanLoci)
	
	barcodeList=pd.read_csv(barcodesFILE, header=None, names=["b"])["b"].tolist()
	
	#Genotypes map creation
	GenotypesDF = creategenotypeDF(vcf)
	
	
	if mode == "matrixgen":
		SparseD = CountsMatrices(CleanSingularLoci, cleanLoci,
			MildcleanLoci, GenotypesDF, 
			barcodeList, vcf, 
			nThreads, bamFile)

		#Save ReadCounts in pickle
		with open(writePath+ '/SparseD.pkl', 'wb') as f:
			pickle.dump(SparseD, f, pickle.HIGHEST_PROTOCOL)
	
	
	elif mode == "deconvolution":
		with open(str(countpath)+ '/SparseD.pkl', 'rb') as f:
			SparseD = pickle.load(f)
		
		DBLmetricsDF = deconvolution(SparseD, vcf, GenotypesDF,barcodeList)
	
	elif mode == "refineDBLs":
		print("ComingSoon")
	
	else:
		
		SparseD = CountsMatrices()
		
		DBLmetricsDF = deconvolution(SparseD, vcf, GenotypesDF,barcodeList)
