
TMPFILE=$(mktemp /hpcnfs/scratch/temporary/Dav_vc/20.4_DemultiMatteo/temp/$USER-XXXXXX-activate_conda.sh)
conda shell.zsh hook > $TMPFILE
set +u; source $TMPFILE; conda activate scansnp; set -u



import argparse
import os

parser = argparse.ArgumentParser(description='Process some integers.')


parser.add_argument('--vcf', dest='vcf', help='joint vcf with all mixed genotypes SNPs',required=True,type=str)


#conditionally optionals
parser.add_argument('--bam', dest='bam', help='for matrixgen mode: BAM file of mixed genotypes',type=str)
parser.add_argument('--barcodes', dest='barcodes', help='for matrixgen mode: FILTERED barcodes file in standard 10x .tsv format',type=str)
parser.add_argument('--pickle', dest='countpath',type=str,help='for deconvolution mode: path to folder containing Counts.npz countBarcodes.tsv and countLoci.tsv')

#Optionals
parser.add_argument('--threads', dest='nthreads', help='threads to be used',default=10,type=str)
parser.add_argument('--outdir', dest='outdir', help='use',default="currentwd",type=str)
parser.add_argument('--flagLowQual', dest='LowQual', help='Force-fit GMM to flag lowquality droplets',default=False,type=bool)
parser.add_argument('--mode', dest='mode', default="deconvolution",
	required=False,
	choices=["matrixgen","deconvolution","skipcount"],
	type=str,
	help='creating pileup matrix with ref and alt supporting reads')

args = parser.parse_args()

if args.mode == "matrixgen":
	if args.bam is None or args.barcodes is None:
		parser.error('please provide --bam and --barcodes for matrixgen mode')

if args.mode == "skipcount":
	if args.countpath is None:
		parser.error('please provide path to Counts.pkl using --pickle for skipcount mode')

if args.mode is None or args.mode == "deconvolution":
	if args.bam is None or args.barcodes is None:
		parser.error('please provide --bam and --barcodes')



mode="skipcount"
bamFile="/hpcnfs/scratch/temporary/Dav_vc/0_CellrangerAlign/Sample_S20272_157/outs/possorted_genome_bam.bam"
vcf="/hpcnfs/scratch/temporary/Dav_vc/2_Combining/Sample_S20273_157/output/combinedVCF/jointGenotype.g.vcf"
barcodesFILE="/hpcnfs/scratch/temporary/Dav_vc/0_CellrangerAlign/Sample_S20272_157/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
nThreads=10
outdir="/hpcnfs/scratch/temporary/Dav_vc/20.4_DemultiMatteo/temp"
countpath="/hpcnfs/scratch/temporary/Dav_vc/20.4_DemultiMatteo/scRNAseq/demultiplexing/organoidMultiplexing-ensembl_human_93/davide/Sample_S20272_157/scansnp/Counts.pkl"
LowQual=True
rawPath="/hpcnfs/scratch/temporary/Dav_vc/0_CellrangerAlign/Sample_S20272_157/outs/raw_feature_bc_matrix/"
filteredPath="/hpcnfs/scratch/temporary/Dav_vc/0_CellrangerAlign/Sample_S20272_157/outs/filtered_feature_bc_matrix/"


if outdir == "currentwd":
	outdir = os.getcwd()
else:
	outdir=outdir


if str(outdir).endswith(' '):
	outdir=str(outdir)[:-1]

if str(outdir).endswith('/'):
	outdir=str(outdir)[:-1]

if str(countpath).endswith('/'):
	countpath=str(countpath)[:-1]




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

sys.path.append('/hpcnfs/scratch/temporary/Dav_vc/20.4_DemultiMatteo/SCanSNP_V1_ssh_EmptyFit/SCanSNP/SCanSNP/')

from VCFUtils import *
from Wrappers import *
from RawBCMatrix_Utils import *



if __name__ == "__main__":
	#Creation of differen loci subsets
	cleanLoci = LociPreClean(vcf)
	MildcleanLoci = LociPreClean_milds(vcf)
	CleanSingularLoci,SingularLoci_Alt,SingularLoci_Ref = SingularLociSCan(vcf,cleanLoci)
	
	#If raw counts cellranger matrix is provided launch RawBCMatrix modules
	if rawPath not None:
		FullDropsKNNseries,EmptyBarcodeList = UnfilteredAdataAdata(filteredPath, rawPath, outdir, nHKgenes=50, raw_to_filt_ratio=1)
		barcodeList=pd.read_csv(barcodesFILE, header=None, names=["b"])["b"].tolist().extend(EmptyBarcodeList)
	else barcodeList = pd.read_csv(barcodesFILE, header=None, names=["b"])["b"].tolist()
	
	#Genotypes map creation
	GenotypesDF = creategenotypeDF(vcf)
	
	
	
	if mode == "matrixgen":
		'''
		Performing only count
		'''
		barcodeList=pd.read_csv(barcodesFILE, header=None, names=["b"])["b"].tolist()
		Counts = CountsMatrices(CleanSingularLoci, cleanLoci,
			MildcleanLoci, GenotypesDF,
			barcodeList, vcf,
			nThreads, bamFile)
		#Save ReadCounts in pickle
		with open(outdir+ '/Counts.pkl', 'wb') as f:
			pickle.dump(Counts, f, pickle.HIGHEST_PROTOCOL)
			
		
		
	elif mode == "skipcount":
		'''
		Performing only deconvolution based on provided counts
		'''
		#Load existing counts
		with open(str(countpath), 'rb') as f:
			Counts = pickle.load(f)
		#Deconvolution
		Cell_IDs = deconvolution(Counts, vcf, GenotypesDF,outdir, LowQual)
		Cell_IDs.to_csv(outdir + "/Cell_IDs.tsv", sep = "\t", header = True, index = True, index_label = "barcode")
		
		
	elif mode == "refineDBLs":
		print("ComingSoon")
		
		
	else:
		'''
		Perform both count and ceconvolution
		'''
		barcodeList=pd.read_csv(barcodesFILE, header=None, names=["b"])["b"].tolist()
		#counting
		Counts = CountsMatrices(CleanSingularLoci, cleanLoci,
			MildcleanLoci, GenotypesDF,
			barcodeList, vcf,
			nThreads, bamFile)
		#Save ReadCounts in pickle
		with open(outdir+ '/Counts.pkl', 'wb') as f:
			pickle.dump(Counts, f, pickle.HIGHEST_PROTOCOL)
		#Deconvolution
		Cell_IDs = deconvolution(Counts, vcf, GenotypesDF,outdir,LowQual)
		Cell_IDs.to_csv(outdir + "/Cell_IDs.tsv", sep = "\t", header = True, index = True, index_label = "barcode")





QualDF = mixtureModel(NoisedIDs, outdir, DBLsList, HighOutlierThreshold = .95,LowOutlierThreshold = .05 )
QualDF.Quality.value_counts()



gm = GaussianMixture(n_components=2, random_state=0, weights_init = [.05,.95],means_init= [[0],[4]]).fit(NoisedIDs["logFC"].to_numpy().reshape(-1,1))
QualLabels = gm.predict(NoisedIDs["logFC"].to_numpy().reshape(-1,1))
FittedMixturePlot(gm, outdir, NoisedIDs)
unique, counts = np.unique(QualLabels, return_counts=True)
dict(zip(unique, counts))

print(list(SNGsmetricsDF.FirstID).extend(list(SNGsmetricsDF.SecondID)))


list(set(list(SNGsmetricsDF.FirstID) + list(SNGsmetricsDF.SecondID)))

it




DBLscomb = list(itertools.combinations(list(set(list(SNGsmetricsDF.FirstID) + list(SNGsmetricsDF.SecondID))), 2))


786+9995+54+3587+21+355







#CountDF reconstruction
# SparseCounts=scipy.sparse.load_npz(str(countpath) + '/Counts.npz')
# Barcodes=pd.read_csv(str(countpath) + '/countBarcodes.tsv', header=None, names=["barcodes"])
# Loci=pd.read_csv(str(countpath) + '/countLoci.tsv', header=None, names=["loci"])

barcodeList = Counts["Barcode"]

#Pre-cleaning loci CHECK!!
#Pre-cleaning loci CHECK!!
cleanLoci = list(set(LociPreClean(vcf)).intersection(list(Counts["Locus"])))
MildcleanLoci = list(set(LociPreClean_milds(vcf)).intersection(list(Counts["Locus"])))

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
	DBLsDF=NoiseRregression(BestInDropDict, barcodeList, vcf, SingularLociScoreDF, BestBarcodeID)
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