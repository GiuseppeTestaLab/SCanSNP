#!/usr/bin/env python

#VCF preprocessing utils


import itertools
import pandas as pd
import io
import re
from collections import Counter
import math


def ExtractSamples(vcf):
	'''
	Extraction of sample names from VCF
	'''
	with open(vcf, "r") as lociList:
		for line in (line for line in lociList if line.startswith("#CHR")):
			#samplePos=range(9,9+len(str(line).split()[9:]),1)
			samples=["_".join(["ID", samp]) for samp in  line.strip().split("\t")[9:]]
	return samples






def readVCF(vcf):
	'''
	Extraction of samples-level columns from VCF
	'''
	with open(vcf, 'r') as f:
		lines = [l for l in f if not l.startswith('#')]
	return pd.read_csv(io.StringIO(''.join(lines)), names = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FIELDS"]+[ sample for sample  in ExtractSamples(vcf)],sep='\t')[[sample for sample  in ExtractSamples(vcf) ]]


def ExtractIndex(vcf):
	'''
	Creates a list of pasted contigs and positions for every locus (e.g: MT_8080)
	'''
	lociIndex=[]
	with open(vcf, "r") as lociList:
		for locus in (locus for locus in lociList if not locus.startswith("#")):
			lociIndex.append("_".join([locus.split("\t")[0],locus.split("\t")[1]]))
	return(lociIndex)


def ExtractInfo(vcf):
	'''
	Extraction of informations about a specific VCF locus (i.e: POS,REF,ALT ...)
	'''
	chromosomes=[];positions=[];id=[];referenceAl=[];altAl=[]
	with open(vcf, "r") as lociList:
		for locus in (locus for locus in lociList if not locus.startswith("#")):
			chromosomes.append(locus.split("\t")[0])
			positions.append(locus.split("\t")[1])
			id.append(locus.split("\t")[2])
			referenceAl.append(locus.split("\t")[3])
			altAl.append(locus.split("\t")[4])
		InfoDict = {k:v for k,v in zip(["CHROM","POS","ID","REF","ALT"],[chromosomes,positions,id,referenceAl,altAl])}
	return InfoDict




def creategenotypeDF(vcf):
	'''
	Here we create a dataframe containing ref alleles and alt alleles count for every sample
	IGNORING DOUBLETS!!
	'''
	#Generating GenotypesDF from Index and Colnames
	colnames=["CHROM","POS","ID","REF","ALT"]+list(itertools.chain(*[ (sample+"_RefAl", sample + "_AltAl") for sample  in ExtractSamples(vcf)]))
	GenotypesDF2 = pd.DataFrame(index=ExtractIndex(vcf), columns=colnames).fillna(0)
	for col in ["CHROM","POS","ID","REF","ALT"]:
		GenotypesDF2[col] = ExtractInfo(vcf)[col]
	#Filling the missing GenotypesDf columns with genotype informations
	ReferenceVCF=readVCF(vcf)
	for column in ReferenceVCF:
		refgenotypes=ReferenceVCF[column].str.split(":", n = 1, expand = True)[0].str.count("0")
		altgenotypes=ReferenceVCF[column].str.split(":", n = 1, expand = True)[0].str.count("1")
		GenotypesDF2[str(column)+"_RefAl"]=list(refgenotypes)
		GenotypesDF2[str(column)+"_AltAl"]=list(altgenotypes)
	return GenotypesDF2



def LociPreClean(vcf):
	'''
	removes MULTIALLELIC LOCI and Loci with NON-POINT alterations
	removes missing cals in g.VCF
	'''
	CleanLociList =[]
	with open(vcf) as variantfile:
		for line in variantfile:
			if line.startswith('#'):
				continue
			else:
				ref=line.strip().split("\t")[3]
				alt=line.strip().split("\t")[4]
				sampleline=line.strip().split("\t")[9:]
				Genotypes=[i.split(':', 1)[0] for i in sampleline]
				GenotypesPhaseRemoved=[re.sub(r'\|', '/', geno) for geno in Genotypes]
				if ('./.' in GenotypesPhaseRemoved or "," in alt or len(ref) > 1 or len(alt) > 1 ):
					continue
				CleanLociList.append("_".join([line.split("\t")[0],line.split("\t")[1]]))
	return CleanLociList



def LociPreClean_milds(vcf):
	'''
	removes MULTIALLELIC LOCI and Loci with NON-POINT alterations
	allows 'holes' in jointvcf, useful for doublets detection
	'''
	MildCleanLociList =[]
	with open(vcf) as variantfile:
		for line in variantfile:
			if line.startswith('#'):
				continue
			else:
				ref=line.strip().split("\t")[3]
				alt=line.strip().split("\t")[4]
				sampleline=line.strip().split("\t")[9:]
				Genotypes=[i.split(':', 1)[0] for i in sampleline]
				GenotypesPhaseRemoved=[re.sub(r'\|', '/', geno) for geno in Genotypes]
				if ("," in alt or len(ref) > 1 or len(alt) > 1 ):
					continue
				MildCleanLociList.append("_".join([line.split("\t")[0],line.split("\t")[1]]))
	return MildCleanLociList




def SingularLociSCan(vcf,cleanLoci):
	'''
	Here we keep loci containing "Singular" information across the cohort: only 1 Ref Genotype (Homoz/Heteroz) or only 1 Alt genotype (Homoz/Heteroz)
	'''
	SingularLoci = []
	with open(vcf) as variantfile:
		SingularLoci_Alt = []
		SingularLoci_Ref = []
		for line in variantfile:
			if line.startswith('#'):
				continue
			else:
				alt=line.strip().split("\t")[4]
				sampleline=line.strip().split("\t")[9:]
				Genotypes=[i.split(':', 1)[0] for i in sampleline]
				nSamples = len(Genotypes)
				GenotypesPhaseRemoved=[re.sub(r'\|', '/', geno) for geno in Genotypes]
				AltGenotypes=["Alt" for genotype in GenotypesPhaseRemoved if "1" in genotype]
				RefGenotypes = ["Ref" for genotype in GenotypesPhaseRemoved if "0" in genotype]
				if (Counter(GenotypesPhaseRemoved)["0/1"] == 1 and len(RefGenotypes) == nSamples):
					SingularLoci_Alt.append("_".join([line.split("\t")[0],line.split("\t")[1]]))
				if (Counter(GenotypesPhaseRemoved)["0/1"] == 1 and len(AltGenotypes) == nSamples):
					SingularLoci_Ref.append("_".join([line.split("\t")[0],line.split("\t")[1]]))
				if (Counter(GenotypesPhaseRemoved)["1/1"] == 1 and Counter(GenotypesPhaseRemoved)["0/0"] == nSamples-1):
					SingularLoci_Alt.append("_".join([line.split("\t")[0],line.split("\t")[1]]))
				if (Counter(GenotypesPhaseRemoved)["0/0"] == 1 and Counter(GenotypesPhaseRemoved)["1/1"] == nSamples-1):
					SingularLoci_Ref.append("_".join([line.split("\t")[0],line.split("\t")[1]]))
	SingularLoci_Alt = list(set(SingularLoci_Alt).intersection(cleanLoci))
	SingularLoci_Ref = list(set(SingularLoci_Ref).intersection(cleanLoci))
	CleanSingularLoci = list(set(SingularLoci_Alt+SingularLoci_Ref))
	return CleanSingularLoci,SingularLoci_Alt,SingularLoci_Ref



def DiffOnlyIndexMaker(vcf,cleanLoci):
	'''
	removes loci containing the same allele across all the cohort of gVCF
	'''
	DiffOnlyIndex=[]
	with open(vcf) as variantfile:
		for line in variantfile:
			if line.startswith('#'):
				continue
			else:
				locusIndex = "_".join([line.split("\t")[0],line.split("\t")[1]])
				if locusIndex not in cleanLoci:
					continue
				alt=line.strip().split("\t")[4]
				sampleline=line.strip().split("\t")[9:]
				Genotypes=[i.split(':', 1)[0] for i in sampleline]
				GenotypesPhaseRemoved=[re.sub(r'\|', '/', geno) for geno in Genotypes]
				if (len(set(GenotypesPhaseRemoved)) == 1 ):
					continue
				DiffOnlyIndex.append("_".join([line.split("\t")[0],line.split("\t")[1]]))
	CleanDiffOnlyIndex = list(set(DiffOnlyIndex).intersection(cleanLoci))
	return CleanDiffOnlyIndex



def splitContig(contig, genotypes, nThreads):
	ChunkDICT = {}
	RemainsDF = pd.DataFrame()
	GenotypeChunk = genotypes[genotypes["CHROM"] == contig].sample(frac = 1)
	if len(GenotypeChunk) >= nThreads:
		chunkSizes = math.floor(len(GenotypeChunk)/int(nThreads))
		LeftOver = len(GenotypeChunk)%nThreads
		ChunkStart = 0
		for Chunk in range(0,nThreads):
			if Chunk < nThreads-1:
				ChunkDICT[Chunk] = GenotypeChunk.iloc[range(ChunkStart, ChunkStart+chunkSizes)]
			elif Chunk == nThreads-1:
				ChunkDICT[Chunk] = GenotypeChunk.iloc[range(ChunkStart, ChunkStart+chunkSizes+LeftOver)]
			ChunkStart = ChunkStart + chunkSizes
	else:
		RemainsDF = RemainsDF.append(GenotypeChunk)
	return ChunkDICT, RemainsDF




def ChunkMaker(GenotypesDF, nThreads, OmniIndex):
	GenotypeChunkDICT = {}
	Remains = pd.DataFrame()
	GenotypeChunkDICT_final = {}
	GenotypesDF_ss = GenotypesDF.loc[OmniIndex].sample(frac = 1)
	for contig in list(GenotypesDF["CHROM"].unique()):
		GenotypeChunkDICT[contig] = splitContig(contig, GenotypesDF_ss, nThreads)[0]
		Remains = Remains.append(splitContig(contig, GenotypesDF_ss, nThreads)[1])
	for Chunk in range(0,nThreads):
		GenotypeChunkDICT_final[Chunk] = pd.DataFrame()
		for key in GenotypeChunkDICT.keys():
			try:
				GenotypeChunkDICT_final[Chunk] = GenotypeChunkDICT_final[Chunk].append(GenotypeChunkDICT[key][Chunk])
			except:
				continue
	minChunk=min(GenotypeChunkDICT_final, key=lambda k: len(GenotypeChunkDICT_final[k]))
	GenotypeChunkDICT_final[minChunk] = GenotypeChunkDICT_final[minChunk].append(Remains)
	return GenotypeChunkDICT_final



def FlattenDict(GenotypeChunkDict):
	'''
	Transform dictionary into List of lists + ranges for efficiency
	'''
	GenotypeChunkList = []
	GenotypeChunkIndexesList = []
	ChunkStart = 0
	for k in GenotypeChunkDict.keys():
		chunk=GenotypeChunkDict[k].loc[:,["CHROM","POS","REF","ALT"]]
		# Create Chunk Index and append to list
		GenotypeChunkIndexesList.append(range(ChunkStart, ChunkStart+chunk.shape[0]))
		ChunkStart=ChunkStart+chunk.shape[0]
		# Re-format chunk to list-like object
		chunk["POS"] = chunk["POS"].astype(int)
		chunk=chunk.values.tolist()
		GenotypeChunkList.append(chunk)
	GenotypeChunkList = [item for sublist in GenotypeChunkList for item in sublist]
	return GenotypeChunkList, GenotypeChunkIndexesList


def ExtractMitoPositions(bamFile):
	mitoContig = ["MT","M"]
	
	bam=pysam.AlignmentFile(bamFile, "rb")
	
	
	
	