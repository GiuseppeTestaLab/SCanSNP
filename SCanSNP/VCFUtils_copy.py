#!/usr/bin/env python

#VCF preprocessing utils


import itertools
import pandas as pd
import io
import re
from collections import Counter
import math
import copy


#vcf = "/group/testa/Project/DemultiplexingSpatial/jointGenotype_modified.g.vcf"

class VCFdata:
	
	def __init__(self, vcfpath):
		self.path = vcfpath
		self.vcfFrame = self.readVCF()
	
	
	def copy(self):
		return copy.copy(self)
	
	def ExtractSamples(self):
		with open(self.path, "r") as lociList:
			for line in (line for line in lociList if line.startswith("#CHR")):
				#samplePos=range(9,9+len(str(line).split()[9:]),1)
				samples=["_".join(["ID", samp]) for samp in  line.strip().split("\t")[9:]]
		return samples
	
	def readVCF(self):
		'''
		Read WHOLE vcf removing following loci
		- Loci with constant genotype (no info)
		- Multiallelic Loci
		'''
		with open(self.path) as variantfile:
			lines = []
			for line in variantfile:
				if line.startswith('#'):
					continue
				locusIndex = "_".join([line.split("\t")[0],line.split("\t")[1]])
				alt=line.strip().split("\t")[4]
				ref=line.strip().split("\t")[3]
				sampleline=line.strip().split("\t")[9:]
				Genotypes=[i.split(':', 1)[0] for i in sampleline]
				GenotypesPhaseRemoved=[re.sub(r'\|', '/', geno) for geno in Genotypes]
				if (len(set(GenotypesPhaseRemoved)) == 1 or "," in alt or len(ref) > 1 or len(alt) > 1):
					continue
				lines.append(line)
			
		VCF =  pd.read_csv(io.StringIO(''.join(lines)), names = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FIELDS"]+[ sample for sample  in self.ExtractSamples()],sep='\t')
		VCF.index = VCF["CHROM"].astype(str)+"_"+VCF["POS"].astype(str)
		return VCF
	
	def ExtractInfo(self):
		'''
		Extraction of informations about a specific VCF locus (i.e: POS,REF,ALT ...)
		'''
		chromosomes=[];positions=[];id=[];referenceAl=[];altAl=[]
		for index, locus in self.vcfFrame.iterrows():
			chromosomes.append(locus["CHROM"])
			positions.append(locus["POS"])
			id.append(locus["ID"])
			referenceAl.append(locus["REF"])
			altAl.append(locus["ALT"])
		InfoDict = {k:v for k,v in zip(["CHROM","POS","ID","REF","ALT"],[chromosomes,positions,id,referenceAl,altAl])}
		return InfoDict
	
	
	def creategenotypeDF(self):
		'''
		Here we create a dataframe containing ref alleles and alt alleles count for every sample
		'''
		vcfFrame = self.readVCF()
		samples = self.ExtractSamples()
		#Generating GenotypesDF from Index and Colnames
		colnames=["CHROM","POS","ID","REF","ALT"]+list(itertools.chain(*[ (sample+"_RefAl", sample + "_AltAl") for sample  in samples]))
		GenotypesDF = pd.DataFrame(index=vcfFrame.index.tolist(), columns=colnames).fillna(0)
		for col in ["CHROM","POS","ID","REF","ALT"]:
			GenotypesDF[col] = vcfFrame[col]
		#Filling the missing GenotypesDf columns with genotype informations
		ReferenceVCF=vcfFrame
		for sample in samples:
			refgenotypes=ReferenceVCF[sample].str.split(":", n = 1, expand = True)[0].str.count("0")
			altgenotypes=ReferenceVCF[sample].str.split(":", n = 1, expand = True)[0].str.count("1")
			GenotypesDF[str(sample)+"_RefAl"]=list(refgenotypes)
			GenotypesDF[str(sample)+"_AltAl"]=list(altgenotypes)
		return GenotypesDF






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
		RemainsDF = pd.concat([RemainsDF,GenotypeChunk])
	return ChunkDICT, RemainsDF




def ChunkMaker(GenotypesDF, nThreads):
	GenotypeChunkDICT = {}
	Remains = pd.DataFrame()
	GenotypeChunkDICT_final = {}
	for contig in list(GenotypesDF["CHROM"].unique()):
		GenotypeChunkDICT[contig] = splitContig(contig, GenotypesDF, nThreads)[0]
		Remains = pd.concat([Remains, splitContig(contig, GenotypesDF, nThreads)[1]])
	for Chunk in range(0,nThreads):
		GenotypeChunkDICT_final[Chunk] = pd.DataFrame()
		for key in GenotypeChunkDICT.keys():
			try:
				GenotypeChunkDICT_final[Chunk] = pd.concat([GenotypeChunkDICT_final[Chunk], GenotypeChunkDICT[key][Chunk]])
			except:
				continue
	minChunk=min(GenotypeChunkDICT_final, key=lambda k: len(GenotypeChunkDICT_final[k]))
	GenotypeChunkDICT_final[minChunk] = pd.concat([GenotypeChunkDICT_final[minChunk],Remains])
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