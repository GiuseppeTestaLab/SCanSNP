
import pysam
from SCanSNP.VCFUtils import *
from scipy.sparse import csr_matrix
import sys

bamFile = "/group/testa/Common/scRNAseq/count/toxo-bis/human_and_toxomod-custom/matteo/ctl_run/multi/outs/per_sample_outs/CTL_rep1/count/sample_alignments.bam"
mitoContig = "cellranger_gex-GRCh38-2020-A_chrM"
barcodesFILE = "/group/testa/Common/scRNAseq/count/toxo-bis/human_and_toxomod-custom/matteo/ctl_run/multi/outs/per_sample_outs/CTL_rep1/count/sample_feature_bc_matrix/barcodes.tsv.gz"
nThreads = 4
barcodeList = pd.read_csv(barcodesFILE, header=None, names=["b"])["b"].astype("string").tolist()
locusT = ['cellranger_gex-GRCh38-2020-A_chrM', 13562]
BarcodeSet = set(barcodeList)


GenotypesDF = ExtractMitoPositions(bamFile, mitoContig)
MildcleanLoci = GenotypesDF.index.tolist()



GenotypeChunkList,GenotypeChunkIndexesList = FlattenDict(ChunkMaker(GenotypesDF, nThreads, MildcleanLoci))



def ReadCounter(chunkIdx, bamFile, barcodeList, GenotypeChunkList, barcodetag, umitag):
	bam=pysam.AlignmentFile(bamFile, "rb")
	BarcodeSet=set(barcodeList)
	readList = []
	Counts = {"sparse_Ref" : csr_matrix((0, len(barcodeList))), "sparse_Alt":csr_matrix((0, len(barcodeList))), "Locus" : pd.Series()  }
	lastPos = chunkIdx[-1]
	for indexPos in chunkIdx:
		
		locus = GenotypeChunkList[indexPos]
		
		readList=Pileupper(bam, locus, lastPos,BarcodeSet, readList,barcodetag, umitag)
		
		if sys.getsizeof(readList) >= 30000000 or indexPos == lastPos:
			try:
				sparse_Ref_TMP,sparse_Alt_TMP,locus_TMP = DFMaker(readList, barcodeList)
			except:
				continue
			Counts["sparse_Ref"] = scipy.sparse.vstack((Counts["sparse_Ref"],sparse_Ref_TMP))
			Counts["sparse_Alt"] = scipy.sparse.vstack((Counts["sparse_Alt"],sparse_Alt_TMP))
			Counts["Locus"] = Counts["Locus"].append(locus_TMP)
			readList = []
	bam.close()
	print("Chunk "+ str(chunkIdx) + " completed")
	return Counts


Counts = ReadCounter(GenotypeChunkIndexesList[0], bamFile, barcodeList, GenotypeChunkList, barcodetag="CB", umitag="UB")
ReadCounter(chunkIdx, bamFile, barcodeList, GenotypeChunkList, barcodetag, umitag)

treadList = Pileupper_noRef(bam, GenotypeChunkList[0],set(barcodeList), [],barcodetag="CB", umitag="UB")


DFMaker_noRef(treadList, barcodeList, locus = GenotypeChunkList[0])

locusIndex = str(GenotypeChunkList[0][0])+"_"+str(GenotypeChunkList[0][1])
ReadsList=pd.DataFrame(treadList, columns = ["Pos","Barcode","Base"])[["Barcode","Base"]]
Counts = ReadsList.groupby(["Barcode"])["Base"].apply(lambda x: "".join(set(np.unique(x)))).to_frame(name=locusIndex)
Counts[locusIndex] =  np.where(Counts[locusIndex].str.len() == 2, Counts[locusIndex], Counts[locusIndex]*2)
return Counts