
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






Counts = ReadCounter(GenotypeChunkIndexesList[0], bamFile, barcodeList, GenotypeChunkList, barcodetag="CB", umitag="UB")

CountsInt = Counts.copy()
CountsInt = CountsInt.replace({"A":"1","C":"2","T":"3","G":"4"}, regex=True).fillna(0)
ham = sklearn.metrics.pairwise_distances(CountsInt, Y=None, metric='hamming')




treadList = []
for l in [['cellranger_gex-GRCh38-2020-A_chrM', 8175],['cellranger_gex-GRCh38-2020-A_chrM', 13562]]:
	treadList = Pileupper_noRef(bam, l,set(barcodeList), treadList,barcodetag="CB", umitag="UB")



treadList = Pileupper_noRef(bam, GenotypeChunkList[0],set(barcodeList), [],barcodetag="CB", umitag="UB")


DFMaker_noRef(treadList, barcodeList, locus = GenotypeChunkList[0])

locusIndex = str(GenotypeChunkList[0][0])+"_"+str(GenotypeChunkList[0][1])
ReadsList=pd.DataFrame(treadList, columns = ["Pos","Barcode","Base"])[["Barcode","Base"]]
Counts = ReadsList.groupby(["Barcode"])["Base"].apply(lambda x: "".join(set(np.unique(x)))).to_frame(name=locusIndex)
Counts[locusIndex] =  np.where(Counts[locusIndex].str.len() == 2, Counts[locusIndex], Counts[locusIndex]*2)
return Counts







distance.hamming(Counts.loc["TTTGGAGCAGGACGAT-1"], Counts.loc["TTTGTTGAGATACAGT-1"])



i = 0
final = len(list(itertools.combinations(Counts.index.tolist(), 2)))
for cmb in itertools.combinations(Counts.index.tolist(), 2):
	t = distance.hamming(Counts.loc[cmb[0]], Counts.loc[cmb[0]])
	print(i/final)
	i = i+1

import sklearn
sklearn.metrics.pairwise_distances(Counts, Y=None, metric='hamming')

Counts.index.tolist()

