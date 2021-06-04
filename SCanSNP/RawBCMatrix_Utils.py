
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
import anndata as ad
from scipy.sparse import csr_matrix
import numpy as np

def UnfilteredAdataAdata(filteredPath, rawPath, outdir, nHKgenes=50, raw_to_filt_ratio=1, n_neighbors=20):
	"""
	RawFilt ratio is the number of raw cells per filtered to retain after the thresholded value selected from cellranger
	KNN will scale bad, be wise...
	"""
	
	UnfilteredAdata=sc.read_10x_mtx(rawPath,var_names='gene_symbols', cache=False)
	FilteredAdata=sc.read_10x_mtx(filteredPath,var_names='gene_symbols', cache=False)
	
	#Store filtered barcodes
	FilteredBCs = FilteredAdata.obs_names.tolist()
	
	sc.pp.calculate_qc_metrics(UnfilteredAdata,inplace=True)
	
	#Get LowDropoutRate genes from unfilt
	LDGs=UnfilteredAdata.var.sort_values(by=["pct_dropout_by_counts"], ascending=False).tail(nHKgenes).index.tolist()
	
	#Select closest raw_to_filt_ratio*#filtered_cells from thresholded cellranger value
	UnfilteredAdata = UnfilteredAdata[~UnfilteredAdata.obs_names.isin( FilteredAdata.obs_names)]
	sc.pp.calculate_qc_metrics(UnfilteredAdata,inplace=True)
	UnfilteredAdata = UnfilteredAdata[UnfilteredAdata.obs.sort_values(by=["total_counts"], ascending=False).head(len(FilteredBCs)*raw_to_filt_ratio).index.tolist(), LDGs]
	
	#Subset also FilteredAdata for LDGss
	sc.pp.calculate_qc_metrics(FilteredAdata,inplace=True)
	FilteredAdata=FilteredAdata[:,LDGs]
	
	#Concat Both adata
	adata = ad.concat([FilteredAdata,UnfilteredAdata], merge="same",label="DropletType", keys=["Full","Empty"])
	del UnfilteredAdata 
	del FilteredAdata
	
	sc.set_figure_params(figsize =[15,10])
	sc.settings.autosave = True
	sc.settings.figdir = outdir
	
	sc.pp.log1p(adata)
	sc.pp.scale(adata)
	
	sc.tl.pca(adata, svd_solver='arpack')
	
	
	#KNN computation
	print("Computing KNNs. may lead to increase memory usage")
	nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm='ball_tree').fit(adata.obsm["X_pca"][:,:2])
	KNN=nbrs.kneighbors_graph(adata.obsm["X_pca"][:,:2]).tocsr()
	
	FullDrops=csr_matrix((np.array(adata.obs["DropletType"]) == "Full" ).astype(int))
	EmptyDrops=csr_matrix((np.array(adata.obs["DropletType"]) == "Empty" ).astype(int))
	
	#Count EmptyDrops in KNN
	EmptyDrops=EmptyDrops.multiply(KNN)
	EmptyDrops=EmptyDrops.sum(axis = 1)
	
	#Count FullDrops in 15NN
	FullDrops=FullDrops.multiply(KNN)
	FullDrops=FullDrops.sum(axis = 1)
	
	#Compute ratio per drop
	FullDropsKNN=np.divide(FullDrops,(EmptyDrops+FullDrops))
	adata.obs["FullDropsKNN"] = FullDropsKNN
	sc.pl.embedding(adata, basis = "pca", color=['DropletType','total_counts','FullDropsKNN'], size = 100, alpha =.3, add_outline = True,  outline_width=(0.3, 0.01), ncols = 1, palette = "viridis", save = "DropletType.svg")
	
	#Prepare the function return
	FullDropsKNNseries = adata.obs.loc[FilteredBCs,"FullDropsKNN" ].to_frame()
	EmptyBarcodeList = adata.obs_names[adata.obs["DropletType"] == "Empty"].tolist()
	
	return FullDropsKNNseries, EmptyBarcodeList



