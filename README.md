# SCanSNP

## __Installation__

<pre>
conda create -n SCanSNP python==3.8
conda activate SCanSNP
conda install -c conda-forge r-base=3.5.1 r-ggplot2 r-mixtools r-dplyr r-robustbase r-cluster r-gridExtra 
git clone git@github.com:GiuseppeTestaLab/SCanSNP.git
pip3 install ./SCanSNP
</pre>

## __Or docker image__

<pre>
docker pull davidecast/scansnp:1.1.0
</pre>

## __Usage__

<pre>
SCanSNP --vcf $multi_sample.vcf --barcodes $barcodes.tsv --bam $bam_file --outdir $output_directory
</pre>

### __additional args__
<pre>
--mode {deconvolution,matrixgen,skipcount}
	anndata with ref and alt reads information [not extensively tested]
 
	deconvolution: default mode will run through the full demultiplexing given 
	input vcf, barcodes list and bamfile
       
	matrixgen: scansnp will only generate anndata with nBarcodes x nLoci 
	counts for each allel will be stored in adata.layers["RefReads"] and adata.layers["AltReads"] 
	the default .X will be same as adata.layers["RefReads"]
	 
	skipcount: pick up the demultiplexing starting from previously saved anndata
	
--counts COUNTPATH: path to previously obtained anndata (--mode matrixgen) mandatory if --mode skipcount
--threads NTHREADS    threads to be used
--segmentation Path to segmentation tsv file with barcodes X number of nuclei. 
	When --platform visium is specified this file can be provided to improve the signal to noise calculation and will unlock the formal assignment of multiple genotypes per pack

--platform {chromium,visium}
</pre>


### 10.1.2022 - Changes

**Modifications in `SCanSNP.py`**

- fixed bug:

	```{python3}
	""" Traceback (most recent call last):
	File "/home/davide.castaldi/git/SCanSNP/SCanSNP/SCanSNP.py", line 177, in
	Counts = CountData(varAdata.layers["sparse_Ref"], varAdata.layers["sparse_Alt"], varAdata.var_names, varAdata.obs_names)
	File "/home/davide.castaldi/.local/lib/python3.6/site-packages/anndata/_core/aligned_mapping.py", line 148, in getitem
	return self._data[key]
	KeyError: `sparse_Ref` """
	```

	by changing `sparse_Ref` and `sparse_Alt` to `RefReads` and `AltReads`

**Modification in `GenUtils.py`**

- in `CountData` class: when initializing the class, `sparseAlt` and `sparseRef` matrices are transposed. This allow the program to load and correctly use a previously saved anndata (varAdata.h5ad) generated with one of the other SCanSNP modalities where both reference and alternative matrices are saved in the object as barcodes x loci matrices.

**Modification in `classifyUtils.py`**

- in `barcodeClassifier_main()`:
	* Added the computation for each ID of the number of negative and positive barcode and of their percentage over the total. Then:
		* if a certain ID has a low number of positive barcodes (in percentage), then the computation of the positive rate and the model fitting for the classification is skipped and instead we add an empty string for all the barcodes. The chosen minimum threshold is 1% of positive barcodes over the total.

		* if a certain ID has a low number of negative barcodes (in percentage), then the computation of the positive rate and the model fitting for the classification is skipped and instead we add that ID to all barcodes. The chosen minimum threshold is 0.5% of positive barcodes over the total.    

__NB:__ The thresholds have not been tested.
