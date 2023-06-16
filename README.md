# SCanSNP
Version prepared for consensus

## __Dependencies manually installed in simulated home__

conda create -n SCanSNP python==3.8

conda activate SCanSNP

conda install -c conda-forge r-base=3.5.1 r-ggplot2 r-mixtools r-dplyr r-robustbase r-cluster r-gridExtra 

pip3 install git/SCanSNP

#conda install -c anaconda readline=6.2
#libreadline.so.6

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
