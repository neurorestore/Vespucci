# README

Vespucci is an R package to prioritize spatial regions involved in the response to an experimental perturbation in spatial transcriptomics. 

The initution behind Vespucci is from a previously validated concept (https://github.com/neurorestore/Augur); namely, cell types undergoing the strongest transcriptional response to a given experimental perturbation are those for which perturbed cells can readily be separated from unperturbed cells by a machine-learning classifier, and that the accuracy (defined as Area-under-the-curve score; AUC) with which of this separation can be used as a quantitative measure of perturbation intensity. 

Vespucci extends this conceptual framework to spatial transcriptomics and identifies spatial regions of a comparative spatial transcriptomics dataset that are undergoing a transcriptional response to a perturbation. Briefly, Vespucci assigns each spatial barcode with an AUC score, measuring the transcriptomic separability of the barcode against *k*-nearest neighbours of contrasting experimental condition. This generates a spatial map of AUC of the coordiante system of spatial transcriptomics, thus allowing the identification of regions that are most perturbed by the experiment.

To mitigate the high computational cost of running a series of classification trials on every spatial barcode, Vespucci utilizes a meta-learning method where it estimates AUC scores using distance metrics between barcodes.  

Extending this, Vespucci also identifies differentially expressed (DE) genes or Gene Ontology (GO) terms that corresponds to the spatial AUC by fitting the feature values into a negative binomial generalized mixed model or generalized linear mixed model. 

## System requirements

Vespucci relies on functions from the following R packages:

```
	Augur (>= 1.0.3),
	dplyr (>= 0.8.0),
	purrr (>= 0.3.2),
	tibble (>= 2.1.3),
	magrittr (>= 1.5),
	tester (>= 0.1.7),
	Matrix (>= 1.2-14),
	sparseMatrixStats (>= 0.1.0),
	recipes (>= 0.1.4),
	rsample (>= 0.0.4),
	yardstick (>= 0.0.3),
	pbmcapply (>= 1.5.0),
	lmtest (>= 0.9-37),
	randomForest (>= 4.6-14),
	tidyselect (>= 0.2.5),
	rlang (>= 0.4.0),
	tidyr (>= 1.1.2),
	RANN (>= 2.6.1),
	dismay (>= 0.0.1),
	nebula (>= 1.4.1),
	Rdpack (>= 0.7),
	Seurat (>= 4.1.3)
```

Vespucci has been tested with R version 4.2.3.

## Installation

To install Vespucci, first install the devtools package, if it is not already installed: 

```r
> install.packages("devtools") 
```

Then, install Vespucci from GitHub: 

```r
> devtools::install_github("neurorestore/Vespucci")
```

This should take no more than a few minutes. 

## Usage

To derive the spatial AUC, run the `calculate_spatial_auc` function which takes as input a preprocessed features-by-cells (e.g., genes-by-cells) matrix, and a data frame containing metadata associated with each cell, minimally including the spatial coordinates, experimental labels, biological replicates and barcodes.
This means that in order to use Vespucci, you should have pre-processed your data (e.g., by read alignment and batch-effect removal) across all experimental labels. 
Furthermore, Vespucci requires the spatial transcriptomics dataset to be registered onto the a common coordinate framework. This can be done via various spatial registration pacakages, for example `RNiftyReg`. 

To run Vespucci with default parameters on a genes-by-cells scRNA-seq matrix `input`, and an accompanying data frame `meta`, with `x`, `y`, `label`, `replicate` and `barcode` columns (the x-coordinate, y-coordinate, the experimental label, biological replicate and barcode respectively), use the `run_vespucci` function:

```r
> vespucci_res = run_vespucci(input, meta)
```

If your columns have different names, you can specify these using the `coord_cols`, `label_col`, `replicate_col`, and `barcode_col` arguments:

```r
> vespucci_res = run_vespucci(
		input, 
		meta, 
		coord_cols = c('x-coord', 'y-coord'),
		label_col = 'condition', 
		replicate_col = 'new_replicate', 
		barcode_col = 'new_barcode'
	)
```

Vespucci can also be run on 3D coordinates by specifying the 3D coordinates.

```r
> vespucci_res = run_vespucci(input, meta, coord_cols = c('x', 'y', 'z'))
```

In the event that you have generated the distance metrics data frame prior to running \code{run_vespucci}, you can specify it as an argument to avoid re-calculating the distance metrics (which takes roughly 6-7 hours for 10,000 barcodes).

```r
> vespucci_res = run_vespucci(input, meta, distance_metrics_df = precalculated_distance_metrics_df)
```


## How to intepret results

Vespucci returns a list of two elements, `spatial_auc_result` and `de_feature_result`, outputs from the `calculate_spatial_auc` and `de_feature_result` functions respectively. 

```r
> spatial_auc_res = vespucci_res$spatial_auc_result
> de_feature_res = vespucci_res$de_feature_result
```

Spatial barcode prioritization values are stored in the `spatial_auc_res$aucs` data frame - for example:

```r
> head(spatial_auc_res$aucs)
 barcode       auc
1 barcode_1 0.5016050
2 barcode_2 0.5137284
3 barcode_3 0.5382223
4 barcode_4 0.4832105
5 barcode_5 0.8231267
6 barcode_6 0.9109776
  ...         ...
```

Differentially expressed (DE) features (genes/GO terms) that corresponds to these prioritization mapping can be obtained from `de_feature_result`

```r
> head(de_feature_res)
# A tibble: 6 × 4
  feature   effect_size       p_val p_val_adj
  <chr>           <dbl>       <dbl>     <dbl>
1 Gene5747         1.24 0.000000136   0.00248
2 Gene8610         1.25 0.00000206    0.0188 
3 Gene19683        1.13 0.00000997    0.0458 
4 Gene19490        1.24 0.0000101     0.0458 
5 Gene21281        1.77 0.0000196     0.0607 
6 Gene21662        1.39 0.0000200     0.0607
  ...         ...
```

By default, `run_vespucci` find DE genes. To find DE GO terms, parse a list of vectors of genes for Gene Ontology; each entry should be a vector of feature names as the argument `go_list`.

```r
> custom_go_list = list('GO_1' = c('Gene1', 'Gene2'), 'GO_2' = c('Gene3', 'Gene4'))
> vespucci_res = run_vespucci(input, meta, go_list = custom_go_list)
```


## Demonstration

To see Vespucci in action, load the simulated Seurat dataset that is bundled with the Vespucci package:

```r
> data("spatial_sim")
```

This dataset consists of 2000 cells with 2 experimental labels, 6 biological replicates, x-y coordinates (in the `label`, `replicate`, `x`, and `y` columns of the metadata). Furthermore, the metadata also has the column `true_label` which denotes if that specific barcode is `Perturbed` or `Control` that is manipulated by the simulation conditions. If you plot the x-y coordinates by the `true_label`, you should see a horizontal stripe in the middle showing the true perturbation pattern.

First, extract the expression and metadata.
```r
> input = GetAssayData(spatial_sim, slot='counts')
> meta = spatial_sim@meta.data
```

Now, we can run Vespucci with the following command,

```r
> vespucci_res = run_vespucci(input = input, meta = meta)
```
Note, this will take roughly 6 hours to run everything from scratch.

To speed up the process for this simulated dataset, we have generated the `distance_metrics_df` and `calculated_auc_df` from a run of this simulated dataset. To test this, run the following lines.

```r
> data('spatial_sim_distance_metrics_df')
> vespucci_res = run_vespucci(
    input = input,
    meta = meta,
    distance_metrics_df = spatial_sim_distance_metrics_df,
    test_sim = T # this only works for the spatial_sim dataset
	)
```

Finally, we can check the results.

```r
> spatial_auc_res = vespucci_res$spatial_auc_result
> de_feature_res = vespucci_res$de_feature_result
> head(spatial_auc_res$aucs)
    barcode       auc
1 barcode_1 0.5016050
2 barcode_2 0.5137284
3 barcode_3 0.5382223
4 barcode_4 0.4832105
5 barcode_5 0.8231267
6 barcode_6 0.9109776
> head(de_feature_res)
# A tibble: 6 × 4
  feature   effect_size       p_val p_val_adj
  <chr>           <dbl>       <dbl>     <dbl>
1 Gene5747         1.24 0.000000136   0.00248
2 Gene8610         1.25 0.00000206    0.0188 
3 Gene19683        1.13 0.00000997    0.0458 
4 Gene19490        1.24 0.0000101     0.0458 
5 Gene21281        1.77 0.0000196     0.0607 
6 Gene21662        1.39 0.0000200     0.0607 
```