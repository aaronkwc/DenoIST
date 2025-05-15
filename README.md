# DenoIST: Denoising Image-based Spatial Transcriptomics data

<!-- badges: start -->

[![R-CMD-check](https://github.com/aaronkwc/DenoIST/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/aaronkwc/DenoIST/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

## Overview

DenoIST is a package for denoising image-based spatial transcriptomics data. It takes a IST count matrix and returns a adjusted count matrix with contamination removed.

The package is designed to be used with the [SpatialExperiment](https://bioconductor.org/packages/release/bioc/html/SpatialExperiment.html) class. If you are using a different format, a matrix input with a data frame of coordinates can also be accepted.

## Installation:

```         
BiocManager::install(c('sparseMatrixStats', 'SpatialExperiment','SummarizedExperiment'))
devtools::install_github("aaronkwc/DenoIST")
```

## Quick start:

In most cases, you will only need to use the `denoist()` wrapper function.

It takes 2-3 inputs:

1.  `mat` : SpatialExperiment object (with the counts in assay() slot) or a count matrix with genes as rows and cells as columns.
2.  `tx`: Transcript data frame (a data frame with each row being an individual transcript, with columns specifying each transcripts' coordinates and qv).
3.  `coords`: If using a count matrix, a data frame (cells x 2) for each cell's centroid 2D coordinate.

The function will return a list with

1.  `adjusted_counts`: The adjusted counts matrix with contamination removed.
2.  `memberships`: A data frame with the inferred identity of each gene in each cell (real or contamination).
3.  `params`: A list with the estimated parameters used in the model.

Some additional parameters you can adjust based on your dataset:

1.  `distance` : The distance (in microns) to use for calculating local neighbourhood effects. Default is 50.
2.  `nbins` : Number of bins to use for calculating ambient background. Default is 200.
3.  `cl` : Number of cores to use for parallel processing. Default is 1.
4.  `out_dir` : An output directory to save the results in. Not mandatory. Default is NULL.

## Examples

With a SpatialExperiment object:

```         
library(DenoIST)
library(SpatialExperiment)

res <- denoist(mat = spe,
              tx = tx,
              coords = NULL,
              distance = 50, nbins = 200, cl = 1,
              out_dir = "denoist_results")
print(res$adjusted_counts)
```

With a count matrix and coordinates:

```         
library(DenoIST)

res <- denoist(mat = mat,
               tx = tx,
               coords = coords,
               distance = 50, nbins = 200, cl = 1,
               out_dir = "denoist_results")
print(res$adjusted_counts)
```
