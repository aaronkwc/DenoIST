# DenoIST: Denoising Image-based Spatial Transcriptomics data

<!-- badges: start -->

[![R-CMD-check](https://github.com/aaronkwc/DenoIST/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/aaronkwc/DenoIST/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

## Overview

DenoIST is a package for denoising image-based spatial transcriptomics data. It takes a IST count matrix and returns a adjusted count matrix with contamination removed.

The package is designed to be used with the [SpatialExperiment](https://bioconductor.org/packages/release/bioc/html/SpatialExperiment.html) class. If you are using a different format, a matrix input with a data frame of coordinates can also be accepted.

This is still very much a work in progress and we are still working on the documentation. Please feel free to open an issue or email at akwok@svi.edu.au if you have any questions or suggestions.

## News:

- **2025-05-20**: Initial release of DenoIST. The package is now available on and GitHub.

## Installation:

```         
BiocManager::install(c('sparseMatrixStats', 'SpatialExperiment','SummarizedExperiment'))
devtools::install_github("aaronkwc/DenoIST")
```

## Quick start:

In most cases, you will only need to use the `denoist()` wrapper function.

It takes 2-3 inputs:

1.  `mat` : SpatialExperiment object (with the counts in assay() slot) or a count matrix with genes as rows and cells as columns.
2.  `tx`: Transcript data frame (a data frame with each row being an individual transcript, with columns specifying each transcripts' coordinates and qv). If your transcript file is not from Xenium and has no qv score, you can set a dummy column of `qv = 20` for all transcripts. This workaround should not be needed in future updates.
3.  `coords`: If using a count matrix, a data frame (cells x 2) for each cell's centroid 2D coordinate.

The function will return a list with

1.  `adjusted_counts`: The adjusted counts matrix with contamination removed.
2.  `memberships`: A data frame with the inferred identity of each gene in each cell (1 for real or 0 for contamination).
3.  `params`: A list with the estimated parameters used in the model. The posterior probabilities of each gene being real or contamination can be found in `params$posterior_probs`, higher means more likely to be contamination.

Some additional parameters you can adjust based on your dataset:

1.  `distance` : The distance (in microns) to use for calculating local neighbourhood effects. Default is 50.
2.  `nbins` : Number of bins to use for calculating ambient background. Default is 200.
3.  `tx_x` and `tx_y`: The column names in the transcript data frame for the x and y coordinates of each transcript. Default is `x` and `y`.
4.  `feature_label`: The column name in the transcript data frame for the gene of each transcript. Default is `gene`. (In Xenium you should change it to `feature_label`)
5.  `cl` : Number of cores to use for parallel processing. Default is 1.
6.  `out_dir` : An output directory to save the results in. Not mandatory. Default is NULL.

You can run `?denoist` for more details on the parameters.

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

## Vignette

A brief vignette:

[Denoising healthy lung Xenium data](https://rawcdn.githack.com/aaronkwc/DenoIST/e6683f326a34bc5f779077e9a0435ec8ec2ce831/vignettes/denoist_spe.html)
