# DenoIST: Denoising Image-based Spatial Transcriptomics data

<!-- badges: start -->

[![R-CMD-check](https://github.com/aaronkwc/DenoIST/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/aaronkwc/DenoIST/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

## Overview

DenoIST is a package for denoising image-based spatial transcriptomics data. It takes a IST count matrix and returns a adjusted count matrix with contamination removed.

The package is designed to be used with the [SpatialExperiment](https://bioconductor.org/packages/release/bioc/html/SpatialExperiment.html) class. If you are using a different format, a matrix input with a data frame of coordinates can also be accepted.

This is still very much a work in progress and we are still working on the documentation. Please feel free to open an issue or email at akwok@svi.edu.au if you have any questions or suggestions.

## News:
- **2026-05-24**: Minor update
  - Added `on_disk` mode to `denoist()` for large datasets that do not fit in memory.

- **2026-03-17**: Major update
  - Added `local_offset_distance_with_background_fast` function which uses `dbscan` for faster neighbour finding. Also fixed critical bug that causes wrong distance calculation.
  - `denoist()` now defaults to fast neighbour finding via the `neighbour_mode` option.

- **2026-02-24**: Minor update
  - Fixed minor bug where background offset cannot be calculated because an entire gene gets filtered out because of low qv. This should not change existing usage as the issue only arises in extremely small toy datasets.
  - `n_inits` can now be tuned in the `denoist()` function for speed.

- **2025-06-25**: Major update
  - Fixed memory usage for parallel processing. Feature is available only on linux/UNIX machines due to dependency on `parallel`.
  - Posterior cutoff is now a tunable parameter.
  - QOL changes including checking input types and error handling.

- **2025-05-20**: Initial release of DenoIST. The package is now available on and GitHub.

## Installation:

From Bioconductor:

```{r install}
if(!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("DenoIST")
```

Or from Github directly:

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

You can run `?denoist` for more details on the extra parameters you can adjust.

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
```

With a count matrix and coordinates:

```         
library(DenoIST)

res <- denoist(mat = mat,
               tx = tx,
               coords = coords,
               distance = 50, nbins = 200, cl = 1,
               out_dir = "denoist_results")
```

## Vignette

Check out the vignette to get started:

```
library(DenoIST)
browseVignettes("DenoIST")
```
