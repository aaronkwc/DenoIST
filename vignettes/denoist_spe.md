---
title: "Denoising Imaged-based Spatial Transcriptomics data with DenoIST"
output: rmarkdown::html_vignette
author: "Aaron Kwok"
date: "2025-05-16"
vignette: >
  %\VignetteIndexEntry{denoist_spe}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



# Introduction

DenoIST (Denoising Image-based Spatial Transcriptomics) is a method for identifying and removing contamination artefacts in image-based single-cell transcriptomics (IST) data. This vignette shows how to use it with a `SpatialExperiment` object or a matrix with coordinates as a separate input.

# Load data

For demonstration, we will use a small Xenium sample from a lung fibrosis study (Vannan & Lyu et. al, 2025). It can be downloaded at [GSE250346](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250346).

From the raw Xenium output, we can then construct a `SpatialExperiment` object with `SpatialExperimentIO`.


``` r
suppressPackageStartupMessages({
  library(DenoIST)
  library(SpatialExperiment)
  library(ggplot2)
  library(patchwork)
})
```


``` r
dir = "/mnt/beegfs/mccarthy/backed_up/general/rlyu/Dataset/LFST_2022/GEO_2025/VUHD116A/relabel_output-XETG00048__0003817__VUHD116A__20230308__003730/outs/"
spe <- readXeniumSXE(dir, returnType = "SPE")
saveRDS(spe, "example_spe.rds")
```

For this vignette, we will use a pre-saved `SpatialExperiment` object generated from the code above.


``` r
spe <- readRDS('example_spe.rds')
spe
#> class: SpatialExperiment 
#> dim: 343 12915 
#> metadata(4): experiment.xenium transcripts cell_boundaries nucleus_boundaries
#> assays(1): counts
#> rownames(343): ABCC2 ACKR1 ... YAP1 ZEB1
#> rowData names(3): ID Symbol Type
#> colnames(12915): aaaaaaab-1 aaaaaaac-1 ... aaaadchc-1 aaaadchd-1
#> colData names(10): cell_id transcript_counts ... nucleus_area sample_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(3): NegControlProbe UnassignedCodeword NegControlCodeword
#> spatialCoords names(2) : x_centroid y_centroid
#> imgData names(0):
```


``` r
tx <- read.csv('VUHD116A_transcripts.csv')
head(tx)
#>   X            transcript_id             cell_id overlaps_nucleus            feature_name
#> 1 1 VUHD116A_281474976711003 VUHD116A_UNASSIGNED                0 NegControlCodeword_0517
#> 2 2 VUHD116A_281474976711398 VUHD116A_UNASSIGNED                0                  COL1A1
#> 3 3 VUHD116A_281474976711406 VUHD116A_UNASSIGNED                0                     LYZ
#> 4 4 VUHD116A_281474976711409 VUHD116A_UNASSIGNED                0                   LAMP3
#> 5 5 VUHD116A_281474976711415 VUHD116A_UNASSIGNED                0                     LYZ
#> 6 6 VUHD116A_281474976711418 VUHD116A_UNASSIGNED                0                     LYZ
#>   x_location y_location z_location        qv fov_name nucleus_distance
#> 1 202.896320   184.3326   17.63912  2.910288      Q15         329.0720
#> 2   7.082226   110.4102   12.68721 40.000000      Q15         537.2145
#> 3  76.882706   196.5678   12.74130 40.000000      Q15         439.8191
#> 4 101.830030   107.3766   13.04323 40.000000      Q15         455.0005
#> 5 149.067670   231.7206   12.96732 10.863165      Q15         360.2229
#> 6 175.447590   171.5916   13.25221  3.295905      Q15         359.2819
```


# Denoising the data

You should only need to use 1 function most of the time, unless you are trying to debug or understand the process. The main function is `denoist()`, which takes a `SpatialExperiment` object (or a matrix with coordinates), plus the transcript data frame as input. It will return a list containing the memberships, adjusted counts, and parameters for each gene.

The `distance` parameter specifies the maximum distance to consider for local background estimation. The `nbins` parameter specifies the number of bins to use for hexagonal binning, which is used for calculating background transcript contamination. They have default values but you can adjust them based on your data. For example if your data is very small in size then perhaps a lower `distance` and `nbins` would be better.

You should also check whether transcript data frame has the correct columns. The `tx_x` and `tx_y` parameters specify the column names for the x and y coordinates, respectively. The `feature_label` parameter specifies the column name for the gene of each transcript. In this example, they are called `x_location`, `y_location` and `feature_name`. You can also speed up the process with more cpus via the `cl` option (which is highly recommended).

Lastly, you can specify an output directory with `out_dir` to save the results automatically. If you don't want to, just leave it empty.















