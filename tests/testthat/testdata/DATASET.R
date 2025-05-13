## code to prepare test dataset goes here
library(SubcellularSpatialData)
library(SpatialExperiment)
library(SingleCellExperiment)

data(tx_small)
spe <- tx2spe(tx_small, bin = "cell")
spe <- spe[,sample(1:ncol(spe), 1000)]
saveRDS(spe, file = "test_spe.rds")

#tx_test <- tx_small[sample(1:nrow(tx_small), 10000), ]

mat <- counts(spe)
# remove rows with rownames starting with BLANK or NegControl
mat <- mat[!grepl("BLANK|NegControl", rownames(mat)), ]
coords <- spatialCoords(spe)

saveRDS(mat, file = "test_mat.rds")
saveRDS(coords, file = "test_coords.rds")
saveRDS(tx_small, file = "test_tx.rds")
# further downsample to increase testing speed

#usethis::use_data(spe, overwrite = TRUE)
#usethis::use_data(mat, overwrite = TRUE)
#usethis::use_data(coords, overwrite = TRUE)
#usethis::use_data(tx_small, overwrite = TRUE)
