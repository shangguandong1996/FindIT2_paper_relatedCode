# Info ---------------------------
## Script name: findIT_MARA.R

## experiment_name:

## Purpose of script:
## 1. XX
## 2. XX

## Author: Guandong Shang
## Date Created: 2021-07-07
## Copyright (c) Guandong Shang, 2021
## Email: shangguandong1996@163.com

# Notes ---------------------------
## 1. XX
## 2. XX   
##

# Prepare -----------------------------------------------------------------

# load up the packages
library(FindIT2)
library(circlize)
library(ComplexHeatmap)

library(dplyr)

# Set Options
options(stringsAsFactors = F)

# load up the data
mergePeak_GR_ATAC <- loadPeakFile("rawdata/ATAC/ATAC_mergePeak.bed")
load("rawdata/ATAC/ATAC_peakNormCount.rda")
load("result/processed_result/wholeMergePeak_motif.rda")
load("rawdata/ATAC/ATAC_LRT_result.rda")

# MARA --------------------------------------------------------------------

motif_GR$feature_score <- motif_GR$score

LRT_result %>% 
    arrange(padj) %>% 
    slice(1:15000) %>% 
    pull(feature_id) -> input_feature_id

set.seed(20160806)
findIT_MARA(input_feature_id = input_feature_id,
            peak_GR = mergePeak_GR_ATAC,
            peakScoreMt = peakCountNorm_ATAC,
            TF_GR_database = motif_GR) -> result_MARA

readr::write_csv(result_MARA, file = "result/influential_TF/MARA.csv")

# intergrate replicate ----------------------------------------------------

colData <- data.frame(row.names = colnames(peakCountNorm_ATAC),
                      type = gsub("_R[0-9]", "",
                                  colnames(peakCountNorm_ATAC)))

MARA_mt <- as.matrix(result_MARA[, -1])
rownames(MARA_mt) <- result_MARA$TF_id
integrate_replicates(mt = MARA_mt,
                     colData = colData,
                     type = "rank_zscore") -> MARA_mt_merge

readr::write_csv(as_tibble(MARA_mt_merge, rownames = "TF_id"), 
                 file = "result/influential_TF/MARA_mt_merge.csv")



# Heatmap -----------------------------------------------------------------
apply(MARA_mt_merge, 1, sd) %>% 
    sort(decreasing = TRUE) %>% 
    names() %>% 
    "["(1:40) -> TF_name



col_fun_matrix <- colorRamp2(seq(-3, 3, length.out = 9), 
                             as.character(BuenColors::jdb_palette("solar_extra")))
Heatmap(MARA_mt_merge[TF_name, ],
        col = col_fun_matrix,
        
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        
        rect_gp = gpar(col = "white", lwd = 2),
        border = TRUE) -> hm

pdf("plot/MARA_top15000.pdf", height = 8)
draw(hm)
dev.off()
