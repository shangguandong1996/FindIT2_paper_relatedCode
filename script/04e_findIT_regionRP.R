# Info ---------------------------
## Script name: findIT_regionRP

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
library(InteractiveFindIT2)
library(dplyr)

library(SummarizedExperiment)

# Set Options
options(stringsAsFactors = F)

# load up the data
Txdb <- AnnotationDbi::loadDb("rawdata/related_annotation/Txdb_gtf_Araport11.sqlite")
load("rawdata/ATAC/ATAC_peakNormCount.rda")
load("result/processed_result/wholeMergePeak_motif.rda")
mergePeak_GR_ATAC <- loadPeakFile("rawdata/ATAC/ATAC_mergePeak.bed")

load("result/processed_result/LEC2_targetTop1000.rda")

# regionRP ----------------------------------------------------------------

mmAnno <- mm_geneScan(peak_GR = mergePeak_GR_ATAC,
                      Txdb = Txdb,
                      upstream = 2e4,
                      downstream = 2e4)

# merge replicate
coldata <- data.frame(row.names = colnames(peakCountNorm_ATAC),
                      type = gsub("_R[0-9]", "", colnames(peakCountNorm_ATAC)))
peakCountNorm_ATAC_merge <- integrate_replicates(peakCountNorm_ATAC,
                                                 colData = coldata,
                                                 type = "value")

regionRP <- calcRP_region(mmAnno = mmAnno,
                          peakScoreMt = peakCountNorm_ATAC_merge,
                          Txdb = Txdb,
                          decay_dist = 1e3,
                          Chrs_included = paste0("Chr", 1:5))

set.seed(20160806)
regionRP_result <- findIT_regionRP(regionRP = regionRP,
                                   Txdb = Txdb,
                                   TF_GR_database = motif_GR,
                                   input_genes = LEC2_target)

TF_pvalue <- assays(regionRP_result)$TF_pvalue

merge_result <- c(regionRP, regionRP_result)
save(merge_result, file = "result/influential_TF/FindIT2_regionRP_merge_result.rda")
shinyParse_findIT_regionRP(merge_result)
shinyParse_findIT_regionRP(merge_result, mode = "TF")
