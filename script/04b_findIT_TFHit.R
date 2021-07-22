# Info ---------------------------
## Script name: findIT_TFHit.R

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
library(dplyr)

# Set Options
options(stringsAsFactors = F)

# load up the data
Txdb <- AnnotationDbi::loadDb("rawdata/related_annotation/Txdb_gtf_Araport11.sqlite")

public_data <- loadPeakFile("rawdata/related_annotation/remap2020_all_macs2_TAIR10_v1_0.bed.gz")
colnames(mcols(public_data))[1] <- "TF_id"
mcols(public_data) <- mcols(public_data)[1]
seqlevels(public_data) <- paste0("Chr", 1:5)

load("rawdata/ChIP/TF_ChIP_bindingSite_GR.rda")
LEC2_bindingSite_GR$TF_id <- "LEC2"
mcols(LEC2_bindingSite_GR) <- mcols(LEC2_bindingSite_GR)[3]

merge_database <- c(public_data, LEC2_bindingSite_GR)

load("result/processed_result/LEC2_targetTop1000.rda")

load("result/processed_result/wholeMergePeak_motif.rda")

# TFHit -------------------------------------------------------------------

set.seed(19960203)
findIT_TFHit(input_genes = LEC2_target,
             TF_GR_database = merge_database,
             Txdb = Txdb,
             scan_dist = 2e4,
             decay_dist = 1e3) -> TFHit_public

set.seed(19960203)
findIT_TFHit(input_genes = LEC2_target,
             TF_GR_database = motif_GR,
             Txdb = Txdb,
             scan_dist = 2e4,
             decay_dist = 1e3) -> TFHit_motif

readr::write_csv(TFHit_motif, file = "result/influential_TF/TFHit_motif.csv")
readr::write_csv(TFHit_public, file = "result/influential_TF/TFHit_public.csv")
