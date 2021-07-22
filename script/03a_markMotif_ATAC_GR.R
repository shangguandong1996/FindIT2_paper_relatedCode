# Info ---------------------------
## Script name: markMotif_ATAC_GR.R

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
library(motifmatchr)
library(BSgenome.Athaliana.TAIR.TAIR9)

# Set Options
options(stringsAsFactors = F)

# load up the data
load("rawdata/related_annotation/motifs_JASPAR2020.rda")
mergePeak_GR_ATAC <- FindIT2::loadPeakFile("rawdata/ATAC/ATAC_mergePeak.bed")

# motif scan --------------------------------------------------------------

# ATAC_GR_resize <- resize(mergePeak_GR_ATAC, width = 500, fix = "center")
# seqinfo(ATAC_GR_resize) <- seqinfo(BSgenome.Athaliana.TAIR.TAIR9)
# ATAC_GR_resize <- trim(ATAC_GR_resize)
motif_ix <- matchMotifs(motifs_JASPAR2020, 
                        mergePeak_GR_ATAC, 
                        genome = BSgenome.Athaliana.TAIR.TAIR9, 
                        out = "positions", 
                        bg = "genome")

motif_GR <- unlist(motif_ix)
motif_GR$TF_id <- gsub("(MA[0-9]{4}\\.[0-9]_)|(AT[0-9]G[0-9]{5}_)", "", names(motif_GR))

save(motif_GR, file = "result/processed_result/wholeMergePeak_motif.rda")
