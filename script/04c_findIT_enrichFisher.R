# Info ---------------------------
## Script name: findIT_enrichFisher.R

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
mergePeak_GR_ATAC <- loadPeakFile("rawdata/ATAC/ATAC_mergePeak.bed")
load("rawdata/ATAC/ATAC_peakNormCount.rda")
load("result/processed_result/wholeMergePeak_motif.rda")
load("result/processed_result/mmAnno_enhancerPromoter.rda")

Txdb <- AnnotationDbi::loadDb("rawdata/related_annotation/Txdb_gtf_Araport11.sqlite")

public_data <- loadPeakFile("rawdata/related_annotation/remap2020_all_macs2_TAIR10_v1_0.bed.gz")
colnames(mcols(public_data))[1] <- "TF_id"
mcols(public_data) <- mcols(public_data)[1]
seqlevels(public_data) <- paste0("Chr", 1:5)

load("rawdata/ChIP/TF_ChIP_bindingSite_GR.rda")
LEC2_bindingSite_GR$TF_id <- "LEC2"
mcols(LEC2_bindingSite_GR) <- mcols(LEC2_bindingSite_GR)[3]

merge_database <- c(public_data, LEC2_bindingSite_GR)

motif_GR$TF_score <- motif_GR$score

load("result/processed_result/LEC2_targetTop1000.rda")

# get related peak --------------------------------------------------------

enhancer_related <- subset(mmAnno_enhancerPromoter, gene_id %in% LEC2_target & 
                               cor > 0.8 & 
                               pvalue < 0.01)$feature_id
promoter_related <- metadata(mmAnno_enhancerPromoter)$mm_promoter_tidy %>% 
    filter(gene_id %in% LEC2_target) %>% 
    pull(promoter_feature) 

gene_BoundRelated <- mm_geneBound(peak_GR = mergePeak_GR_ATAC,
                                  Txdb = Txdb,
                                  input_genes = LEC2_target)$feature_id

related_feature_id <- unique(c(enhancer_related, promoter_related, 
                               gene_BoundRelated))

save(related_feature_id, 
     file = "result/processed_result/LEC2_targetTop1000_relatedATACPeak.rda")

# enrichFisher -------------------------------------------------------------

findIT_enrichFisher(input_feature_id = related_feature_id,
                    peak_GR = mergePeak_GR_ATAC,
                    TF_GR_database = motif_GR) -> enrichFisher_motif                               

findIT_enrichFisher(input_feature_id = related_feature_id,
                    peak_GR = mergePeak_GR_ATAC,
                    TF_GR_database = merge_database) -> enrichFisher_public

readr::write_csv(enrichFisher_motif, file = "result/influential_TF/enrichFisher_motif.csv")
readr::write_csv(enrichFisher_public, file = "result/influential_TF/enrichFisher_public.csv")


