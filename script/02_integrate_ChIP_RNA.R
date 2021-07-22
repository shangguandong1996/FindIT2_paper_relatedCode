# Info ---------------------------
## Script name: integrate_ChIP_RNA.R

## experiment_name:

## Purpose of script:
## 1. XX
## 2. XX

## Author: Guandong Shang
## Date Created: 2021-07-06
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

load("rawdata/RNA/LEC2_RNA_diffResult.rda")
load("rawdata/ChIP/TF_ChIP_bindingSite_GR.rda")

gene_alias <- readr::read_csv("rawdata/related_annotation/gene_aliases_20190630_tidy.csv")

dir.create("result/integrate_ChIP_RNA")


# calcRP_TFHit ------------------------------------------------------------

# First, we check TF's category, promoter or enhancer
mmAnno_nearest <- mm_nearestGene(LEC2_bindingSite_GR,
                                 Txdb = Txdb)
plot_annoDistance(mmAnno_nearest)

# scan 10000 bp to get related peak
mmAnno <- mm_geneScan(LEC2_bindingSite_GR,
                      Txdb = Txdb,
                      upstream = 10000,
                      downstream = 10000)

calcRP_TFHit(mmAnno = mmAnno,
             Txdb = Txdb,
             decay_dist = 1000) -> RP_TFHit

integrate_ChIP_RNA(result_geneRP = RP_TFHit,
                   result_geneDiff = diff_RNA) -> mergeResult

mergeResult %>% 
    left_join(gene_alias, by = c("gene_id" = "name")) -> mergeResult

readr::write_csv(mergeResult, file = "result/integrate_ChIP_RNA/mergeResult.csv")
#readr::write_delim(mergeResult, file = "~/mergeResult.txt")

# output LEC2 top 1000 target gene
LEC2_target <- mergeResult %>% 
  dplyr::slice(1:1000) %>% 
  pull(gene_id)  
save(LEC2_target, file = "result/processed_result/LEC2_targetTop1000.rda")
