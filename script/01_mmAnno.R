# Info ---------------------------
## Script name: mmAnno.R

## experiment_name:
experiment_name <- "ATAC"

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
library(ggplot2)

library(dplyr)

# Set Options
options(stringsAsFactors = F)

# load up the data
mergePeak_GR <- loadPeakFile("rawdata/ATAC/ATAC_mergePeak.bed")

load("rawdata/ATAC/ATAC_peakNormCount.rda")

Txdb <- AnnotationDbi::loadDb("rawdata/relate_annotation/Txdb_gtf_Araport11.sqlite")

gene_alias <- readr::read_csv("rawdata/relate_annotation/gene_aliases_20190630_tidy.csv")

dir.create("plot/mmAnno")
dir.create("result/mmAnno")
dir.create("result/processed_result")

# mmAnno ------------------------------------------------------------------

mmAnno_nearest <- mm_nearestGene(peak_GR = mergePeak_GR,
                                 Txdb = Txdb)

# mmAnno_geneScan <- mm_geneScan(peak_GR = mergePeak_GR,
#                                Txdb = Txdb,
#                                upstream = 20000,
#                                downstream = 20000)

#subset(mmAnno_geneScan, gene_id == "AT2G05940")

# mmAnno summary ----------------------------------------------------------

getAssocPairNumber(mmAnno_nearest) %>% 
    arrange(desc(peakNumber)) %>% 
    readr::write_csv(file = "result/mmAnno/associatedPeakN_nearest.csv")

# getAssocPairNumber(mmAnno_geneScan) %>% 
#     arrange(desc(peakNumber)) %>% 
#     readr::write_csv(file = "result/mmAnno/associatedPeakN_geneScan.csv")

# mmAnno plot -------------------------------------------------------------

plot_peakGeneAlias_summary(mmAnno = mmAnno_nearest) +
    theme_bw() +
    ggtitle("nearest") -> p1

# plot_peakGeneAlias_summary(mmAnno = mmAnno_geneScan) +
#     theme_bw() +
#     ggtitle("geneScan") -> p2

pdf("plot/mmAnno/mmAnno_peakGeneAlias_summary.pdf")
p1
dev.off()


# mmAnno Cor --------------------------------------------------------------

enhancerPromoterCor(peak_GR = mergePeak_GR,
                    Txdb = Txdb,
                    up_scanPromoter = 1000,
                    down_scanPromoter = 1000,
                    up_scanEnhancer = 2e4,
                    down_scanEnhacner = 2e4,
                    peakScoreMt = peakCountNorm_ATAC) -> mmAnno_enhancerPromoter

save(mmAnno_enhancerPromoter,
     file = "result/processed_result/mmAnno_enhancerPromoter.rda")

mmAnno_enhancerPromoter_filter <- subset(mmAnno_enhancerPromoter, cor > 0.8 & pvalue < 0.01)
names(mmAnno_enhancerPromoter_filter) <- NULL
mmAnno_ePLink_df <- as.data.frame(mmAnno_enhancerPromoter_filter)
readr::write_csv(mmAnno_ePLink_df,
                 file = "result/mmAnno/mmAnno_enhancerPromoter_filter.csv")

pdf("plot/mmAnno/mmAnno_enhancerPromoterCor_summary.pdf")
plot_peakGeneAlias_summary(mmAnno_enhancerPromoter,
                           mmAnno_enhancerPromoter_filter,
                           output_type = "gene_id")
plot_peakGeneAlias_summary(mmAnno_enhancerPromoter,
                           mmAnno_enhancerPromoter_filter,
                           output_type = "feature_id")
dev.off()

getAssocPairNumber(mmAnno_enhancerPromoter_filter) %>% 
    arrange(desc(peakNumber)) %>% 
    left_join(gene_alias, by = c("gene_id" = "name")) -> a

pdf("plot/mmAnno/mmAnno_enhancerPromoterCor_singlePlot.pdf",
    width = 12)
plot_peakGeneCor(mmAnnoCor = mmAnno_enhancerPromoter_filter,
                 select_gene = "AT1G80840")
dev.off()
