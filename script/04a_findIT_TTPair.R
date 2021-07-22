# Info ---------------------------
## Script name: findIT_TTPair.R

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
public_database <- readr::read_delim("rawdata/related_annotation/iGRN_network_full.txt",
                                     delim = "\t") %>% 
    select(1:2) %>% 
    rename(TF_id = `#TF`)

load("result/processed_result/LEC2_targetTop1000.rda")

gene_alias <- readr::read_csv("rawdata/related_annotation/gene_aliases_20190630_tidy.csv")

dir.create("result/influential_TF")

# TTpair ------------------------------------------------------------------

TTPair_result <- findIT_TTPair(input_genes = LEC2_target,
                               TF_target_database = public_database)

TTPair_result %>% 
    left_join(gene_alias, by = c("TF_id" = "name")) -> TTPair_result
readr::write_csv(TTPair_result, "result/influential_TF/TTPair_public.csv")
