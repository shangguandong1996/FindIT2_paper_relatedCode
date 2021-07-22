# Info ---------------------------
## Script name: download_database.R

## experiment_name:

## Purpose of script:
## 1. XX
## 2. XX

## Author: Guandong Shang
## Date Created: 2021-07-21
## Copyright (c) Guandong Shang, 2021
## Email: shangguandong1996@163.com

# Notes ---------------------------
## 1. XX
## 2. XX   
##

# Prepare -----------------------------------------------------------------


# download FindIT_TTpari database
# http://bioinformatics.psb.ugent.be/webtools/iGRN/pages/download
iGRN_database <- "http://bioinformatics.psb.ugent.be/webtools/iGRN/files/igrn_data/iGRN_network_full.txt"
download.file(iGRN_database,
              destfile = "rawdata/related_annotation/iGRN_network_full.txt")

# remap At 2020
# https://remap.univ-amu.fr/download_page
At_TF_database <- "https://remap.univ-amu.fr/storage/remap2020/tair10/tf/MACS2/remap2020_all_macs2_TAIR10_v1_0.bed.gz"
download.file(At_TF_database,
              destfile = "rawdata/related_annotation/remap2020_all_macs2_TAIR10_v1_0.bed.gz")
