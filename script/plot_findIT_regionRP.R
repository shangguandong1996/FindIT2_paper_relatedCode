# Info ---------------------------
## Script name: plot_findIT_regionRP.R

## experiment_name:

## Purpose of script:
## 1. XX
## 2. XX

## Author: Guandong Shang
## Date Created: 2021-07-12
## Copyright (c) Guandong Shang, 2021
## Email: shangguandong1996@163.com

# Notes ---------------------------
## 1. XX
## 2. XX   
##

# Prepare -----------------------------------------------------------------

# load up the packages
library(ggplot2)
library(dplyr)

library(FindIT2)

# Set Options
options(stringsAsFactors = F)

# load up the data
load("result/FindIT/FindIT2_regionRP_merge_result.rda")


# plot point --------------------------------------------------------------

TF_pvalue <- assays(merge_result)$TF_pvalue
TF_pvalue[, 1] %>% 
    as_tibble(rownames = "TF_id") %>% 
    mutate(rank = rank(value)) -> plot_data
plot_data %>% 
    ggplot(aes(x = rank, y = -log10(value))) +
    geom_point() +
    ggrepel::geom_text_repel(data = filter(plot_data, TF_id == "LEC2"),
                             label = "LEC2",
                             min.segment.length = unit(0, 'lines'), 
                             nudge_y = .2) +
    theme_bw() +
    ggsave("plot/findTF/findIT_regionRP_point.pdf")


FindIT2::shinyParse_findIT_regionRP(merge_result)
