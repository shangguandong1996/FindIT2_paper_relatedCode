# Info ---------------------------
## Script name: plot_findIT_summary.R

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
library(dplyr)
library(ggplot2)

# Set Options
options(stringsAsFactors = F)


# merge_result ------------------------------------------------------------

purrr::map(list.files("result/influential_TF", pattern = "_(motif|public).csv",
                      full.names = TRUE), function(i){    
    readr::read_csv(i) %>% 
        select(TF_id, rank) %>% 
        mutate(LEC2 = TF_id %in% c("LEC2","AT1G28300"),
               name = gsub(".csv", "", basename(i)),
               rank_zscore = FindIT2:::INT(-rank)) %>% 
        select(rank_zscore, LEC2, rank, name)
}) %>% 
    bind_rows() -> merge_data

ggplot(merge_data, aes(x = name, y = rank_zscore)) +
    ggrastr::geom_quasirandom_rast() +
    geom_point(data = filter(merge_data, LEC2),
               color = "#d94d40") +
    ggrepel::geom_text_repel(data = filter(merge_data, LEC2),
                             label = "LEC2",
                             min.segment.length = unit(0, 'lines'),
                             nudge_y = 0.5) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    ggsave("plot/LEC2_rank_quasirandom.pdf",
           width = 6, height = 4)


ggplot(filter(merge_data, LEC2), aes(x = name, y = rank, group = LEC2)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    geom_text(aes(label = rank), vjust = -0.5) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    ggsave("plot/LEC2_rank_trend.pdf",
           width = 6, height = 5)
