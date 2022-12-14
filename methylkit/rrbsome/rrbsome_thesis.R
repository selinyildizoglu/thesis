library(tidyverse)

metf_fg_fd <- '../../data/fix_fragments/'

metf_fg_fs <- list.files(metf_fg_fd, pattern='10x')

metf_fg <- bind_rows(
  lapply(X=metf_fg_fs,
         FUN=function(x) { 
           read_csv(
             paste0(metf_fg_fd, x),
             col_types=cols()) %>%
           mutate(sample=gsub('.csv', '', x)) } ))

metf_fg_sum <- metf_fg %>%
  filter(! grepl("BK", chr)) %>%
  group_by(chr, start_f, end_f) %>%
  summarise(num_samples=n_distinct(sample),
            .groups='drop')


dir.create('../../data/rrbsome/', showWarnings = F)


#### SET A MINIMUM NUMBER OF SAMPLES REQUIREMENT FOR YOUR DATASET

min_samples <- 

metf_fg_sum %>%
  filter(num_samples >= min_samples) %>%
  select(-num_samples) %>%
  write_tsv('../../data/rrbsome/rrbsome_10x_24min.bed',
            col_names=F)

