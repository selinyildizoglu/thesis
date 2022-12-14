.libPaths('/data/Blizard-Rakyan/Fran/R/packages_4.1.1/')

library(tidyverse)

rrbsome_fg_fd <- '/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/fix_fragments/'

rrbsome_fg_fd_samples <- list.files(rrbsome_fg_fd)

#rrbsome_fg_fs <- list.files(paste0(rrbsome_fg_fd, '/', rrbsome_fg_fd_samples, '/'), pattern='10x')

rrbsome_fg <- bind_rows(
  lapply(X=rrbsome_fg_fd_samples,
         FUN=function(x) { 
           read_csv(
             paste0(rrbsome_fg_fd, x, '/', x, '_filt_10x.csv'),
             col_types=cols()) %>%
           mutate(sample=gsub('.csv', '', x)) } ))


rrbsome_fg_sum <- rrbsome_fg %>%
  filter(! grepl("BK", chr)) %>%
  group_by(chr, start_f, end_f) %>%
  summarise(num_samples=n_distinct(sample),
            .groups='drop')


dir.create('/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/rrbsome/', showWarnings = F)


#### SET A MINIMUM NUMBER OF SAMPLES REQUIREMENT FOR YOUR DATASET

min_samples <- 24

rrbsome_fg_sum %>%
  filter(num_samples == min_samples) %>%
  select(-num_samples) %>%
  write_tsv('/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/rrbsome/rrbsome_10x_24min.bed',
            col_names=F)

