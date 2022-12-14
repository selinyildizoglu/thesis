.libPaths('/data/Blizard-Rakyan/Fran/R/packages_4.1.1/')

library(tidyverse)

ifd <- '../../data/bismark/'

md <- read_csv('../../data/metadata/metf_md.csv',
               col_types=cols()) %>%
  mutate(
    cond_num=case_when(
      cond=='CON' ~ 0,
      cond=='OBESE' ~ 1,
      TRUE ~ 2)) %>%
  filter(worked) %>%
  rowwise() %>%
  mutate(
    path=paste0(
      ifd, sample, '/', list.files(paste0(ifd, sample, '/'),
                                   pattern=".bam$")))

methRaw <- methylKit::processBismarkAln( 
  location = as.list(md$paths),
  sample.id = as.list(md$samples),
  assembly = "MM10", 
  read.context = "CpG", 
  save.folder = "../../data/mk_wg/",
  treatment = md$cond_num)
