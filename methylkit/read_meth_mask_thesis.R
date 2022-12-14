.libPaths('/data/Blizard-Rakyan/Fran/R/packages_4.1.1/')

library(tidyverse)

ifd <- '/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/bam_mask_mouse/'

md <- read_csv('/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/metadata/metadata.csv',
               col_types=cols()) %>%
  mutate(
    cond_num=case_when(
      Treatment=='CON' ~ 0,
      Treatment=='CON_0' ~ 1,
      Treatment=='S' ~ 2,
Treatment=='VEH' ~ 3,
Treatment=='DEC' ~ 4,
Treatment=='S_DEC' ~ 5)) %>%
  rowwise() %>%
  mutate(
    path=paste0(
      ifd, sample, '/', list.files(paste0(ifd, sample, '/'),
                                   pattern=".bam$")))

methRaw <- methylKit::processBismarkAln( 
  location = as.list(md$path),
  sample.id = as.list(md$sample),
  assembly = "MM10", 
  read.context = "CpG", 
  save.folder = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/mk_wg_mask/",
  treatment = md$cond_num)
