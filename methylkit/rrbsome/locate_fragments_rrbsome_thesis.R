.libPaths('/data/Blizard-Rakyan/Fran/R/packages_4.1.1/')

library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)


the_sample <- args[1]

sample_summary <- read_tsv(
    paste0('/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/bamtobed_rrbsome/', the_sample, '/', the_sample, '.bed'),
    col_types=cols(),
    col_names=c('chr', 'start', 'end', 'read', 'qual', 'strand')) %>%
  group_by(chr, start, end) %>%
  summarise(num_reads=n_distinct(read), .groups='drop')

cat('Sample information summarised\n')

breaks <- read_tsv(
    paste0(
    '/data/Blizard-Rakyan/Rakyan_Lab_Files/',
    'Genomes/Mouse/WG/MM10/rDNA/WG_rDNA/mm10.rDNA.mspI.bed'),
    col_types=cols(),
    col_names=c('chr', 'start', 'end')) %>%
    group_by(chr) %>%
    mutate(tile_num=row_number()) %>%
    ungroup()

cat('Breaks read\n')

bind_rows(
    lapply(
        X=unique(sample_summary$chr),
        FUN=function(x) {
            aux <- filter(breaks, chr==x)
            the_breaks <- unique(c(aux$start, aux$end))
            the_tiles <- sample_summary %>%
                filter(chr==x) %>%
                mutate(
                    tile_s=findInterval(start, the_breaks),
                    tile_e=findInterval(end, the_breaks))
        } ) ) %>%
  inner_join(breaks, by=c('chr', 'tile_s'='tile_num'), suffix=c('', '_s')) %>%
  inner_join(breaks, by=c('chr', 'tile_e'='tile_num'), suffix=c('', '_e')) %>%
  write_csv(paste0('/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/locate_fragments/',
                   the_sample, '.csv'), quote='none')

cat('Tiles located\n')
