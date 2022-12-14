.libPaths('/data/Blizard-Rakyan/Selin/R/packages/')

library(tidyverse)

the_sample <- commandArgs(trailingOnly=TRUE)[1]

mouse_rdna_mspI_fragments <- read_csv(
    paste0('/data/Blizard-Rakyan/Rakyan_Lab_Files/Genomes/Mouse/rDNA/BK000964/Looped/',
           'BK000964_looped_3008_mspI.csv'),
        col_types=cols())

fd <- paste0('/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/bamtobed_rdna/', the_sample, '/')

beds <- read_tsv(
            paste0(fd, the_sample, '.bed'),
                col_names=c('chr', 'start', 'end', 'read_id', 'score', 'strand'),
                col_types=cols()) %>%
            separate(read_id, into=c('read_id', 'mate'), sep='/') %>%
            group_by(read_id) %>%
            summarise(start=min(start), end=max(end)) %>%
            mutate(read_num=row_number(read_id)) %>%
            select(-read_id) %>%
            mutate(sample=the_sample) %>%
            relocate(sample)

distances <- beds %>%
    mutate(start_g=cut(start+1,
            breaks = c(0, mouse_rdna_mspI_fragments$end, Inf),
                                    include.lowest = F),
    end_g=cut(end-1,
                    breaks = c(0, mouse_rdna_mspI_fragments$end, Inf),
                                        include.lowest = F)) %>%
    mutate(start_g=gsub('\\(', '', gsub(']', '', start_g))) %>%
    mutate(end_g=gsub('\\(', '', gsub(']', '', end_g))) %>%
    separate(start_g, into=c('start_g_s', 'start_g_e'), sep=',', convert=T) %>%
    separate(end_g, into=c('end_g_s', 'end_g_e'), sep=',', convert=T) %>%
    rowwise() %>%
    mutate(dist_s=min(abs(start-start_g_s), abs(start-end_g_s)),
    dist_e=min(abs(end-start_g_e), abs(end-end_g_e))) %>%
    ungroup()

dir.create('.tmp', showWarnings=F)

distances %>%
    mutate(susp=(dist_s > 6 | dist_e > 6)) %>%
    group_by(sample) %>%
    summarise(prop_susp=mean(susp)) %>%
    write_csv(paste0('.tmp/', the_sample ,'.csv'))


