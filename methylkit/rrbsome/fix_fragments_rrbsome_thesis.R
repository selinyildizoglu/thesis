.libPaths('/data/Blizard-Rakyan/Fran/R/packages_4.1.1/')

library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

the_sample <- args[1]

locations <- read_csv(
  paste0('/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/locate_fragments/',
  the_sample, '.csv'), col_types=cols())


find_likely <- function(
  start_read, end_read,
  start_st, end_st,
  start_et, end_et, threshold=10) {
  if (start_st == start_et) return(1)
  ov <- c()
  ov[1] <- min(end_read, end_st) - max(start_read, start_st)
  ov[2] <- min(end_read, end_et) - max(start_read, start_et)
  ov[3] <- ifelse(end_st != start_et,
                  min(end_read, start_et) - max(start_read, end_st), 0)
  max_pos <- as.numeric(which(max(ov)==ov))
  case_when(
    sum(ov > threshold) != 1 ~ 0,
    length(max_pos) > 1 ~ 0,
    TRUE ~ max_pos[1]
  )
}

locations$maj <- locations %>%
  select(start_read=start, end_read=end,
         start_st=start_s, end_st=end_s,
         start_et=start_e, end_et=end_e) %>%
  purrr::pmap_dbl(find_likely, threshold=10)

ofd <- paste0('/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/fix_fragments/', the_sample, '/')
dir.create(ofd, showWarnings = FALSE, recursive = TRUE)

locations_read <- locations %>% 
  filter(maj != 0, maj != 3 | tile_s - tile_e < 3) %>%
  mutate(
    start_f=case_when(
      maj==3 ~ end_s,
      maj==2 ~ start_e,
      TRUE ~ start_s),
    end_f=case_when(
      maj==3 ~ start_e,
      maj==2 ~ end_e,
      TRUE ~ end_s )) 
  
locations_read %>%
  select(chr, start, end, start_f, end_f) %>%
  write_csv(paste0(ofd, the_sample, '_read.csv'), quote='none')

locations_sum <- locations_read  %>%
  group_by(chr, start_f, end_f) %>%
  summarise(num_reads=sum(num_reads), .groups='drop')

locations_sum %>%
  write_csv(paste0(ofd, the_sample, '_unfilt.csv'), quote='none')

locations_sum %>%
  filter(end_f - start_f >= 25, end_f - start_f <= 250) %>%
  filter(num_reads >= 5) %>%
  write_csv(paste0(ofd, the_sample, '_filt_5x.csv'), quote='none')

locations_sum %>%
  filter(end_f - start_f >= 25, end_f - start_f <= 250) %>%
  filter(num_reads >= 10) %>%
  write_csv(paste0(ofd, the_sample, '_filt_10x.csv'), quote='none')
