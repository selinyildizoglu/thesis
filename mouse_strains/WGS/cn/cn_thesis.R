.libPaths('/data/Blizard-Rakyan/Selin/R/packages')

rdna_fd <- '/data/Blizard-Rakyan/Selin/projects/mouse_strains/WGS_thesis/data/covs_rdna/'
exons_fd <- '/data/Blizard-Rakyan/Selin/projects/mouse_strains/WGS_thesis/data/covs_exons/'

library(tidyverse)

adjust_pos <- function(pos, loop_size) {
      if_else(pos<=loop_size, pos-loop_size-1, pos-loop_size)
}

max_length <- 101

regions <- tribble(
        ~name, ~start, ~end,
        "18S", 4008, 5877,
        "28S", 8123, 12852) 
loop_size <- 3008

depths <- data.frame()

for (sample in dir(exons_fd)) {

    sample <- gsub('.cov', '', sample)
    cat(sample)
    
    exons_depth <- read.csv(paste0(exons_fd, sample, '.cov'),
        sep ='\t', header=F, 
        col.names = c('name', 'pos', 'depth'),
        stringsAsFactors = F) %>%
        group_by(name) %>%
        filter(pos > max_length,
        pos < (max(pos) - max_length)) %>%
        ungroup() %>%
        summarise(mean_depth = mean(depth)) %>%
        pull(mean_depth)

    rdna <- read.csv(paste0(rdna_fd, sample, '.cov'),
        sep ='\t', header=F, 
        col.names = c('name', 'pos', 'depth'),
        stringsAsFactors = F) %>%
        mutate(pos = adjust_pos(pos, loop_size))

    rdna_depth <- rdna %>%
        filter(
            pos >= (
                regions %>% filter(name == '18S'))$start,
            pos <= (
                regions %>% filter(name == '28S'))$end
        ) %>%
        summarise(mean_depth = mean(depth)) %>%
        pull(mean_depth)    

    rdna_18s_depth <- rdna %>%
        filter(
            pos >= (
                regions %>% filter(name == '18S'))$start,
            pos <= (
                regions %>% filter(name == '18S'))$end
        ) %>%
        summarise(mean_depth = mean(depth)) %>%
        pull(mean_depth)

    rdna_28s_depth <- rdna %>%
        filter(
            pos >= (
                regions %>% filter(name == '28S'))$start,
            pos <= (
                regions %>% filter(name == '28S'))$end
        ) %>%
        summarise(mean_depth = mean(depth)) %>%
        pull(mean_depth)

    depths <- rbind(depths,
        data.frame(sample = sample,
            exons = exons_depth,
            rdna = rdna_depth,
            rdna_18s = rdna_18s_depth,
            rdna_28s = rdna_28s_depth))    

    cat('\n')

}

depths_cn <- depths %>%
    mutate(cn_18s = rdna_18s / exons,
        cn_28s = rdna_28s / exons)

dir.create('/data/Blizard-Rakyan/Selin/projects/mouse_strains/WGS_thesis/data/cn/', showWarnings = FALSE)

write_csv(depths_cn, '/data/Blizard-Rakyan/Selin/projects/mouse_strains/WGS_thesis/data/cn/cn_WGS_thesis.csv')

