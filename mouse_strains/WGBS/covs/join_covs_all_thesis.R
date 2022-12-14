.libPaths('/data/Blizard-Rakyan/Selin/R/packages/')

library(tidyverse)
fd <- '/data/Blizard-Rakyan/Selin/projects/mouse_strains/WGBS_thesis/data/covs_rdna_all/'
fns <- list.files(fd)

covs <- do.call(
        rbind,
        lapply(
            X = fns,
            FUN = function(x) {
                read_tsv(
                    paste0(fd, x),
                    col_names=c('ref', 'pos', 'cov'),
                    col_types=cols()) %>%
                mutate(sample=gsub('.cov', '', x)) %>%
                select(sample, pos, cov)
                }))

write_csv(covs, 'covs_all.csv')
