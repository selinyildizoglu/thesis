.libPaths('/data/Blizard-Rakyan/Selin/R/packages/')

library(tidyverse)
fd <- '../../data/covs_hoxa/'
fns <- dir(fd)

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

write_csv(covs, 'covs_hoxa.csv')
