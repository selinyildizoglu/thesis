.libPaths('/data/Blizard-Rakyan/Fran/R/packages_4.0.2/')

library(tidyverse)

fd <- '../../data/bme/'

bind_rows(
    lapply(
    X=list.files(fd),
    FUN=function(x) {
        read_tsv(
            paste0(fd, x, '/', x, '.bismark.cov.gz'),
            col_names=c('contig', 'start', 'end', 'meth_prop', 'meth', 'unmeth'),
            col_types=cols()) %>%
        filter(meth + unmeth > 100, start >= 4008+3009, end <=5877+3009) %>%
        summarise(meth_prop=sum(meth)/(sum(meth) + sum(unmeth))) %>%
        mutate(sample=x) %>%
        select(sample, meth_prop)})) %>%
write_csv('meth_18s.csv')
