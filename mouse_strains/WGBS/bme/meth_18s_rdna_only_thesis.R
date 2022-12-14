.libPaths('/data/Blizard-Rakyan/Selin/R/packages/')

library(tidyverse)

fd <- '/data/Blizard-Rakyan/Selin/projects/mouse_strains/WGBS_thesis/data/bme/'

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
