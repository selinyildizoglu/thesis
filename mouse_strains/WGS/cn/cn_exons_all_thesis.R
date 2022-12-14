.libPaths('/data/Blizard-Rakyan/Selin/R/packages')

library(tidyverse)

exome_fd='/data/Blizard-Rakyan/Selin/projects/mouse_strains/WGS_thesis/data/covs_exons_all/'

bind_rows(
    lapply(X=list.files(exome_fd),
           FUN=function(x) {
               read_tsv(
                    paste0(exome_fd, x),
                           col_types=cols(), col_names=c('chr', 'pos', 'cov')) %>%
                   summarise(sample=gsub(".cov", "", x),
                             exome_cov_all=mean(cov))
           }
)) %>% write_csv('cov_exome_all_WGS_thesis.csv')
