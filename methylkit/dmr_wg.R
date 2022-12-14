.libPaths('/data/Blizard-Rakyan/Fran/R/packages_4.1.1/')

if(any(grepl("package:methylKit", search()))) 
    detach("package:methylKit") else message("methylKit not loaded")

if(any(grepl("package:tidyverse", search()))) 
  detach("package:tidyverse") else message("tidyverse not loaded")

library(tidyverse)

ifd <- '../../data/mk_wg/'

md <- read_csv('../../data/metadata/metf_md.csv',
        col_types=cols()) %>%
    mutate(
        cond=case_when(
            cond=='CON' ~ 'c',
            cond=='OBESE' ~ 'o',
            TRUE ~ 'm'),
        cond_num=case_when(
            cond=='c' ~ 0,
            cond=='o' ~ 1,
            TRUE ~ 2) +
            if_else(age=='8w', 3, 0)) %>%
    filter(worked) %>%
    rowwise() %>%
    mutate(
        path=paste0(ifd, sample, '_CpG.txt'),
        name=paste(age, substr(cond, 1, 1),
                   substr(sex, 1, 1), sample, sep='_'))

methRaw <- methylKit::methRead( 
    location = as.list(md$path),
    sample.id = as.list(md$name),
    assembly="MM10", 
    treatment = md$cond_num)

meth_all <- methylKit::unite(methRaw, destrand=T, min.per.group = 8L)

# methylKit::PCASamples(meth_all)

# methylKit::clusterSamples(
#     meth_all, dist="correlation", method="ward", plot=TRUE)

meth_tiled <- methylKit::tileMethylCounts(
    meth_all, win.size=1000, step.size=1000, cov.bases = 10)

# methylKit::PCASamples(meth_tiled)


######

## Remove outliers

# a10, d4, f6

# md_no_out <- md %>% filter(! sample %in% c('a10', 'd4', 'f6'))

# meth_tiled_no_out <- 
#     methylKit::reorganize(
#     meth_tiled,
#     sample.ids = md_no_out$name,
#     treatment = md_no_out$cond_num)

# methylKit::PCASamples(meth_tiled_no_out)


md_no_out_2 <- md %>% filter(! sample %in% c(
    'a10', 'd4', 'f6', 'd3', 'e5', 'b3', 'a3', 'c2'))

meth_tiled_no_out_2 <- methylKit::reorganize(
    meth_tiled,
    sample.ids = md_no_out_2$name,
    treatment = md_no_out_2$cond_num)

methylKit::PCASamples(meth_tiled_no_out_2)

###### 

calc_diff <- function(the_age='E19', cond_t='o', cond_c='c') {

    md_aux <- md_no_out_2 %>%
        mutate(cond_num=if_else(cond==cond_t, 1, 0))

    md_filt <- md_aux %>%
        filter(age==the_age, cond %in% c(cond_c, cond_t))

    as(
        methylKit::calculateDiffMeth(
            methylKit::reorganize(
                meth_tiled,
                sample.ids = md_filt$name,
                treatment = md_filt$cond_num),
            mc.cores=8),
        'GRanges') %>%
    as_tibble() %>%
    arrange(qvalue, desc(meth.diff)) %>%
    write_csv(
        paste0('../../data/dmr/', the_age, '_',
               cond_t, '_', cond_c, '.csv'))

    as(
        methylKit::calculateDiffMeth(
            methylKit::reorganize(
                meth_tiled,
                sample.ids = md_filt$name,
                treatment = md_filt$cond_num),
            covariates = md_filt %>% select(sex),
            mc.cores=8),
        'GRanges') %>%
    as_tibble() %>% arrange(qvalue, desc(meth.diff)) %>%
    write_csv(
        paste0('../../data/dmr/', the_age, '_',
               cond_t, '_', cond_c, '_sex.csv'))

    md_males <- md_aux %>%
        filter(age==the_age, cond %in% c(cond_c, cond_t), sex=='male')

    as(
       methylKit::calculateDiffMeth(
            methylKit::reorganize(
                meth_tiled,
                sample.ids = md_males$name,
                treatment = md_males$cond_num),
            mc.cores=8),
        'GRanges') %>%
        as_tibble() %>% arrange(qvalue, desc(meth.diff)) %>%
        write_csv(
            paste0('../../data/dmr/', the_age, '_', 
               cond_t, '_', cond_c, '_males.csv'))

    md_females <- md_aux %>%
        filter(age==the_age, cond %in% c(cond_c, cond_t), sex=='female')

    as(
       methylKit::calculateDiffMeth(
            methylKit::reorganize(
                meth_tiled,
                sample.ids = md_females$name,
                treatment = md_females$cond_num),
            mc.cores=8),
        'GRanges') %>%
        as_tibble() %>% arrange(qvalue, desc(meth.diff)) %>%
        write_csv(
            paste0('../../data/dmr/', the_age, '_',
                cond_t, '_', cond_c, '_females.csv'))

}

# E19

## Obese - Control

calc_diff()

## Obese - Metformin

calc_diff(cond_t='o', cond_c='m')

## Metformin - Control

calc_diff(cond_t='m', cond_c='c')

# 8w

## Obese - Control

calc_diff(the_age='8w')

## Obese - Metformin

calc_diff(the_age='8w', cond_t='o', cond_c='m')

## Metformin - Control

calc_diff(the_age='8w', cond_t='m', cond_c='c')
