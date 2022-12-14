.libPaths('/data/Blizard-Rakyan/Fran/R/packages_4.1.1/')

if(any(grepl("package:methylKit", search()))) 
    detach("package:methylKit") else message("methylKit not loaded")

if(any(grepl("package:tidyverse", search()))) 
  detach("package:tidyverse") else message("tidyverse not loaded")


#... %>% mutate(fdr=p.adjust(pval, method='fdr'))

library(tidyverse)

ifd <- 'MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/mk_wg_mask/'

md <- read_csv('MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/metadata/metadata.csv',
        col_types=cols()) %>%
    mutate(
        cond=case_when(
            Treatment=='CON_0' ~ 'c_0',
            Treatment=='CON' ~ 'c',
            Treatment=='S' ~ 's',
        Treatment=='VEH'~ 'v',
        Treatment=='DEC'~'d',
        Treatment=='S_DEC'~'s_d'),
        cond_num=case_when(
            cond=='c_0' ~ 0,
            cond=='c' ~ 1,
            cond=='s' ~ 2,
            cond=='v'~3,
            cond=='d'~4,
            cond=='s_d'~5)) %>%
    rowwise() %>%
    mutate(path=paste0(ifd, sample, '_CpG.txt'))
    

methRaw <- methylKit::methRead( 
    location = as.list(md$path),
    sample.id = as.list(md$sample),
    assembly="MM10", 
    treatment = md$cond_num)

meth_all <- methylKit::unite(methRaw, destrand = T, min.per.group = 4L)

methylKit::PCASamples(meth_all, screeplot=TRUE)



methylKit::clusterSamples(
meth_all, dist="correlation", method="ward", plot=TRUE)


meth_tiled <- methylKit::tileMethylCounts(
    meth_all, win.size=500, step.size=500, cov.bases = 10)

# methylKit::PCASamples(meth_tiled)

methylKit::PCASamples(meth_tiled, screeplot=TRUE)

methylKit::PCASamples(meth_tiled)

methylKit::clusterSamples(
meth_tiled, dist="correlation", method="ward", plot=TRUE)



#class(md$library_batch)

#md_2 <- md %>% mutate(library_batch=as.character(library_batch)) %>% select(library_batch)

#sampleAnnotation = md_2 %>% select(library_batch)

#sampleAnnotation

#as=methylKit::assocComp(mBase=meth_all, sampleAnnotation)
#as

#25.5.22 new strategy - remove PC2 to correct for batch effect
meth_2=methylKit::removeComp(meth_all,comp=2)
meth_2

methylKit::PCASamples(meth_2, screeplot=TRUE)

meth_tiled_no_pc2 <- methylKit::tileMethylCounts(
  meth_2, win.size=500, step.size=500, cov.bases = 10)

methylKit::PCASamples(meth_tiled_no_pc2, screeplot=TRUE)

methylKit::clusterSamples(
meth_2, dist="correlation", method="ward", plot=TRUE)

methylKit::clusterSamples(
meth_tiled_no_pc2, dist="correlation", method="ward", plot=TRUE)

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


#md_no_out_2 <- md %>% filter(! sample %in% c(
    'a10', 'd4', 'f6', 'd3', 'e5', 'b3', 'a3', 'c2'))

#meth_tiled_no_out_2 <- methylKit::reorganize(
    #meth_tiled,
    #sample.ids = md_no_out_2$name,
    #treatment = md_no_out_2$cond_num)

methylKit::reorganize(meth_tiled_no_pc2,
                      sample.ids = md$sample,
                      treatment = md$cond_num,
                      batch_no = md$library_batch)


#methylKit::PCASamples(meth_tiled_no_out_2)

######DMR_CSVs 

calc_diff <- function(cond_t='s', cond_c='c') {

    #md_aux <- md_no_out_2 %>%
        #mutate(cond_num=if_else(cond==cond_t, 1, 0))

    md_filt <- md %>%
        filter(cond %in% c(cond_c, cond_t))

    as(
        methylKit::calculateDiffMeth(
            methylKit::reorganize(
                meth_tiled_no_pc2,
                sample.ids = md_filt$sample,
                treatment = md_filt$cond_num),
            mc.cores=8),
        'GRanges') %>%
    as_tibble() %>%
    arrange(qvalue, desc(meth.diff)) %>%
    write_csv(
        paste0('MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/dmr_mask/',
               cond_t, '_', cond_c, '_mask_no_pc2.csv'))

    as(
        methylKit::calculateDiffMeth(
            methylKit::reorganize(
                meth_tiled_no_pc2,
                sample.ids = md_filt$sample,
                treatment = md_filt$cond_num),
            covariates = md_filt %>% select(library_batch),
            mc.cores=8),
        'GRanges') %>%
    as_tibble() %>% arrange(qvalue, desc(meth.diff)) %>%
    write_csv(
        paste0('MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/dmr_mask/',
               cond_t, '_', cond_c, '_batch_no_mask_no_pc2.csv'))

    as(
      methylKit::calculateDiffMeth(
        methylKit::reorganize(
          meth_tiled_no_pc2,
          sample.ids = md_filt$sample,
          treatment = md_filt$cond_num),
        covariates = md_filt %>% select(Treatment),
        mc.cores=8),
      'GRanges') %>%
      as_tibble() %>% arrange(qvalue, desc(meth.diff)) %>%
      write_csv(
        paste0('MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/dmr_mask/',
               cond_t, '_', cond_c, '_treatment_mask_no_pc2.csv'))
  
}

# E19

## Serum deprivation day5 - Control  day5

calc_diff()

## Control day5 - Control day 0

calc_diff(cond_t='c', cond_c='c_0')

## Serum deprivation day 5 - Control day 0

calc_diff(cond_t='s', cond_c='c_0')

## Serum deprivation_dec day 5 - serum deprviation day 5

calc_diff(cond_t='s_d', cond_c='s')

## Decitabine - Vehicle

calc_diff(cond_t='d', cond_c='v')


## Vehicle - Control day 5

calc_diff(cond_t='v', cond_c='c')

## Vehicle - control day 0 

calc_diff(cond_t='v', cond_c='c_0')

## Decitabine - Control day 5

calc_diff(cond_t='d', cond_c='c')


## serum deprivation and decitabine - decitabine 

calc_diff(cond_t='s_d', cond_c='d')

## serum deprivation and decitabine - control 

calc_diff(cond_t='s_d', cond_c='c')

## serum deprivation and decitabine - control_0 

calc_diff(cond_t='s_d', cond_c='c_0')


## Decitabine - Control day 0

calc_diff(cond_t='d', cond_c='c_0')
