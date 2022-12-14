.libPaths('/data/Blizard-Rakyan/Selin/R/packages/')

library(tidyverse)


fd <- '/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/bme_30_51/'


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
write_csv('meth_18s_30_51.csv')



library(readr)
meth_18s_30_51 <- read_csv("meth_18s_30_51.csv")
meth_18s_30_51


meth_18s_30_51_replicate <- meth_18s_30_51 %>%
  mutate(Replicate = case_when(
    startsWith(sample, "N1") ~ "1",
    startsWith(sample, "N2") ~ "2",
    startsWith(sample, "N3") ~ "3",
    startsWith(sample, "N4") ~ "4"
  ))
meth_18s_30_51_replicate

meth_18s_30_51_replicate_treatment <- meth_18s_30_51_replicate %>%
  mutate(Treatment = case_when(
    endsWith(sample, "_0") ~ "CON_0",
    endsWith(sample, "_CON") ~ "CON",
    endsWith(sample, "_S") ~ "S_DEP",
    endsWith(sample, "_DEC") ~ "DEC",
    endsWith(sample, "_S_D") ~ "S_DEP_DEC",
    endsWith(sample, "_VEH") ~ "VEH",
    
    
  ))

meth_18s_140_160_replicate_treatment

meth_18s_140_160_replicate_treatment <- rename(meth_18s_140_160_replicate_treatment, meth_prop_140_160 = "meth_prop")


library(forcats)
library(ggplot2)
rDNA_meth_S_C <- meth_18s_140_160_replicate_treatment %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% ggplot() +
  aes(x = Treatment, y = meth_prop_140_160) +
  geom_point(aes(colour = factor(Treatment)),  size = 2, alpha = 0.8) +
  ggpubr::stat_compare_means(aes(label = ..p.signif..),
                             paired = TRUE, ref.group = "CON_0", hide.ns=T,
                             bracket.size=1, label.y=0.3) +
  theme_bw() +
  labs(x="Treatment condition", y="Proportion of 18S methylation") +
  coord_cartesian(ylim=c(0, 0.3)) +
  ggtitle("18S rDNA methylation in serum deprived MEFs")+
  scale_colour_discrete('Treatment Condition')
rDNA_meth_S_C

library(readr)
meth_18s <- read_csv("MEF_SD_DEC/RRBS/Data_Novogene/scripts/bme/meth_18s.csv")
meth_18s

meth_18_both <- inner_join(meth_18s, meth_18s_140_160_replicate_treatment, by = "sample")



library(ggplot2)
meth_corr <- meth_18_both %>% 
  ggplot() +
  aes(x = meth_prop, y = meth_prop_140_160) + 
  geom_point(aes(colour = Treatment), alpha = 1) +
  theme_bw()+
  labs(x = "rDNA methylation before window filter", y = "rDNA methylation after 140-160 window filter") + 
  ggpubr::stat_cor() +
  ggtitle("Methylation data correlation")
meth_corr
