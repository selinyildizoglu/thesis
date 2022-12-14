.libPaths('/data/Blizard-Rakyan/Fran/R/packages_4.1.1/')

if(any(grepl("package:methylKit", search()))) 
  detach("package:methylKit") else message("methylKit not loaded")

detach("package:biomaRt")

###haplotypes mouse strains
common_path_m = "/data/Blizard-Rakyan/Selin/projects/mouse_strains/WGBS/data/blink/"
sub_dirs = list.files(common_path_m)
sub_dirs

link_table_m <- bind_rows(lapply(
  X=sub_dirs,
  FUN=function(x){read_csv(
    paste0(common_path_m, x, '/', x, '_link.csv'),
    col_types=cols()) %>%
      mutate(sample=x)}))
link_table_m

write.csv(link_table_m, "/data/Blizard-Rakyan/Selin/projects/mouse_strains/WGBS/csv_tables/link_table_m.csv", row.names = FALSE)

library(readxl)
samples_ID_strain <- read_excel("mouse_strains/WGBS/csv_tables/samples_ID_strain.xlsx")
samples_ID_strain

link_table_m_strain <- inner_join(link_table_m, samples_ID_strain, by = "sample")

write.csv(link_table_m_strain, "/data/Blizard-Rakyan/Selin/projects/mouse_strains/WGBS/csv_tables/link_table_m_strain.csv", row.names = FALSE)

#haplotype mouse strain table

link_table_BL6J <-  link_table_m_strain %>% filter(Strain == "BL6J")

library(readr)
haplotype_specific_alleles <- read_csv("MEF_SD_DEC/RRBS/Data_Novogene/scripts/blink/haplotype_specific_alleles.csv")

haplotypes_mef <- haplotype_specific_alleles %>%
  group_by(pos_ref, allele, dominant) %>%
  mutate(num_hap=n(), haplotype=if_else(num_hap==1, haplotype, 'other')) %>%
  ungroup() %>%
  select(pos_ref, allele, haplotype) %>%
  distinct()


adjust_haplo <- function(haplotype)  {
  case_when(
    haplotype=='A_G_T_A' ~ 'ATA',
    haplotype=='A_G_T_G' ~ 'ATG',
    haplotype=='C_T_C_A' ~ 'CCA',
    haplotype=='C_T_T_A' ~ 'CTA',
    TRUE ~ haplotype)
}

library(dplyr)
link_table_BL6J_adj <- link_table_BL6J %>%
  mutate(adj_snp_pos=if_else(snp_pos<= 3008, snp_pos-3009, snp_pos-3008))
link_table_BL6J_adj

link_table_BL6J_adj_ = subset(link_table_BL6J_adj, select = -snp_pos )

link_table_BL6J_adj_hap <- link_table_BL6J_adj_ %>%
  inner_join(haplotypes_mef %>%
               mutate(haplotype=adjust_haplo(haplotype)),
             by=c('adj_snp_pos'='pos_ref',
                  'snp_allele'='allele'))


link_table_BL6J_filtered <- link_table_BL6J_adj_hap %>%
  group_by(sample, adj_snp_pos, snp_allele, meth_pos, haplotype) %>%
  filter(n()>=10) %>%
  group_by(sample, adj_snp_pos, meth_pos) %>%
  filter(n()>=50, n_distinct(haplotype)>1) %>%
  group_by(sample, haplotype, ID) %>%
  summarise(meth=sum(meth_state=='+')/n(),
            .groups='drop')

write.csv(link_table_BL6J_filtered, "/data/Blizard-Rakyan/Selin/projects/mouse_strains/WGBS/csv_tables/link_table_BL6J_filtered.csv", row.names = FALSE)

##BL6J only snps table
common_path2_m = "/data/Blizard-Rakyan/Selin/projects/mouse_strains/WGBS/data/blink/"
sub_dirs2 = list.files(common_path2_m)

snps_table_m <- bind_rows(lapply(
  X=sub_dirs2,
  FUN=function(x){
    read_csv(
      paste0(common_path2_m, x, '/', x, '_snps.csv'),
      col_types=cols()) %>%
      mutate(sample=x, alt_prop_exp = 1 - ref_prop_exp,alt_prop_obs = 1 - ref_prop_obs)})) 
snps_table_m


snps_table_m_strain <- inner_join(snps_table_m, samples_ID_strain, by = "sample")

write.csv(snps_table_m_strain, "/data/Blizard-Rakyan/Selin/projects/mouse_strains/WGBS/csv_tables/snps_table_m_strain.csv", row.names = FALSE)

snps_table_BL6J <- snps_table_m_strain %>% filter(Strain == "BL6J")

library(dplyr)
snps_table_BL6J_adj <- snps_table_BL6J %>%
  mutate(adj_snp_pos=if_else(snp_pos<= 3008, snp_pos-3009, snp_pos-3008))
snps_table_BL6J_adj

write.csv(snps_table_BL6J_adj, "/data/Blizard-Rakyan/Selin/projects/mouse_strains/WGBS/csv_tables/snps_table_BL6J_adj.csv", row.names = FALSE)


library(dplyr)
link_table_BL6J_adj <- link_table_BL6J %>%
  mutate(adj_snp_pos=if_else(snp_pos<= 2120, snp_pos-2121, snp_pos-2120))
link_table_BL6J_adj

#######ch3 human

common_path_h = "/data/Blizard-Rakyan/Robert/Data/Human_Mandinka/Real_Data/WGBS/data/blink/"
sub_dirs_h = list.files(common_path_h)
sub_dirs_h

link_table_h <- bind_rows(lapply(
  X=sub_dirs_h,
  FUN=function(x){read_csv(
    paste0(common_path_h, x, '/', x, '_link.csv'),
    col_types=cols()) %>%
      mutate(sample=x)}))
link_table_h


mandinka_meth <- link_table_h %>%
  group_by(sample, snp_pos, snp_allele, meth_pos) %>%
  filter(n()>=10) %>%
  group_by(sample, snp_pos, meth_pos) %>%
  filter(n()>=50, n_distinct(snp_allele)>1) %>%
  group_by(sample, snp_pos, snp_allele, meth_pos) %>%
  summarise(meth=sum(meth_state=='+')/n(), .groups='drop')

write.csv(mandinka_meth, "/data/Blizard-Rakyan/Selin/projects/mouse_strains/WGBS/csv_tables/human/mandinka_meth.csv", row.names = FALSE)

common_path_2h = "/data/Blizard-Rakyan/Robert/Data/Human_Mandinka/Real_Data/WGBS/blink/"
sub_dirs_2h = list.files(common_path_h)

snps_table_h <- bind_rows(lapply(
  X=sub_dirs2,
  FUN=function(x){
    read_csv(
      paste0(common_path2_m, x, '/', x, '_snps.csv'),
      col_types=cols()) %>%
      mutate(sample=x, alt_prop_exp = 1 - ref_prop_exp,alt_prop_obs = 1 - ref_prop_obs)})) 
snps_table_h




#2120





#################################################### MEF_SD_DEC

common_path = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/blink_nd/"
sub_dirs = list.files(common_path)
sub_dirs

link_table_nd <- bind_rows(lapply(
  X=sub_dirs,
  FUN=function(x){read_csv(
    paste0(common_path, x, '/', x, '_link.csv'),
    col_types=cols()) %>%
      mutate(sample=x)}))
link_table_nd

write.csv(link_table_nd, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/link_table_nd.csv", row.names = FALSE)

common_path_i = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/blink/"

link_table_i <- bind_rows(lapply(
  X=sub_dirs,
  FUN=function(x){read_csv(
    paste0(common_path_i, x, '/', x, '_link.csv'),
    col_types=cols()) %>%
      mutate(sample=x)}))
link_table_i

write.csv(link_table_i, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/link_table_i.csv", row.names = FALSE)


#8.7.22 colour blind palette

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# To use for fills, add
scale_fill_manual(values=cbPalette)

# To use for line and point colors, add
scale_colour_manual(values=cbPalette)



#7.7.22


library(readr)
haplotype_specific_alleles <- read_csv("MEF_SD_DEC/RRBS/Data_Novogene/scripts/blink/haplotype_specific_alleles.csv")

haplotypes_mef <- haplotype_specific_alleles %>%
  group_by(pos_ref, allele, dominant) %>%
  mutate(num_hap=n(), haplotype=if_else(num_hap==1, haplotype, 'other')) %>%
  ungroup() %>%
  select(pos_ref, allele, haplotype) %>%
  distinct()


adjust_haplo <- function(haplotype)  {
  case_when(
    haplotype=='A_G_T_A' ~ 'ATA',
    haplotype=='A_G_T_G' ~ 'ATG',
    haplotype=='C_T_C_A' ~ 'CCA',
    haplotype=='C_T_T_A' ~ 'CTA',
    TRUE ~ haplotype)
}

library(dplyr)
link_table_i_adj <- link_table_i %>%
  mutate(adj_snp_pos=if_else(snp_pos<= 3008, snp_pos-3009, snp_pos-3008))
link_table_i_adj

link_table_i_adj = subset(link_table_i_adj, select = -snp_pos )

link_table_i_adj <- link_table_i_adj %>%
  inner_join(haplotypes_mef %>%
               mutate(haplotype=adjust_haplo(haplotype)),
             by=c('adj_snp_pos'='pos_ref',
                  'snp_allele'='allele'))

link_i_Treatment <- link_table_i_adj %>%
  mutate(Treatment = case_when(
    endsWith(sample, "_0") ~ "CON_0",
    endsWith(sample, "_CON") ~ "CON",
    endsWith(sample, "_VEH") ~ "VEH",
    endsWith(sample, "_DEC") ~ "DEC",
    endsWith(sample, "_S") ~ "S_DEP",
    endsWith(sample, "_S_D") ~ "S_DEP_DEC"
  ))

link_i_Treatment

link_i_Treatment_replicate <- link_i_Treatment %>%
  mutate(Replicate = case_when(
    startsWith(sample, "N1") ~ "1",
    startsWith(sample, "N2") ~ "2",
    startsWith(sample, "N3") ~ "3",
    startsWith(sample, "N4") ~ "4"
  ))

link_i_Treatment_replicate



link_table_i_filtered <- link_i_Treatment_replicate %>%
  group_by(sample, adj_snp_pos, snp_allele, meth_pos, haplotype) %>%
  filter(n()>=10) %>%
  group_by(sample, adj_snp_pos, meth_pos, Treatment, Replicate) %>%
  filter(n()>=50, n_distinct(haplotype)>1) %>%
  group_by(sample, haplotype, Treatment, Replicate) %>%
  summarise(meth=sum(meth_state=='+')/n(),
            .groups='drop')

library(forcats)
library(ggplot2)
rDNA_DEC_ATA_meth <- link_table_i_filtered %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% 
  filter(haplotype == 'ATA') %>% ggplot() +
  aes(x = Treatment, y = meth) +
  geom_jitter(aes(colour = Treatment, shape = Replicate),  size = 2) +
  theme_bw() +
  labs(x="Treatment", y= "Methylation level") +
  ggpubr::stat_cor(size = 3) +
  ggtitle("ATA methylation level in Decitabine-exposed MEFs") +
  scale_colour_discrete('Treatment Condition') +
  scale_shape_discrete('Replicate')
rDNA_DEC_ATA_meth

library(forcats)
library(ggplot2)
rDNA_DEC_ATG_meth <- link_table_i_filtered %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% 
  filter(haplotype == 'ATG') %>% ggplot() +
  aes(x = Treatment, y = meth) +
  geom_point(aes(colour = Treatment, shape = Replicate),  size = 2) +
  theme_bw() +
  labs(x="Treatment", y= "Methylation level") +
  ggpubr::stat_cor(size = 3) +
  ggtitle("ATG methylation level in Decitabine-exposed MEFs") +
  scale_colour_discrete('Treatment Condition') +
  scale_shape_discrete('Replicate') 
rDNA_DEC_ATG_meth

library(forcats)
library(ggplot2)
rDNA_DEC_CTA_meth <- link_table_i_filtered %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% 
  filter(haplotype == 'CTA') %>% ggplot() +
  aes(x = Treatment, y = meth) +
  geom_point(aes(colour = Treatment, shape = Replicate),  size = 2) +
  theme_bw() +
  labs(x="Treatment", y= "Methylation level") +
  ggpubr::stat_cor(size = 3) +
  ggtitle("CTA methylation level in Decitabine-exposed MEFs") +
  scale_colour_discrete('Treatment Condition') +
  scale_shape_discrete('Replicate') 
rDNA_DEC_CTA_meth

library(forcats)
library(ggplot2)
rDNA_DEC_CTA_meth <- link_table_i_filtered %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% 
  filter(haplotype == 'CTA') %>% ggplot() +
  aes(x = Treatment, y = meth) +
  geom_point(aes(colour = Treatment, shape = Replicate),  size = 2) +
  theme_bw() +
  labs(x="Treatment", y= "Methylation level") +
  ggpubr::stat_cor(size = 3) +
  ggtitle("CTA methylation level in Decitabine-exposed MEFs") +
  scale_colour_discrete('Treatment Condition') +
  scale_shape_discrete('Replicate') 
rDNA_DEC_CTA_meth

library(forcats)
library(ggplot2)
rDNA_DEC_CCA_meth <- link_table_i_filtered %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% 
  filter(haplotype == 'CCA') %>% ggplot() +
  aes(x = Treatment, y = meth) +
  geom_point(aes(colour = Treatment, shape = Replicate),  size = 2) +
  theme_bw() +
  labs(x="Treatment", y= "Methylation level") +
  ggpubr::stat_cor(size = 3) +
  ggtitle("CCA methylation level in Decitabine-exposed MEFs") +
  scale_colour_discrete('Treatment Condition') +
  scale_shape_discrete('Replicate') 
rDNA_DEC_CCA_meth

write.csv(link_table_i_filtered, "MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/important/link_table_i_filtered.csv", row.names = FALSE)





#8.7.22 - blink nd

library(readr)
haplotype_specific_alleles <- read_csv("MEF_SD_DEC/RRBS/Data_Novogene/scripts/blink/haplotype_specific_alleles.csv")

haplotypes_mef <- haplotype_specific_alleles %>%
  group_by(pos_ref, allele, dominant) %>%
  mutate(num_hap=n(), haplotype=if_else(num_hap==1, haplotype, 'other')) %>%
  ungroup() %>%
  select(pos_ref, allele, haplotype) %>%
  distinct()


adjust_haplo <- function(haplotype)  {
  case_when(
    haplotype=='A_G_T_A' ~ 'ATA',
    haplotype=='A_G_T_G' ~ 'ATG',
    haplotype=='C_T_C_A' ~ 'CCA',
    haplotype=='C_T_T_A' ~ 'CTA',
    TRUE ~ haplotype)
}

library(dplyr)
link_table_nd_adj <- link_table_nd %>%
  mutate(adj_snp_pos=if_else(snp_pos<= 3008, snp_pos-3009, snp_pos-3008))
link_table_nd_adj

link_table_nd_adj = subset(link_table_nd_adj, select = -snp_pos )

link_table_nd_adj <- link_table_nd_adj %>%
  inner_join(haplotypes_mef %>%
               mutate(haplotype=adjust_haplo(haplotype)),
             by=c('adj_snp_pos'='pos_ref',
                  'snp_allele'='allele'))

link_nd_Treatment <- link_table_nd_adj %>%
  mutate(Treatment = case_when(
    endsWith(sample, "_0") ~ "CON_0",
    endsWith(sample, "_CON") ~ "CON",
    endsWith(sample, "_VEH") ~ "VEH",
    endsWith(sample, "_DEC") ~ "DEC",
    endsWith(sample, "_S") ~ "S_DEP",
    endsWith(sample, "_S_D") ~ "S_DEP_DEC"
  ))

link_nd_Treatment

link_nd_Treatment_replicate <- link_nd_Treatment %>%
  mutate(Replicate = case_when(
    startsWith(sample, "N1") ~ "1",
    startsWith(sample, "N2") ~ "2",
    startsWith(sample, "N3") ~ "3",
    startsWith(sample, "N4") ~ "4"
  ))

link_nd_Treatment_replicate



link_table_nd_filtered <- link_nd_Treatment_replicate %>%
  group_by(sample, adj_snp_pos, snp_allele, meth_pos, haplotype) %>%
  filter(n()>=10) %>%
  group_by(sample, adj_snp_pos, meth_pos, Treatment, Replicate) %>%
  filter(n()>=50, n_distinct(haplotype)>1) %>%
  group_by(sample, haplotype, Treatment, Replicate) %>%
  summarise(meth=sum(meth_state=='+')/n(),
            .groups='drop')

library(forcats)
library(ggplot2)
rDNA_DEC_ATA_meth <- link_table_nd_filtered %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% 
  filter(haplotype == 'ATA') %>% ggplot() +
  aes(x = Treatment, y = meth) +
  geom_jitter(aes(colour = Treatment, shape = Replicate),  size = 2) +
  theme_bw() +
  labs(x="Treatment", y= "Methylation level") +
  ggpubr::stat_cor(size = 3) +
  ggtitle("ATA methylation level in Decitabine-exposed MEFs") +
  scale_colour_discrete('Treatment Condition') +
  scale_shape_discrete('Replicate') 
rDNA_DEC_ATA_meth

library(forcats)
library(ggplot2)
rDNA_DEC_ATG_meth <- link_table_nd_filtered %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% 
  filter(haplotype == 'ATG') %>% ggplot() +
  aes(x = Treatment, y = meth) +
  geom_point(aes(colour = Treatment, shape = Replicate),  size = 2) +
  theme_bw() +
  labs(x="Treatment", y= "Methylation level") +
  ggpubr::stat_cor(size = 3) +
  ggtitle("ATG methylation level in Decitabine-exposed MEFs") +
  scale_colour_discrete('Treatment Condition') +
  scale_shape_discrete('Replicate') 
rDNA_DEC_ATG_meth

library(forcats)
library(ggplot2)
rDNA_DEC_CTA_meth <- link_table_nd_filtered %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% 
  filter(haplotype == 'CTA') %>% ggplot() +
  aes(x = Treatment, y = meth) +
  geom_point(aes(colour = Treatment, shape = Replicate),  size = 2) +
  theme_bw() +
  labs(x="Treatment", y= "Methylation level") +
  ggpubr::stat_cor(size = 3) +
  ggtitle("CTA methylation level in Decitabine-exposed MEFs") +
  scale_colour_discrete('Treatment Condition') +
  scale_shape_discrete('Replicate') 
rDNA_DEC_CTA_meth

library(forcats)
library(ggplot2)
rDNA_DEC_CTA_meth <- link_table_nd_filtered %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% 
  filter(haplotype == 'CTA') %>% ggplot() +
  aes(x = Treatment, y = meth) +
  geom_point(aes(colour = Treatment, shape = Replicate),  size = 2) +
  theme_bw() +
  labs(x="Treatment", y= "Methylation level") +
  ggpubr::stat_cor(size = 3) +
  ggtitle("CTA methylation level in Decitabine-exposed MEFs") +
  scale_colour_discrete('Treatment Condition') +
  scale_shape_discrete('Replicate') 
rDNA_DEC_CTA_meth

library(forcats)
library(ggplot2)
rDNA_DEC_CCA_meth <- link_table_nd_filtered %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% 
  filter(haplotype == 'CCA') %>% ggplot() +
  aes(x = Treatment, y = meth) +
  geom_point(aes(colour = Treatment, shape = Replicate),  size = 2) +
  theme_bw() +
  labs(x="Treatment", y= "Methylation level") +
  ggpubr::stat_cor(size = 3) +
  ggtitle("CCA methylation level in Decitabine-exposed MEFs") +
  scale_colour_discrete('Treatment Condition') +
  scale_shape_discrete('Replicate') 
rDNA_DEC_CCA_meth




# old below

#blink_nd
common_path2 = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/blink_nd/"
sub_dirs2 = list.files(common_path2)

snps_table_nd <- bind_rows(lapply(
  X=sub_dirs2,
  FUN=function(x){
    read_csv(
      paste0(common_path2, x, '/', x, '_snps.csv'),
      col_types=cols()) %>%
      mutate(sample=x, alt_prop_exp = 1 - ref_prop_exp,alt_prop_obs = 1 - ref_prop_obs)})) 
snps_table_nd


write.csv(snps_table_nd, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/snps_table_nd.csv", row.names = FALSE)

#blink_i
common_path_i_2 = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/blink/"
sub_dirs2_i = list.files(common_path_i_2)

snps_table_i <- bind_rows(lapply(
  X=sub_dirs2_i,
  FUN=function(x){
    read_csv(
      paste0(common_path_i_2, x, '/', x, '_snps.csv'),
      col_types=cols()) %>%
      mutate(sample=x, alt_prop_exp = 1 - ref_prop_exp,alt_prop_obs = 1 - ref_prop_obs)})) 
snps_table_i

write.csv(snps_table_i, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/snps_table_i.csv", row.names = FALSE)




fig_snps_nd <- ggplot(snps_table_nd) +
  geom_point() +
  aes(x = alt_prop_obs, y = alt_prop_exp) +
  labs(x="AAF in RRBS", y="AAF in WGS") +
  theme_minimal() +
  ggpubr::stat_cor()
fig_snps_nd



fig_snps_nd <- ggplot(snps_nd_Treatment_replicate) +
  geom_point() +
  aes(x = alt_prop_obs, y = alt_prop_exp) +
  labs(x="AAF in RRBS", y="AAF in WGS") +
  theme_minimal()
fig_snps_nd


fig_snps_nd_variables <- fig_snps_nd + facet_grid(vars(snps_nd_Treatment_replicate$Replicate), vars(snps_nd_Treatment_replicate$Treatment), scales = "free") +
  geom_point(aes(colour = Treatment), alpha = 0.5) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels =c('', 0.5, 1)) +
  scale_x_continuous(breaks =c(0, 0.5, 1), labels =c('', 0.5, 1))
fig_snps_nd_variables

library(ggpubr)
fig_snps_nd_variables + stat_cor(method = "pearson", size = 2)


snps_nd_Treatment <- snps_table_nd %>%
  mutate(Treatment = case_when(
    endsWith(sample, "_0") ~ "CON_0",
    endsWith(sample, "_CON") ~ "CON",
    endsWith(sample, "_VEH") ~ "VEH",
    endsWith(sample, "_DEC") ~ "DEC",
    endsWith(sample, "_S") ~ "S",
    endsWith(sample, "_S_D") ~ "S_D"
  ))

snps_nd_Treatment

snps_nd_Treatment_replicate <- snps_nd_Treatment %>%
  mutate(Replicate = case_when(
    startsWith(sample, "N1") ~ "1",
    startsWith(sample, "N2") ~ "2",
    startsWith(sample, "N3") ~ "3",
    startsWith(sample, "N4") ~ "4"
  ))

snps_nd_Treatment_replicate

link_nd_Treatment <- link_table_nd %>%
  mutate(Treatment = case_when(
    endsWith(sample, "_0") ~ "CON_0",
    endsWith(sample, "_CON") ~ "CON",
    endsWith(sample, "_VEH") ~ "VEH",
    endsWith(sample, "_DEC") ~ "DEC",
    endsWith(sample, "_S") ~ "S",
    endsWith(sample, "_S_D") ~ "S_D"
  ))

link_nd_Treatment

link_nd_Treatment_replicate <- link_nd_Treatment %>%
  mutate(Replicate = case_when(
    startsWith(sample, "N1") ~ "1",
    startsWith(sample, "N2") ~ "2",
    startsWith(sample, "N3") ~ "3",
    startsWith(sample, "N4") ~ "4"
  ))

link_nd_Treatment_replicate

#methylation level per CpG site associated with each snp
meth_level_CpG_nd <- link_nd_Treatment_replicate %>%
  group_by(sample, snp_pos, snp_allele, meth_pos, Treatment, Replicate) %>%
  summarise(meth = mean(meth_state=="+"))
meth_level_CpG_nd

write.csv(meth_level_CpG_nd, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/meth_level_CpG_nd.csv", row.names = FALSE)

#sample, snp position and allele. average across all CpG sites #20.7.21

meth_level_snp_nd <- meth_level_CpG_nd %>%
  group_by(sample, snp_pos, snp_allele, Treatment, Replicate) %>%
  summarise(meth = mean(meth))
meth_level_snp_nd

write.csv(meth_level_snp_nd, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/meth_level_snp_nd.csv", row.names = FALSE)

#3008 is -1, 3009 is 1. if position is below or equal to 3008, subtract 3009, if it is above 3008 subtract 3008.
library(dplyr)
meth_level_snp_nd_adj <- meth_level_snp_nd %>%
  mutate(adj_snp_pos=if_else(snp_pos<= 3008, snp_pos-3009, snp_pos-3008))
meth_level_snp_nd_adj


#plot
library(ggplot2)
meth_level_nd_graph <- ggplot(meth_level_snp_nd_adj) +
  aes(x = factor(adj_snp_pos), y = meth) +
  geom_point(aes(colour = snp_allele), size = 1, alpha = 0.5, position = position_dodge(width = 0.9)) +
  labs(x="SNP position in rDNA unit", y="Methylation level") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =c('', 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, size = 5)) +
  scale_color_discrete('Allele')
  theme_minimal()

meth_level_nd_1 <- meth_level_nd_graph + facet_grid(vars(meth_level_snp_nd_adj$Treatment))
meth_level_nd_1

#plot for ATA
library(ggplot2)
meth_level_nd_graph_ATA <- meth_level_snp_nd_adj %>% filter(adj_snp_pos %in% c(6007, 6832)) %>% 
  ggplot() + 
  aes(x = factor(adj_snp_pos), y = meth) +
  geom_point(aes(colour = snp_allele), size = 2, alpha = 0.5, position = position_dodge(width = 0.9)) +
  labs(x="SNP position in rDNA unit", y="Methylation level") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =c('', 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, size = 5)) +
  facet_grid(vars(Treatment)) +
  scale_color_discrete('Allele')
theme_minimal()
meth_level_nd_graph_ATA

#improved plot for ATA 
library(forcats)
library(ggplot2)
meth_level_nd_graph_ATA <- meth_level_snp_nd_adj %>% filter(adj_snp_pos %in% c(6007, 6832)) %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% 
  ggplot() + 
  aes(x = Treatment, y = meth) +
  geom_point(aes(colour = snp_allele), size = 2, alpha = 0.5, position = position_dodge(width = 0.2)) +
  labs(x="SNP position in rDNA unit", y="Methylation level") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =c('', 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_minimal() +
  facet_grid(vars(adj_snp_pos)) +
  scale_color_discrete('Allele')
theme_minimal()
meth_level_nd_graph_ATA


#plot for -104
library(forcats)
library(ggplot2)
meth_level_nd_graph_ATA_104 <- meth_level_snp_nd_adj %>% filter(adj_snp_pos %in% -104) %>%  mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>%
  ggplot() + 
  aes(x = Treatment, y = meth) +
  geom_point(aes(colour = snp_allele), size = 2, alpha = 0.5, position = position_dodge(width = 0.2)) +
  labs(x="Treatment condition", y="Methylation level") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =c('', 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_minimal() +
  facet_grid(vars(adj_snp_pos)) +
  scale_color_discrete('Allele')
theme_minimal()
meth_level_nd_graph_ATA_104



#ATA table 
meth_level_nd_6832_table <- meth_level_snp_nd_adj %>% filter(adj_snp_pos %in% c(6007, 6832))
meth_level_nd_6832_table




meth_level_snp_nd_adj %>% filter(adj_snp_pos %in% 6007) #no values given

library(dplyr)
meth_6832_VEH_DEC <- meth_level_nd_6832_table %>% filter(Treatment == 'DEC' | Treatment == 'VEH')
meth_6832_VEH_DEC

meth_level_snp_nd_adj %>% filter(adj_snp_pos %in% -104) %>% filter()


#ATG - 12736 - 28S

#plot for ATG
library(ggplot2)
meth_level_nd_graph_ATG <- meth_level_snp_nd_adj %>% filter(adj_snp_pos %in% c(12736)) %>% 
  ggplot() + 
  aes(x = factor(adj_snp_pos), y = meth) +
  geom_point(aes(colour = snp_allele), size = 2, alpha = 0.5, position = position_dodge(width = 0.9)) +
  labs(x="SNP position in rDNA unit", y="Methylation level") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =c('', 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, size = 5)) +
  facet_grid(vars(Treatment)) +
  scale_color_discrete('Allele')
theme_minimal()
meth_level_nd_graph_ATG

meth_level_snp_nd_adj %>% filter(adj_snp_pos %in% c(12736))

#CTA
library(ggplot2)
meth_level_nd_graph_CTA <- meth_level_snp_nd_adj %>% filter(adj_snp_pos %in% 1091)%>% 
  ggplot() + 
  aes(x = factor(adj_snp_pos), y = meth) +
  geom_point(aes(colour = snp_allele), size = 2, alpha = 0.5, position = position_dodge(width = 0.9)) +
  labs(x="SNP position in rDNA unit", y="Methylation level") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =c('', 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, size = 5)) +
  facet_grid(vars(Treatment)) +
  scale_color_discrete('Allele')
theme_minimal()
meth_level_nd_graph_CTA


#not helpful
library(tidyr)
meth_level_snp_nd_adj_pivot<- meth_level_snp_nd_adj %>% pivot_wider(c(sample, snp_pos, snp_allele, Replicate, adj_snp_pos), names_from=Treatment, values_from=meth)
meth_level_snp_nd_adj_pivot



#snps_table again
library(readr)
snps_table_nd <- read_csv("MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/snps_table_nd.csv")
snps_table_nd

snps_nd_Treatment <- snps_table_nd %>%
  mutate(Treatment = case_when(
    endsWith(sample, "_0") ~ "CON_0",
    endsWith(sample, "_CON") ~ "CON",
    endsWith(sample, "_VEH") ~ "VEH",
    endsWith(sample, "_DEC") ~ "DEC",
    endsWith(sample, "_S") ~ "S",
    endsWith(sample, "_S_D") ~ "S_D"
  ))

snps_nd_Treatment

snps_nd_Treatment_replicate <- snps_nd_Treatment %>%
  mutate(Replicate = case_when(
    startsWith(sample, "N1") ~ "1",
    startsWith(sample, "N2") ~ "2",
    startsWith(sample, "N3") ~ "3",
    startsWith(sample, "N4") ~ "4"
  ))

snps_nd_Treatment_replicate

write.csv(snps_nd_Treatment_replicate,"/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/snps_nd_Treatment_replicate.csv", row.names = FALSE)

#3008 is -1, 3009 is 1. if position is below or equal to 3008, subtract 3009, if it is above 3008 subtract 3008.
library(dplyr)
snps_nd_adj <- snps_nd_Treatment_replicate %>%
  mutate(adj_snp_pos=if_else(snp_pos<= 3008, snp_pos-3009, snp_pos-3008))
snps_nd_adj

write.csv(snps_nd_adj,"/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/snps_nd_adj.csv", row.names = FALSE)

#ATA
snps_nd_adj %>% filter(adj_snp_pos %in% 6832) %>% filter(Treatment == 'DEC'| Treatment == 'VEH')

snps_nd_adj %>% filter(adj_snp_pos %in% -104) %>% filter(Treatment == 'DEC'| Treatment == 'VEH')

snps_nd_adj %>% filter(adj_snp_pos %in% 6832) %>% filter(Treatment == 'CON'| Treatment == 'VEH')

snps_nd_adj %>% filter(adj_snp_pos %in% 6832) %>% filter(Treatment == 'CON'| Treatment == 'S')

#CCA
snps_nd_adj %>% filter(adj_snp_pos %in% 174) %>% filter(Treatment == 'DEC'| Treatment == 'VEH')

snps_nd_adj %>% filter(adj_snp_pos %in% 174) %>% filter(Treatment == 'CON'| Treatment == 'VEH')

#CTA
snps_nd_adj %>% filter(adj_snp_pos %in% 1091) %>% filter(Treatment == 'DEC'| Treatment == 'VEH')

snps_nd_adj %>% filter(adj_snp_pos %in% 1091) %>% filter(Treatment == 'CON'| Treatment == 'VEH')

#### back to meth level - statistical analysis

meth_level_snp_nd_adj

meth_level_snp_nd_adj %>%
  group_by(sample, Treatment, Replicate, adj_snp_pos) %>%
  filter(n_distinct(snp_allele)==2) %>%
  mutate(group=if_else(snp_allele==min(snp_allele), 'm', 'M')) %>%
  do(w = wilcox.test(meth ~ group, data=., exact=FALSE)) %>%
  summarise(sample, Treatment, Replicate, adj_snp_pos, pval=w$p.value, .groups='drop') %>%
  group_by(sample) %>%
  mutate(padj=p.adjust(pval, method = 'fdr')) %>%
  group_by(Treatment, adj_snp_pos) %>%
  summarise(num_sig=sum(padj<0.01), num_samples=n()) %>%
  filter(num_sig>0)


#Correct wilcoxon test - multiple CpG sites per position 
#3008 is -1, 3009 is 1. if position is below or equal to 3008, subtract 3009, if it is above 3008 subtract 3008.
library(dplyr)
meth_level_CpG_nd_adj <- meth_level_CpG_nd %>%
  mutate(adj_snp_pos=if_else(snp_pos<= 3008, snp_pos-3009, snp_pos-3008))
meth_level_CpG_nd_adj

meth_level_CpG_wilcoxon <- meth_level_CpG_nd_adj %>%
  group_by(sample, Treatment, Replicate, adj_snp_pos) %>%
  filter(n_distinct(snp_allele)==2) %>%
  mutate(group=if_else(snp_allele==min(snp_allele), 'm', 'M')) %>%
  do(w = wilcox.test(meth ~ group, data=., exact=FALSE)) %>%
  summarise(sample, Treatment, Replicate, adj_snp_pos, pval=w$p.value, .groups='drop') %>%
  group_by(sample) %>%
  mutate(padj=p.adjust(pval, method = 'fdr')) %>%
  group_by(Treatment, adj_snp_pos) %>%
  summarise(num_sig=sum(padj<0.01), num_samples=n()) %>%
  filter(num_sig>0)



meth_level_for_join <- meth_level_CpG_nd_adj %>% 
  group_by(Treatment, Replicate, adj_snp_pos, snp_allele) %>%
  summarise(meth = mean(meth))
meth_level_for_join

meth_level_wilcoxon_nd_forfig <- inner_join(meth_level_cpg_wilcoxon, meth_level_for_join, by = c('Treatment', 'adj_snp_pos'))
meth_level_wilcoxon_nd_forfig



library(ggplot2)
meth_level_wilcoxon_graph_nd <- ggplot(meth_level_wilcoxon_nd_forfig) +
  aes(x = factor(adj_snp_pos), y = meth) +
  geom_point(aes(colour = snp_allele), size = 2, alpha = 0.5, position = position_dodge(width = 0.9)) +
  labs(x="SNP position in rDNA unit", y="Methylation level") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =c('', 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, size = 5)) +
  facet_grid(rows = vars(Treatment)) +
  scale_color_discrete('Allele')
theme_minimal()



#ATA filter(adj_snp_pos %in% 6832) - A (C in other haplotypes)

library(ggplot2)

meth_level_wilcoxon_graph_nd  <- meth_level_wilcoxon_nd_forfig %>% filter(adj_snp_pos %in% c(-104, 6832)) %>% ggplot() +
  aes(x = factor(adj_snp_pos), y = meth) +
  geom_point(aes(colour = snp_allele), size = 2, alpha = 0.5, position = position_dodge(width = 0.9)) +
  labs(x="SNP position in rDNA unit", y="Methylation level") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =c('', 0.2, 0.4, 0.6, 0.8, 1)) +
  facet_grid(rows = vars(Treatment)) +
  scale_color_discrete('Allele') + 
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, size = 5))
meth_level_wilcoxon_graph_nd

#6832 only
meth_level_wilcoxon_graph_nd  <- meth_level_wilcoxon_nd_forfig %>% filter(adj_snp_pos %in% 6832) %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% ggplot() +
  aes(x = Treatment, y = meth) +
  geom_point(aes(colour = snp_allele), size = 2, alpha = 0.5, position = position_dodge(width = 0.2)) +
  labs(x="SNP position in rDNA unit", y="Methylation level") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =c('', 0.2, 0.4, 0.6, 0.8, 1)) +
  facet_grid(rows = vars(adj_snp_pos)) +
  scale_color_discrete('Allele') + 
  theme_minimal()
meth_level_wilcoxon_graph_nd

#104
library(forcats)
meth_level_wilcoxon_graph_nd  <- meth_level_wilcoxon_nd_forfig %>% filter(adj_snp_pos %in% -104) %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>%
  ggplot() +
  aes(x = Treatment, y = meth) +
  geom_point(aes(colour = snp_allele), size = 2, alpha = 0.5, position = position_dodge(width = 0.2)) +
  labs(x="SNP position in rDNA unit", y="Methylation level") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =c('', 0.2, 0.4, 0.6, 0.8, 1)) +
  facet_grid(rows = vars(adj_snp_pos)) +
  scale_color_discrete('Allele') + 
  theme_minimal()
meth_level_wilcoxon_graph_nd


#ATG filter(adj_snp_pos %in% 12736) - G (A in other haplotypes)

#CCA filter(adj_snp_pos %in% 174) - C (T in other haplotypes)

library(ggplot2)

meth_level_wilcoxon_graph_nd  <- meth_level_wilcoxon_nd_forfig %>% filter(adj_snp_pos %in% 174) %>% ggplot() +
  aes(x = factor(adj_snp_pos), y = meth) +
  geom_point(aes(colour = snp_allele), size = 2, alpha = 0.5, position = position_dodge(width = 0.9)) +
  labs(x="SNP position in rDNA unit", y="Methylation level") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =c('', 0.2, 0.4, 0.6, 0.8, 1)) +
  facet_grid(rows = vars(Treatment)) +
  scale_color_discrete('Allele') + 
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, size = 5))
meth_level_wilcoxon_graph_nd



#CTA filter(adj_snp_pos %in% 1091) - C (T in othe haplotypes)

library(ggplot2)

meth_level_wilcoxon_graph_nd  <- meth_level_wilcoxon_nd_forfig %>% filter(adj_snp_pos %in% 1091) %>% ggplot() +
  aes(x = factor(adj_snp_pos), y = meth) +
  geom_point(aes(colour = snp_allele), size = 2, alpha = 0.5, position = position_dodge(width = 0.9)) +
  labs(x="SNP position in rDNA unit", y="Methylation level") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =c('', 0.2, 0.4, 0.6, 0.8, 1)) +
  facet_grid(rows = vars(Treatment)) +
  scale_color_discrete('Allele') + 
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, size = 5))
meth_level_wilcoxon_graph_nd


#snps AAFs

ATA_AAF_table <- snps_nd_adj %>% filter(adj_snp_pos %in% 6832) 
ATA_AAF_table

write.csv(ATA_AAF_table, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/ATA_AAF_table.csv", row.names = FALSE)

library(ggplot2)
ATA_AAF <-ATA_AAF_table %>%  ggplot() +
  aes(x = Treatment, y= alt_prop_obs)+
  geom_point(aes(color = Replicate), alpha = 1,) +
  theme_bw()+
  labs(x = "Treatment condition", y = "AAF of 'A' at position 6832")
  ggtitle("AAF of rDNA haplotype: ATA")
ATA_AAF

ggplot(ATA_AAF_table, aes(x=Treatment , y=alt_prop_obs)) +
  geom_jitter(aes(color= Replicate), width=0.1, alpha=0.5) +
  ggpubr::stat_compare_means(aes(label = ..p.signif..),
                             ref.group = "CON", hide.ns=T,
                             bracket.size=1, label.y=0.2) +
  theme_bw() +
  labs(x="Treatment Condition", y= "AAF of 'A' at position 6832") +
  coord_cartesian(ylim=c(0, 0.2))+
  ggtitle("AAF of rDNA haplotype: ATA")


#CCA - C is the associated allele and alternative allele in table 

CCA_table <- snps_nd_adj %>% filter(adj_snp_pos %in% 174)
CCA_table

ggplot(CCA_table, aes(x=Treatment , y=alt_prop_obs)) +
  geom_jitter(aes(color= Replicate), width=0.1, alpha=0.5) +
  ggpubr::stat_compare_means(aes(label = ..p.signif..),
                             ref.group = "CON", hide.ns=T,
                             bracket.size=1, label.y=0.2) +
  theme_bw() +
  labs(x="Treatment Condition", y= "AAF of 'C' at position 174") +
  coord_cartesian(ylim=c(0, 0.2))+
  ggtitle("AAF of rDNA haplotype: CCA")


#CTA
CTA_AAF_table <- snps_nd_adj %>% filter(adj_snp_pos %in% 1091) 
CTA_AAF_table

ggplot(CTA_AAF_table, aes(x=Treatment , y=alt_prop_obs)) +
  geom_jitter(aes(color= Replicate), width=0.1, alpha=0.5) +
  ggpubr::stat_compare_means(aes(label = ..p.signif..),
                             ref.group = "CON", hide.ns=T,
                             bracket.size=1, label.y=0.7) +
  theme_bw() +
  labs(x="Treatment Condition", y= "AAF of 'C' at position 1091") +
  coord_cartesian(ylim=c(0, 0.7))+
  ggtitle("AAF of rDNA haplotype: CTA")

#ATG

ATG_AAF_table <- snps_nd_adj %>% filter(adj_snp_pos %in% 12736) 
ATG_AAF_table

###############################################################################
###Analysis: informed mode

library(readr)
snps_table_i <- read_csv("MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/snps_table_i.csv")
snps_table_i

snps_i_Treatment <- snps_table_i %>%
  mutate(Treatment = case_when(
    endsWith(sample, "_0") ~ "CON_0",
    endsWith(sample, "_CON") ~ "CON",
    endsWith(sample, "_VEH") ~ "VEH",
    endsWith(sample, "_DEC") ~ "DEC",
    endsWith(sample, "_S") ~ "S_DEP",
    endsWith(sample, "_S_D") ~ "S_DEP_DEC"
  ))

snps_i_Treatment

snps_i_Treatment_replicate <- snps_i_Treatment %>%
  mutate(Replicate = case_when(
    startsWith(sample, "N1") ~ "1",
    startsWith(sample, "N2") ~ "2",
    startsWith(sample, "N3") ~ "3",
    startsWith(sample, "N4") ~ "4"
  ))

snps_i_Treatment_replicate


fig_snps_i <- ggplot(snps_table_i) +
  geom_point() +
  aes(x = alt_prop_obs, y = alt_prop_exp) +
  labs(x="AAF in RRBS", y="AAF in WGS") +
  theme_minimal()
fig_snps_i


fig_snps_i <- ggplot(snps_i_Treatment_replicate) +
  geom_point() +
  aes(x = alt_prop_obs, y = alt_prop_exp) +
  labs(x="AAF in RRBS", y="AAF in WGS") +
  theme_minimal()
fig_snps_i


fig_snps_i_variables <- fig_snps_i + facet_grid(vars(snps_i_Treatment_replicate$Replicate), vars(snps_i_Treatment_replicate$Treatment), scales = "free") +
  geom_point(aes(colour = Treatment), alpha = 0.5) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels =c('', 0.5, 1)) +
  scale_x_continuous(breaks =c(0, 0.5, 1), labels =c('', 0.5, 1))
fig_snps_i_variables

library(ggpubr)
fig_snps_i_variables + stat_cor(method = "pearson", size = 2)

#3008 is -1, 3009 is 1. if position is below or equal to 3008, subtract 3009, if it is above 3008 subtract 3008.
library(dplyr)
snps_i_adj <- snps_i_Treatment_replicate %>%
  mutate(adj_snp_pos=if_else(snp_pos<= 3008, snp_pos-3009, snp_pos-3008))
snps_i_adj

write.csv(snps_i_adj,"/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/snps_i_adj.csv", row.names = FALSE)

#snps (i) AAFs

ATA_AAF_table_i <- snps_i_adj %>% filter(adj_snp_pos %in% 6832) 
ATA_AAF_table_i

write.csv(ATA_AAF_table_i, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/ATA_AAF_table_i.csv", row.names = FALSE)

library(forcats)
library(ggplot2)
ATA_AAF_i <-ATA_AAF_table_i %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% ggplot() +
  aes(x = Treatment, y= alt_prop_obs)+
  geom_point(aes(color = Replicate), alpha = 1,) +
  theme_bw()+
  labs(x = "Treatment condition", y = "Frequency of 'A' at position 6832") +
  theme(axis.text = element_text(size = 7)) +
ggtitle("AAF of rDNA haplotype: ATA (blink informed mode)")
ATA_AAF_i

ggplot(ATA_AAF_table, aes(x=Treatment , y=alt_prop_obs)) +
  geom_jitter(aes(color= Replicate), width=0.1, alpha=0.5) +
  ggpubr::stat_compare_means(aes(label = ..p.signif..),
                             ref.group = "CON", hide.ns=T,
                             bracket.size=1, label.y=0.2) +
  theme_bw() +
  labs(x="Treatment Condition", y= "AAF of 'A' at position 6832") +
  coord_cartesian(ylim=c(0, 0.2))+
  ggtitle("AAF of rDNA haplotype: ATA")


#CCA - C is the associated allele and alternative allele in table 

CCA_table_i <- snps_i_adj %>% filter(adj_snp_pos %in% 174)
CCA_table_i

library(forcats)
library(ggplot2)
CCA_i <-CCA_table_i %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% ggplot() +
  aes(x = Treatment, y= alt_prop_obs)+
  geom_point(aes(color = Replicate), alpha = 1,) +
  theme_bw()+
  labs(x = "Treatment condition", y = "Frequency of 'A' at position 6832") +
  theme(axis.text = element_text(size = 7)) +
  ggtitle("AAF of rDNA haplotype: ATA (blink informed mode)")
ATA_AAF_i

ggplot(CCA_table_i, aes(x=Treatment , y=alt_prop_obs)) +
  geom_jitter(aes(color= Replicate), width=0.1, alpha=0.5) +
  ggpubr::stat_compare_means(aes(label = ..p.signif..),
                             ref.group = "CON", hide.ns=T,
                             bracket.size=1, label.y=0.2) +
  theme_bw() +
  labs(x="Treatment Condition", y= "Frequency of 'C' at position 174") +
  coord_cartesian(ylim=c(0, 0.2))+
  ggtitle("AAF of rDNA haplotype: CCA (blink informed mode)")


#CTA
CTA_AAF_table_i <- snps_i_adj %>% filter(adj_snp_pos %in% 1091) 
CTA_AAF_table_i

ggplot(CTA_AAF_table_i, aes(x=Treatment , y=alt_prop_obs)) +
  geom_jitter(aes(color= Replicate), width=0.1, alpha=0.5) +
  ggpubr::stat_compare_means(aes(label = ..p.signif..),
                             ref.group = "CON", hide.ns=T,
                             bracket.size=1, label.y=0.7) +
  theme_bw() +
  labs(x="Treatment Condition", y= "AAF of 'C' at position 1091") +
  coord_cartesian(ylim=c(0, 0.7))+
  ggtitle("AAF of rDNA haplotype: CTA")

#ATG

ATG_AAF_table_i <- snps_i_adj %>% filter(adj_snp_pos %in% 12736) 
ATG_AAF_table_i

##link - i 

library(readr)
link_table_i <- read_csv("MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/link_table_i.csv")
link_table_i

link_i_Treatment <- link_table_i %>%
  mutate(Treatment = case_when(
    endsWith(sample, "_0") ~ "CON_0",
    endsWith(sample, "_CON") ~ "CON",
    endsWith(sample, "_VEH") ~ "VEH",
    endsWith(sample, "_DEC") ~ "DEC",
    endsWith(sample, "_S") ~ "S_DEP",
    endsWith(sample, "_S_D") ~ "S_DEP_DEC"
  ))

link_i_Treatment

link_i_Treatment_replicate <- link_i_Treatment %>%
  mutate(Replicate = case_when(
    startsWith(sample, "N1") ~ "1",
    startsWith(sample, "N2") ~ "2",
    startsWith(sample, "N3") ~ "3",
    startsWith(sample, "N4") ~ "4"
  ))

link_i_Treatment_replicate

#methylation level per CpG site associated with each snp
meth_level_CpG_i <- link_i_Treatment_replicate %>%
  group_by(sample, snp_pos, snp_allele, meth_pos, Treatment, Replicate) %>%
  summarise(meth = mean(meth_state=="+"))
meth_level_CpG_i

write.csv(meth_level_CpG_i, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/meth_level_CpG_i.csv", row.names = FALSE)

#sample, snp position and allele. average across all CpG sites #20.7.21

meth_level_snp_i <- meth_level_CpG_i %>%
  group_by(sample, snp_pos, snp_allele, Treatment, Replicate) %>%
  summarise(meth = mean(meth))
meth_level_snp_i

write.csv(meth_level_snp_i, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/meth_level_snp_i.csv", row.names = FALSE)

#3008 is -1, 3009 is 1. if position is below or equal to 3008, subtract 3009, if it is above 3008 subtract 3008.
library(dplyr)
meth_level_snp_i_adj <- meth_level_snp_i %>%
  mutate(adj_snp_pos=if_else(snp_pos<= 3008, snp_pos-3009, snp_pos-3008))
meth_level_snp_i_adj


#plot
library(forcats)
library(ggplot2)
meth_level_i_graph <- meth_level_snp_i_adj %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% ggplot() +
  aes(x = factor(adj_snp_pos), y = meth) +
  geom_point(aes(colour = snp_allele), size = 1, alpha = 0.5, position = position_dodge(width = 0.9)) +
  labs(x="SNP position in rDNA unit", y="Methylation level") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =c('', 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, size = 5)) +
  scale_color_discrete('Allele')
theme_minimal()

meth_level_i_1 <- meth_level_i_graph + facet_grid(vars(meth_level_snp_i_adj$Treatment))
meth_level_i_1

#plot for ATA
library(ggplot2)
meth_level_nd_graph_ATA <- meth_level_snp_nd_adj %>% filter(adj_snp_pos %in% c(6007, 6832)) %>% 
  ggplot() + 
  aes(x = factor(adj_snp_pos), y = meth) +
  geom_point(aes(colour = snp_allele), size = 2, alpha = 0.5, position = position_dodge(width = 0.9)) +
  labs(x="SNP position in rDNA unit", y="Methylation level") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =c('', 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, size = 5)) +
  facet_grid(vars(Treatment)) +
  scale_color_discrete('Allele')
theme_minimal()
meth_level_nd_graph_ATA

#improved plot for ATA 
library(forcats)
library(ggplot2)
meth_level_nd_graph_ATA <- meth_level_snp_nd_adj %>% filter(adj_snp_pos %in% c(6007, 6832)) %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% 
  ggplot() + 
  aes(x = Treatment, y = meth) +
  geom_point(aes(colour = snp_allele), size = 2, alpha = 0.5, position = position_dodge(width = 0.2)) +
  labs(x="SNP position in rDNA unit", y="Methylation level") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =c('', 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_minimal() +
  facet_grid(vars(adj_snp_pos)) +
  scale_color_discrete('Allele')
theme_minimal()
meth_level_nd_graph_ATA


#plot for -104
library(forcats)
library(ggplot2)
meth_level_nd_graph_ATA_104 <- meth_level_snp_nd_adj %>% filter(adj_snp_pos %in% -104) %>%  mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>%
  ggplot() + 
  aes(x = Treatment, y = meth) +
  geom_point(aes(colour = snp_allele), size = 2, alpha = 0.5, position = position_dodge(width = 0.2)) +
  labs(x="Treatment condition", y="Methylation level") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =c('', 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_minimal() +
  facet_grid(vars(adj_snp_pos)) +
  scale_color_discrete('Allele')
theme_minimal()
meth_level_nd_graph_ATA_104



#ATA table 
meth_level_nd_6832_table <- meth_level_snp_nd_adj %>% filter(adj_snp_pos %in% c(6007, 6832))
meth_level_nd_6832_table




meth_level_snp_nd_adj %>% filter(adj_snp_pos %in% 6007) #no values given

library(dplyr)
meth_6832_VEH_DEC <- meth_level_nd_6832_table %>% filter(Treatment == 'DEC' | Treatment == 'VEH')
meth_6832_VEH_DEC

meth_level_snp_nd_adj %>% filter(adj_snp_pos %in% -104) %>% filter()


#ATG - 12736 - 28S

#plot for ATG
library(ggplot2)
meth_level_nd_graph_ATG <- meth_level_snp_nd_adj %>% filter(adj_snp_pos %in% c(12736)) %>% 
  ggplot() + 
  aes(x = factor(adj_snp_pos), y = meth) +
  geom_point(aes(colour = snp_allele), size = 2, alpha = 0.5, position = position_dodge(width = 0.9)) +
  labs(x="SNP position in rDNA unit", y="Methylation level") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =c('', 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, size = 5)) +
  facet_grid(vars(Treatment)) +
  scale_color_discrete('Allele')
theme_minimal()
meth_level_nd_graph_ATG

meth_level_snp_nd_adj %>% filter(adj_snp_pos %in% c(12736))

#CTA
library(ggplot2)
meth_level_nd_graph_CTA <- meth_level_snp_nd_adj %>% filter(adj_snp_pos %in% 1091)%>% 
  ggplot() + 
  aes(x = factor(adj_snp_pos), y = meth) +
  geom_point(aes(colour = snp_allele), size = 2, alpha = 0.5, position = position_dodge(width = 0.9)) +
  labs(x="SNP position in rDNA unit", y="Methylation level") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =c('', 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, size = 5)) +
  facet_grid(vars(Treatment)) +
  scale_color_discrete('Allele')
theme_minimal()
meth_level_nd_graph_CTA



###rerun snps table analysis without informed mode so we change axis labels. 

library(readr)
snps_table_nd <- read_csv("MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/snps_table_nd.csv")
snps_table_nd


snps_nd_Treatment <- snps_table_nd %>%
  mutate(Treatment = case_when(
    endsWith(sample, "_0") ~ "CON_0",
    endsWith(sample, "_CON") ~ "CON",
    endsWith(sample, "_VEH") ~ "VEH",
    endsWith(sample, "_DEC") ~ "DEC",
    endsWith(sample, "_S") ~ "S_DEP",
    endsWith(sample, "_S_D") ~ "S_DEP_DEC"
  ))

snps_nd_Treatment

snps_nd_Treatment_replicate <- snps_nd_Treatment %>%
  mutate(Replicate = case_when(
    startsWith(sample, "N1") ~ "1",
    startsWith(sample, "N2") ~ "2",
    startsWith(sample, "N3") ~ "3",
    startsWith(sample, "N4") ~ "4"
  ))

snps_nd_Treatment_replicate



fig_snps_nd <- ggplot(snps_i_Treatment_replicate) +
  geom_point() +
  aes(x = alt_prop_obs, y = alt_prop_exp) +
  labs(x="AAF in RRBS", y="AAF in WGS") +
  theme_minimal()
fig_snps_i


fig_snps_i_variables <- fig_snps_i + facet_grid(vars(snps_i_Treatment_replicate$Replicate), vars(snps_i_Treatment_replicate$Treatment), scales = "free") +
  geom_point(aes(colour = Treatment), alpha = 0.5) +
  scale_y_continuous(breaks=c(0, 0.5, 1), labels =c('', 0.5, 1)) +
  scale_x_continuous(breaks =c(0, 0.5, 1), labels =c('', 0.5, 1))
fig_snps_i_variables

library(ggpubr)
fig_snps_i_variables + stat_cor(method = "pearson", size = 2)

#3008 is -1, 3009 is 1. if position is below or equal to 3008, subtract 3009, if it is above 3008 subtract 3008.
library(dplyr)
snps_i_adj <- snps_i_Treatment_replicate %>%
  mutate(adj_snp_pos=if_else(snp_pos<= 3008, snp_pos-3009, snp_pos-3008))
snps_i_adj

write.csv(snps_i_adj,"/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/snps_i_adj.csv", row.names = FALSE)

#snps () AAFs

ATA_AAF_table_i <- snps_i_adj %>% filter(adj_snp_pos %in% 6832) 
ATA_AAF_table_i

write.csv(ATA_AAF_table_i, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/blink/ATA_AAF_table_i.csv", row.names = FALSE)

library(ggplot2)
ATA_AAF_i <-ATA_AAF_table_i %>%  ggplot() +
  aes(x = Treatment, y= alt_prop_obs)+
  geom_point(aes(color = Replicate), alpha = 1,) +
  theme_bw()+
  labs(x = "Treatment condition", y = "Frequency of 'A' at position 6832") +
  theme(axis.text = element_text(size = 7)) +
  ggtitle("AAF of rDNA haplotype: ATA (blink informed mode)")
ATA_AAF_i

ggplot(ATA_AAF_table, aes(x=Treatment , y=alt_prop_obs)) +
  geom_jitter(aes(color= Replicate), width=0.1, alpha=0.5) +
  ggpubr::stat_compare_means(aes(label = ..p.signif..),
                             ref.group = "CON", hide.ns=T,
                             bracket.size=1, label.y=0.2) +
  theme_bw() +
  labs(x="Treatment Condition", y= "AAF of 'A' at position 6832") +
  coord_cartesian(ylim=c(0, 0.2))+
  ggtitle("AAF of rDNA haplotype: ATA")


#CCA - C is the associated allele and alternative allele in table 

CCA_table <- snps_nd_adj %>% filter(adj_snp_pos %in% 174)
CCA_table

ggplot(CCA_table, aes(x=Treatment , y=alt_prop_obs)) +
  geom_jitter(aes(color= Replicate), width=0.1, alpha=0.5) +
  ggpubr::stat_compare_means(aes(label = ..p.signif..),
                             ref.group = "CON", hide.ns=T,
                             bracket.size=1, label.y=0.2) +
  theme_bw() +
  labs(x="Treatment Condition", y= "AAF of 'C' at position 174") +
  coord_cartesian(ylim=c(0, 0.2))+
  ggtitle("AAF of rDNA haplotype: CCA")


#CTA
CTA_AAF_table <- snps_nd_adj %>% filter(adj_snp_pos %in% 1091) 
CTA_AAF_table

ggplot(CTA_AAF_table, aes(x=Treatment , y=alt_prop_obs)) +
  geom_jitter(aes(color= Replicate), width=0.1, alpha=0.5) +
  ggpubr::stat_compare_means(aes(label = ..p.signif..),
                             ref.group = "CON", hide.ns=T,
                             bracket.size=1, label.y=0.7) +
  theme_bw() +
  labs(x="Treatment Condition", y= "AAF of 'C' at position 1091") +
  coord_cartesian(ylim=c(0, 0.7))+
  ggtitle("AAF of rDNA haplotype: CTA")

#ATG

ATG_AAF_table <- snps_nd_adj %>% filter(adj_snp_pos %in% 12736) 
ATG_AAF_table

