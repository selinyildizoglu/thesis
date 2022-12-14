.libPaths('/data/Blizard-Rakyan/Fran/R/packages_4.1.1/')

if(any(grepl("package:methylKit", search()))) 
  detach("package:methylKit") else message("methylKit not loaded")

detach("package:biomaRt")



library(readr)
reads_per_gene_mRNAseq <- read_csv("MEF_SD_DEC/mRNAseq/data/star/reads_per_gene_mRNAseq.csv")
reads_per_gene_mRNAseq




if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

library(edgeR)

library(dplyr)
library(stringr)
mrna_seq_treatment <- reads_per_gene_mRNAseq %>%
  mutate(Treatment = case_when(
    str_detect(sample, "N20") ~ "CON_0",
    str_detect(sample, "N2C") ~ "CON",
    str_detect(sample, "N2V") ~ "VEH",
    str_detect(sample, "N2D") ~ "DEC",
    str_detect(sample, "N2S") ~ "S_DEP",
    str_detect(sample, "N40") ~ "CON_0",
    str_detect(sample, "N4C") ~ "CON",
    str_detect(sample, "N4V") ~ "VEH",
    str_detect(sample, "N4D") ~ "DEC",
    str_detect(sample, "N4S") ~ "S_DEP",
    ))
mrna_seq_treatment

mrna_seq_treatment_replicate <- mrna_seq_treatment %>%
  mutate(Replicate = case_when(
    startsWith(sample, "N2") ~ "2",
    startsWith(sample, "N4") ~ "4"
  ))


design_table <- mrna_seq_treatment_replicate %>% select(sample, Treatment, Replicate)

design_table <- distinct(design_table)


library(tidyverse)

##Fran's code

#sample_info <- read_csv(file.path(fd, 'sample_ids_tidy.csv')) %>%
  #arrange(sample)
#tissues <- as.factor(sample_info$tissue)
#design <- model.matrix(~tissues)

sample_info <- design_table %>%
  arrange(sample)
Treatment <- as.factor(sample_info$Treatment)
Replicate <- as.factor(design_table$Replicate)
design <- model.matrix(~Treatment + ~Replicate)

reads_per_gene_mRNAseq

#read_counts <- read_csv(file.path(fd, 'read_counts.csv')) %>%
  #select(sample, gene_id, count=second_strand) %>%
  #spread(sample, count) %>%
  #column_to_rownames('gene_id')

read_counts <- reads_per_gene_mRNAseq %>%
  select(sample, gene_ID, count=second_read_strand) %>%
  spread(sample, count) %>%
  column_to_rownames('gene_ID')

dge_list <- calcNormFactors(DGEList(
  counts = as.matrix(read_counts),
  group = Treatment))
dge_est <- estimateDisp(dge_list, design)
fit <- glmQLFit(dge_est, design)
qlf <- glmQLFTest(fit, coef=2)
topTags(qlf)

#compare CON and S_DEP only


design_table_CON_S_DEP <- design_table %>% filter(Treatment %in% c('CON', 'S_DEP'))

sample_info_CON_S_DEP <- design_table_CON_S_DEP %>%
  arrange(sample)
Treatment <- as.factor(sample_info_CON_S_DEP$Treatment)
Replicate <- as.factor(sample_info_CON_S_DEP$Replicate)
design_CON_S_DEP <- model.matrix(~Treatment + ~Replicate)

reads_per_gene_mRNAseq_CON_S <- mrna_seq_treatment_replicate %>% filter(Treatment %in% c('CON', 'S_DEP')) %>% select(gene_ID, second_read_strand, sample)


read_counts_CON_S <- reads_per_gene_mRNAseq_CON_S  %>%
  select(sample, gene_ID, count=second_read_strand) %>%
  spread(sample, count) %>%
  column_to_rownames('gene_ID')

dge_list <- calcNormFactors(DGEList(
  counts = as.matrix(read_counts_CON_S),
  group = Treatment))
dge_est <- estimateDisp(dge_list, design_CON_S_DEP)
fit <- glmQLFit(dge_est, design_CON_S_DEP)
qlf_CON_S <- glmQLFTest(fit, coef=2)
topTags(qlf_CON_S)

toptags_CON_S_table <- topTags(qlf_CON_S) 

write.csv(toptags_CON_S_table, "MEF_SD_DEC/mRNAseq/data/csvs/toptags_CON_S_DEP.csv")

write.csv(qlf_CON_S, "MEF_SD_DEC/mRNAseq/data/csvs/qlf_CON_S.csv")

##GO analysis for CON_S using most significant - geenerating table, copyiing and pasting names and using the online tool 

CON_S_qlf_gene_ID <- topTags(qlf_CON_S, 1e9)[[1]] %>%
  rownames_to_column('gene_id')  %>%
  as_tibble()

#whole table 
CON_S_qlf_gene_ID_ <- separate(CON_S_qlf_gene_ID, col=gene_id, into=c('gene_id'), extra='drop')

#write.csv
write.csv(CON_S_qlf_gene_ID_, "MEF_SD_DEC/mRNAseq/data/csvs/CON_S_qlf_gene_ID_.csv")

#FDR < 0.01

CON_S_qlf_gene_ID_sig <- CON_S_qlf_gene_ID %>% filter(FDR < 0.01)

CON_S_qlf_gene_id_sig_ <- separate(CON_S_qlf_gene_ID_sig, col=gene_id, into=c('gene_id'), extra='drop')

gene_id_sig <- CON_S_qlf_gene_id_sig_$gene_id


write.csv(CON_S_qlf_gene_id_sig_, "MEF_SD_DEC/mRNAseq/data/csvs/CON_S_sig.csv", row.names = FALSE)


BiocManager::install("biomaRt")


##make upregulated and ddownregulated tables and then find top 10, make two lists, then run through script below 


CON_S_qlf_gene_ID #table to use

CON_S_qlf_negative <- CON_S_qlf_gene_ID %>% filter(logFC < 0) # downregulation

CON_S_qlf_pos <- CON_S_qlf_gene_ID %>% filter(logFC > 0) #upregulation


CON_S_qlf_negative <- CON_S_qlf_negative %>% arrange(FDR)

CON_S_qlf_pos <- CON_S_qlf_pos %>% arrange(FDR)

CON_S_qlf_negative_ <- separate(CON_S_qlf_negative, col=gene_id, into=c('gene_id'), extra='drop')

CON_S_qlf_pos_ <- separate(CON_S_qlf_pos, col=gene_id, into=c('gene_id'), extra='drop')

CON_S_qlf_pos %>% head(10)

CON_S_qlf_pos_top10 <- CON_S_qlf_pos %>% head(10)

CON_S_qlf_pos_top10 <- separate(CON_S_qlf_pos_top10, col=gene_id, into=c('gene_id'), extra='drop')

CON_S_qlf_negative_top10 <- CON_S_qlf_negative %>% head(10)

CON_S_qlf_neg_top10 <- separate(CON_S_qlf_negative_top10, col=gene_id, into=c('gene_id'), extra='drop')

#make gene lists

gene_list_neg_top10_CON_S <- CON_S_qlf_neg_top10$gene_id

gene_list_pos_top10_CON_S <- CON_S_qlf_pos_top10$gene_id

#run for CON_S neg 


ensembl_mm10 <- biomaRt::useEnsembl(biomart = "ensembl",
                                    dataset = 'mmusculus_gene_ensembl')

go_results_neg_CON_S <- biomaRt::getBM(
  attributes=c("ensembl_gene_id", "external_gene_name",
               "go_id", "name_1006", "namespace_1003"),
  filters="ensembl_gene_id",
  values=gene_list_neg_top10_CON_S,
  mart=ensembl_mm10) %>%
  as_tibble()


ensembl_mm10 <- biomaRt::useEnsembl(biomart = "ensembl",
                                    dataset = 'mmusculus_gene_ensembl', version =106)

go_results_pos_CON_0_DEC <- biomaRt::getBM(
  attributes=c("ensembl_gene_id", "external_gene_name",
               "go_id", "name_1006", "namespace_1003"),
  filters="ensembl_gene_id",
  values=gene_list_pos_con_0_dec_10,
  mart=ensembl_mm10) %>%
  as_tibble()




#plot


go_results_neg_CON_S_ %>% filter(!name_1006 %in% c('regulation of bone mineralization', 'female pregnancy', 'angiogenesis')) %>% group_by(name_1006) %>% filter(name_1006!="")  %>% filter(n()>1)  %>% ggplot() + geom_bar(aes(y=forcats::fct_infreq(name_1006))) +
  facet_grid(rows=vars(namespace_1003), scales='free_y') +
  labs(x='Number of genes', y='GO Term') +
  theme_bw() +
  theme(panel.grid=element_blank(), aspect.ratio = 1/2, plot.title=element_text(size=10)) +
  ggtitle("Genes that show downregulated expression in Serum deprivation 
          relative to Control Day 5")

  
go_results_pos_CON_S_ <- go_results_pos_CON_S %>% mutate(namespace_1003=gsub("_"," ", namespace_1003))

go_results_pos_CON_S_ %>% group_by(name_1006) %>% filter(name_1006!="")  %>% filter(n()>1)  %>% ggplot() + geom_bar(aes(y=forcats::fct_infreq(name_1006))) +
  facet_grid(rows=vars(namespace_1003), scales='free_y') +
  labs(x='Number of genes', y='GO Term') +
  theme_bw() +
  theme(panel.grid=element_blank(), aspect.ratio = 1/2, plot.title=element_text(size=10)) +
  ggtitle("Genes that show upregulated expression in Serum deprivation 
          relative to Control Day 5")


library(readr)
go_results_neg_CON_S <- read_csv("MEF_SD_DEC/mRNAseq/data/csvs/go_results_neg_CON_S.csv")
go_results_neg_CON_S

go_results_neg_CON_S_ <- go_results_neg_CON_S %>% mutate(namespace_1003=gsub("_"," ", namespace_1003))

go_results_neg_CON_DEC_ %>% filter(!name_1006 %in% c('neuron differentiation', 'nervous system development', 'in utero embryonic development', 'heart development')) %>% group_by(name_1006) %>% filter(name_1006!="")  %>% filter(n()>1)  %>% ggplot() + geom_bar(aes(y=forcats::fct_infreq(name_1006))) +
  facet_grid(rows=vars(namespace_1003), scales='free_y') +
  labs(x='Number of genes', y='GO Term') +
  theme_bw() +
  theme(panel.grid=element_blank(), aspect.ratio = 1/2, plot.title=element_text(size=10)) + 
  ggtitle("Genes that show downregulated expression 
          in Decitabine-exposed MEFs")


####ccode
go_results_pos_CON_DEC_ <- go_results_pos_CON_DEC %>% mutate(namespace_1003=gsub("_"," ", namespace_1003))


go_results_pos_CON_DEC_ %>% group_by(name_1006) %>% filter(!name_1006 %in% c('calcium-dependent cell-cell adhesion via plasma membrane cell adhesion molecules', 'in utero embryonic development', 'positive regulation of angiogenesis', 'spermatogenesis')) %>% filter(name_1006!="")  %>% filter(n()>1)  %>% ggplot() + geom_bar(aes(y=forcats::fct_infreq(name_1006))) +
  facet_grid(rows=vars(namespace_1003), scales='free_y') +
  labs(x='Number of genes', y='GO Term') +
  scale_x_continuous(limits=c(0, 8)) +
  theme_bw() +
  theme(panel.grid=element_blank(), aspect.ratio = 1/2, plot.title=element_text(size=10)) +
  ggtitle("Genes that show upregulated expression 
          in Decitabine-exposed MEFs")
####code



#write csvs

write.csv(go_results_neg_CON_S, "MEF_SD_DEC/mRNAseq/data/csvs/go_results_neg_CON_S.csv", row.names = FALSE)

write.csv(go_results_pos_CON_S, "MEF_SD_DEC/mRNAseq/data/csvs/go_results_pos_CON_S.csv", row.names = FALSE)

write.csv(CON_S_qlf_pos_top10, "MEF_SD_DEC/mRNAseq/data/csvs/CON_S_qlf_pos_top10.csv", row.names = FALSE)

write.csv(CON_S_qlf_neg_top10, "MEF_SD_DEC/mRNAseq/data/csvs/CON_S_qlf_neg_top10.csv", row.names = FALSE)

CON_S_qlf_negative_ <- separate(CON_S_qlf_negative, col=gene_id, into=c('gene_id'), extra='drop')

CON_S_qlf_pos_ <- separate(CON_S_qlf_pos, col=gene_id, into=c('gene_id'), extra='drop')


write.csv(CON_S_qlf_negative_, "MEF_SD_DEC/mRNAseq/data/csvs/CON_S_qlf_negative_.csv", row.names = FALSE)

write.csv(CON_S_qlf_pos_, "MEF_SD_DEC/mRNAseq/data/csvs/CON_S_qlf_pos_.csv", row.names = FALSE)


#run for CON_S pos

ensembl_mm10 <- biomaRt::useEnsembl(biomart = "ensembl",
                                    dataset = 'mmusculus_gene_ensembl')

go_results_pos_CON_S <- biomaRt::getBM(
  attributes=c("ensembl_gene_id", "external_gene_name",
               "go_id", "name_1006", "namespace_1003"),
  filters="ensembl_gene_id",
  values=gene_list_pos_top10_CON_S,
  mart=ensembl_mm10) %>%
  as_tibble()



go_results_pos_CON_S %>% filter(ensembl_gene_id == 'ENSMUSG00000027750')

go_results_pos_CON_S %>% filter(ensembl_gene_id == 'ENSMUSG00000036390')

#S_DEP and CON_0


design_table_CON_0_S_DEP <- design_table %>% filter(Treatment %in% c('CON_0', 'S_DEP'))

design_table_CON_0_S_DEP 

sample_info_CON_0_S_DEP <- design_table_CON_0_S_DEP %>%
  arrange(sample)
Treatment <- as.factor(sample_info_CON_0_S_DEP$Treatment)
Replicate <- as.factor(sample_info_CON_0_S_DEP$Replicate)
design_CON_0_S_DEP <- model.matrix(~Treatment + ~Replicate)

reads_per_gene_mRNAseq_CON_0_S <- mrna_seq_treatment_replicate %>% filter(Treatment %in% c('CON_0', 'S_DEP')) %>% select(gene_ID, second_read_strand, sample)

read_counts_CON_0_S <- reads_per_gene_mRNAseq_CON_0_S  %>%
  select(sample, gene_ID, count=second_read_strand) %>%
  spread(sample, count) %>%
  column_to_rownames('gene_ID')

dge_list <- calcNormFactors(DGEList(
  counts = as.matrix(read_counts_CON_0_S),
  group = Treatment))
dge_est <- estimateDisp(dge_list, design_CON_0_S_DEP)
fit <- glmQLFit(dge_est, design_CON_0_S_DEP)
qlf_CON_0_S <- glmQLFTest(fit, coef=2)
topTags(qlf_CON_0_S)

#add gene_id to column
CON_0_S_qlf_gene_ID <- topTags(qlf_CON_0_S, 1e9)[[1]] %>%
  rownames_to_column('gene_id')  %>%
  as_tibble()

Con_0_S_top10 <- topTags(qlf_CON_0_S)

#remove version to enable GO recognition - make top 10 plot

CON_0_S_qlf_gene_ID_ <- separate(CON_0_S_qlf_gene_ID, col=gene_id, into=c('gene_id'), extra='drop')

CON_0_S_top10 <- CON_0_S_qlf_gene_ID_ %>% head(10)


write.csv(CON_0_S_qlf_gene_ID_, "MEF_SD_DEC/mRNAseq/data/csvs/CON_0_S_qlf_gene_ID_.csv", row.names = FALSE)

write.csv(CON_0_S_top10, "MEF_SD_DEC/mRNAseq/data/csvs/CON_0_S_top10.csv", row.names = FALSE)

#FDR < 0.01

CON_0_S_FDR <- CON_0_S_qlf_gene_ID_ %>% filter(FDR < 0.01)

write.csv(CON_0_S_FDR, "MEF_SD_DEC/mRNAseq/data/csvs/CON_0_S_sig.csv", row.names = FALSE)

#upreg and downreg tables

CON_0_S_pos <- CON_0_S_qlf_gene_ID_ %>% filter(logFC > 0)

CON_0_S_neg <- CON_0_S_qlf_gene_ID_ %>% filter(logFC < 0)

write.csv(CON_0_S_pos, "MEF_SD_DEC/mRNAseq/data/csvs/CON_0_S_pos.csv", row.names = FALSE)

write.csv(CON_0_S_neg, "MEF_SD_DEC/mRNAseq/data/csvs/CON_0_S_top10.csv", row.names = FALSE)

#arrange by FDR

CON_0_S_pos <- CON_0_S_pos %>% arrange(FDR)

CON_0_S_neg <- CON_0_S_neg  %>% arrange(FDR)

#top 10 pos and neg

CON_0_S_pos_10 <- CON_0_S_pos %>% head(10)

CON_0_S_neg_10 <- CON_0_S_neg %>% head(10)



##GO

#gene lists

gene_list_neg_con_0_S_10 <- CON_0_S_neg_10$gene_id

gene_list_pos_con_0_S_10 <- CON_0_S_pos_10$gene_id

#run for CON_0_S neg

ensembl_mm10 <- biomaRt::useEnsembl(biomart = "ensembl",
                                    dataset = 'mmusculus_gene_ensembl')

go_results_neg_CON_0_S <- biomaRt::getBM(
  attributes=c("ensembl_gene_id", "external_gene_name",
               "go_id", "name_1006", "namespace_1003"),
  filters="ensembl_gene_id",
  values=gene_list_neg_con_0_S_10,
  mart=ensembl_mm10) %>%
  as_tibble()

#run for Con_0_S pos

ensembl_mm10 <- biomaRt::useEnsembl(biomart = "ensembl",
                                    dataset = 'mmusculus_gene_ensembl')

go_results_pos_CON_0_S <- biomaRt::getBM(
  attributes=c("ensembl_gene_id", "external_gene_name",
               "go_id", "name_1006", "namespace_1003"),
  filters="ensembl_gene_id",
  values=gene_list_pos_con_0_S_10,
  mart=ensembl_mm10) %>%
  as_tibble()


#plot


go_results_neg_CON_0_S %>% group_by(name_1006) %>% filter(name_1006!="")  %>% filter(n()>1)  %>% ggplot() + geom_bar(aes(y=forcats::fct_infreq(name_1006))) +
  facet_grid(rows=vars(namespace_1003), scales='free_y') +
  labs(x='Number of genes', y='GO Term') +
  theme_bw() +
  ggtitle("Genes that show downregulated expression in Serum deprivation relative to Control Day 0")




go_results_pos_CON_0_S %>% group_by(name_1006) %>% filter(name_1006!="")  %>% filter(n()>1)  %>% ggplot() + geom_bar(aes(y=forcats::fct_infreq(name_1006))) +
  facet_grid(rows=vars(namespace_1003), scales='free_y') +
  labs(x='Number of genes', y='GO Term') +
  theme_bw() +
  ggtitle("Genes that show upregulated expression in Serum deprivation relative to Control Day 0")


#write csvs

write.csv(go_results_pos_CON_0_S, "MEF_SD_DEC/mRNAseq/data/csvs/go_results_pos_CON_0_S.csv", row.names = FALSE)

write.csv(go_results_neg_CON_0_S, "MEF_SD_DEC/mRNAseq/data/csvs/go_results_neg_CON_0_S.csv", row.names = FALSE)

###con_0 and con

library(forcats)

design_table_CON_CON_0 <- design_table %>% filter(Treatment %in% c('CON_0', 'CON')) %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0'))

design_table_CON_CON_0_old <- design_table %>% filter(Treatment %in% c('CON', 'CON_0'))


sample_info_CON_CON_0_old <- design_table_CON_CON_0_old %>%
  arrange(sample)
Treatment <- as.factor(sample_info_CON_CON_0_old$Treatment)
Replicate <- as.factor(sample_info_CON_CON_0_old$Replicate)
design_CON_CON_0_old <- model.matrix(~Treatment + ~Replicate)



sample_info_CON_CON_0 <- design_table_CON_CON_0 %>%
  arrange(sample)
Treatment <- as.factor(sample_info_CON_CON_0$Treatment)
Replicate <- as.factor(sample_info_CON_CON_0$Replicate)
design_CON_CON_0 <- model.matrix(~Treatment + ~Replicate)

reads_per_gene_mRNAseq_CON_CON_0 <- mrna_seq_treatment_replicate %>% filter(Treatment %in% c('CON', 'CON_0')) %>% select(gene_ID, second_read_strand, sample)



read_counts_CON_CON_0 <- reads_per_gene_mRNAseq_CON_CON_0 %>%
  select(sample, gene_ID, count=second_read_strand) %>%
  spread(sample, count) %>%
  column_to_rownames('gene_ID')

dge_list <- calcNormFactors(DGEList(
  counts = as.matrix(read_counts_CON_CON_0),
  group = Treatment))
dge_est <- estimateDisp(dge_list, design_CON_CON_0 )
fit <- glmQLFit(dge_est, design_CON_CON_0 )
qlf_CON_CON_0 <- glmQLFTest(fit, coef=2)
topTags(qlf_CON_CON_0)


#add gene_id to column
CON_CON_0_qlf_gene_ID <- topTags(qlf_CON_CON_0, 1e9)[[1]] %>%
  rownames_to_column('gene_id')  %>%
  as_tibble()



#remove version to enable GO recognition - make top 10 plot

CON_CON_0_qlf_gene_ID_2 <- separate(CON_CON_0_qlf_gene_ID, col=gene_id, into=c('gene_id'), extra='drop')

CON_CON_0_top10_2 <- CON_CON_0_qlf_gene_ID_2 %>% head(10)


write.csv(CON_CON_0_qlf_gene_ID_2, "MEF_SD_DEC/mRNAseq/data/csvs/CON_CON_0_qlf_gene_ID_2.csv", row.names = FALSE)

write.csv(CON_CON_0_top10_2, "MEF_SD_DEC/mRNAseq/data/csvs/CON_CON_0_top10_2.csv", row.names = FALSE)

#FDR < 0.01

CON_CON_0_FDR_2 <- CON_CON_0_qlf_gene_ID_2 %>% filter(FDR < 0.01)

write.csv(CON_CON_0_FDR_2, "MEF_SD_DEC/mRNAseq/data/csvs/CON_CON_0_sig_2.csv", row.names = FALSE)

#upreg and downreg tables

CON_CON_0_pos_2 <- CON_CON_0_qlf_gene_ID_2 %>% filter(logFC > 0)

CON_CON_0_neg_2 <- CON_CON_0_qlf_gene_ID_2 %>% filter(logFC < 0)

write.csv(CON_CON_0_pos_2, "MEF_SD_DEC/mRNAseq/data/csvs/CON_CON_0_pos_2.csv", row.names = FALSE)

write.csv(CON_CON_0_neg_2, "MEF_SD_DEC/mRNAseq/data/csvs/CON_CON_0_neg_2.csv", row.names = FALSE)

#arrange by FDR

CON_CON_0_pos_2 <- CON_CON_0_pos_2 %>% arrange(FDR)

CON_CON_0_neg_2 <- CON_CON_0_neg_2 %>% arrange(FDR)

#top 10 pos and neg

CON_CON_0_pos_10_2 <- CON_CON_0_pos_2 %>% head(10)

CON_CON_0_neg_10_2 <- CON_CON_0_neg_2 %>% head(10)



##GO

#gene lists

gene_list_neg_con_con_0_10_2 <- CON_CON_0_neg_10_2$gene_id

gene_list_pos_con_con_0_10_2 <- CON_CON_0_pos_10_2$gene_id

#run for CON_CON_0 neg

ensembl_mm10 <- biomaRt::useEnsembl(biomart = "ensembl",
                                    dataset ='mmusculus_gene_ensembl', version =106)

go_results_neg_CON_CON_0_2 <- biomaRt::getBM(
  attributes=c("ensembl_gene_id", "external_gene_name",
               "go_id", "name_1006", "namespace_1003"),
  filters="ensembl_gene_id",
  values=gene_list_neg_con_con_0_10_2,
  mart=ensembl_mm10) %>%
  as_tibble()



#run for Con_con_0 pos

ensembl_mm10 <- biomaRt::useEnsembl(biomart = "ensembl",
                                    dataset = 'mmusculus_gene_ensembl')

go_results_pos_CON_CON_0 <- biomaRt::getBM(
  attributes=c("ensembl_gene_id", "external_gene_name",
               "go_id", "name_1006", "namespace_1003"),
  filters="ensembl_gene_id",
  values=gene_list_pos_con_con_0_10,
  mart=ensembl_mm10) %>%
  as_tibble()


#plot


go_results_neg_CON_CON_0 %>% group_by(name_1006) %>% filter(name_1006!="")  %>% filter(n()>1)  %>% ggplot() + geom_bar(aes(y=forcats::fct_infreq(name_1006))) +
  facet_grid(rows=vars(namespace_1003), scales='free_y') +
  labs(x='Number of genes', y='GO Term') +
  theme_bw() +
  ggtitle("Genes that show downregulated expression in Control Day 0 relative to Control Day 5")


go_results_pos_CON_CON_0 %>% group_by(name_1006) %>% filter(name_1006!="") %>% filter(n()>1) %>% ggplot() + geom_bar(aes(y=forcats::fct_infreq(name_1006))) +
  facet_grid(rows=vars(namespace_1003), scales='free_y') +
  labs(x='Number of genes', y='GO Term') +
  theme_bw() +
  ggtitle("Genes that show upregulated expression in Control day 0 relative to Control Day 5")


#write csvs

write.csv(go_results_neg_CON_CON_0, "MEF_SD_DEC/mRNAseq/data/csvs/go_results_neg_CON_CON_0.csv", row.names = FALSE)

write.csv(go_results_pos_CON_CON_0, "MEF_SD_DEC/mRNAseq/data/csvs/go_results_pos_CON_CON_0.csv", row.names = FALSE)






