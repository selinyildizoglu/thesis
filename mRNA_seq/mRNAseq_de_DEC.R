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

#DEC and CON


design_table_CON_DEC <- design_table %>% filter(Treatment %in% c('CON', 'DEC'))

design_table_CON_DEC 

sample_info_CON_DEC <- design_table_CON_DEC %>%
  arrange(sample)
Treatment <- as.factor(sample_info_CON_DEC$Treatment)
Replicate <- as.factor(sample_info_CON_DEC$Replicate)
design_CON_DEC <- model.matrix(~Treatment + ~Replicate)

reads_per_gene_mRNAseq_CON_DEC <- mrna_seq_treatment_replicate %>% filter(Treatment %in% c('CON', 'DEC')) %>% select(gene_ID, second_read_strand, sample)

read_counts_CON_DEC <- reads_per_gene_mRNAseq_CON_DEC  %>%
  select(sample, gene_ID, count=second_read_strand) %>%
  spread(sample, count) %>%
  column_to_rownames('gene_ID')

dge_list <- calcNormFactors(DGEList(
  counts = as.matrix(read_counts_CON_DEC),
  group = Treatment))
dge_est <- estimateDisp(dge_list, design_CON_DEC)
fit <- glmQLFit(dge_est, design_CON_DEC)
qlf_CON_DEC <- glmQLFTest(fit, coef=2)
topTags(qlf_CON_DEC)

#add gene_id to column
CON_DEC_qlf_gene_ID <- topTags(qlf_CON_DEC, 1e9)[[1]] %>%
  rownames_to_column('gene_id')  %>%
  as_tibble()

CON_DEC_top10 <- topTags(qlf_CON_DEC)

#remove version to enable GO recognition - make top 10 plot

CON_DEC_qlf_gene_ID_ <- separate(CON_DEC_qlf_gene_ID, col=gene_id, into=c('gene_id'), extra='drop')

CON_DEC_top10 <- CON_DEC_qlf_gene_ID_ %>% head(10)


write.csv(CON_DEC_qlf_gene_ID_, "MEF_SD_DEC/mRNAseq/data/csvs/CON_DEC_qlf_gene_ID_.csv", row.names = FALSE)

write.csv(CON_DEC_top10, "MEF_SD_DEC/mRNAseq/data/csvs/CON_DEC_top10.csv", row.names = FALSE)

#FDR < 0.01

CON_DEC_FDR <- CON_DEC_qlf_gene_ID_ %>% filter(FDR < 0.01)

write.csv(CON_DEC_FDR, "MEF_SD_DEC/mRNAseq/data/csvs/CON_DEC_FDR.csv", row.names = FALSE)

CON_DEC_FDR_pos <- CON_DEC_FDR %>% filter(logFC > 0)

write.csv(CON_DEC_FDR, "MEF_SD_DEC/mRNAseq/data/csvs/CON_DEC_FDR_pos.csv", row.names = FALSE)

CON_DEC_FDR_neg <- CON_DEC_FDR %>% filter(logFC < 0)


write.csv(CON_DEC_FDR, "MEF_SD_DEC/mRNAseq/data/csvs/CON_DEC_FDR_neg.csv", row.names = FALSE)


#upreg and downreg tables for plot

CON_DEC_pos <- CON_DEC_qlf_gene_ID_ %>% filter(logFC > 0)

CON_DEC_neg <- CON_DEC_qlf_gene_ID_ %>% filter(logFC < 0)

write.csv(CON_DEC_pos, "MEF_SD_DEC/mRNAseq/data/csvs/CON_DEC_pos.csv", row.names = FALSE)

write.csv(CON_DEC_neg, "MEF_SD_DEC/mRNAseq/data/csvs/CON_DEC_neg.csv", row.names = FALSE)

#arrange by FDR

CON_DEC_pos <- CON_DEC_pos %>% arrange(FDR)

CON_DEC_neg <- CON_DEC_neg  %>% arrange(FDR)

#top 10 pos and neg

CON_DEC_pos_10 <- CON_DEC_pos %>% head(10)

CON_DEC_neg_10 <- CON_DEC_neg %>% head(10)



##GO

#gene lists

gene_list_neg_con_dec_10 <- CON_DEC_neg_10$gene_id

gene_list_pos_con_dec_10 <- CON_DEC_pos_10$gene_id

#run for CON_DEC neg

ensembl_mm10 <- biomaRt::useEnsembl(biomart = "ensembl",
                                    dataset = 'mmusculus_gene_ensembl')

go_results_neg_CON_DEC <- biomaRt::getBM(
  attributes=c("ensembl_gene_id", "external_gene_name",
               "go_id", "name_1006", "namespace_1003"),
  filters="ensembl_gene_id",
  values=gene_list_neg_con_dec_10,
  mart=ensembl_mm10) %>%
  as_tibble()

#run for Con_dec pos

ensembl_mm10 <- biomaRt::useEnsembl(biomart = "ensembl",
                                    dataset = 'mmusculus_gene_ensembl')

go_results_pos_CON_DEC <- biomaRt::getBM(
  attributes=c("ensembl_gene_id", "external_gene_name",
               "go_id", "name_1006", "namespace_1003"),
  filters="ensembl_gene_id",
  values=gene_list_pos_con_dec_10,
  mart=ensembl_mm10) %>%
  as_tibble()


#plot

library(readr)
go_results_neg_CON_DEC <- read_csv("MEF_SD_DEC/mRNAseq/data/csvs/go_results_neg_CON_DEC.csv")
go_results_neg_CON_DEC

library(readr)
go_results_pos_CON_DEC <- read_csv("MEF_SD_DEC/mRNAseq/data/csvs/go_results_pos_CON_DEC.csv")
go_results_pos_CON_DEC


go_results_neg_CON_DEC_ <- go_results_neg_CON_DEC %>% mutate(namespace_1003=gsub("_"," ", namespace_1003))

go_results_neg_CON_DEC_ %>% filter(!name_1006 %in% c('neuron differentiation', 'nervous system development', 'in utero embryonic development', 'heart development')) %>% group_by(name_1006) %>% filter(name_1006!="")  %>% filter(n()>1)  %>% ggplot() + geom_bar(aes(y=forcats::fct_infreq(name_1006))) +
  facet_grid(rows=vars(namespace_1003), scales='free_y') +
  labs(x='Number of genes', y='GO Term') +
  theme_bw() +
  theme(panel.grid=element_blank(), aspect.ratio = 1/2, plot.title=element_text(size=10)) + 
  ggtitle("Genes that show downregulated expression 
          in Decitabine-exposed MEFs")



go_results_pos_CON_DEC_ <- go_results_pos_CON_DEC %>% mutate(namespace_1003=gsub("_"," ", namespace_1003))


go_results_pos_CON_DEC_ %>% group_by(name_1006) %>% filter(!name_1006 %in% c('calcium-dependent cell-cell adhesion via plasma membrane cell adhesion molecules', 'in utero embryonic development', 'positive regulation of angiogenesis', 'spermatogenesis')) %>% filter(name_1006!="")  %>% filter(n()>1)  %>% ggplot() + geom_bar(aes(y=forcats::fct_infreq(name_1006))) +
  facet_grid(rows=vars(namespace_1003), scales='free_y') +
  labs(x='Number of genes', y='GO Term') +
  scale_x_continuous(limits=c(0, 8)) +
  theme_bw() +
  theme(panel.grid=element_blank(), aspect.ratio = 1/2, plot.title=element_text(size=10)) +
  ggtitle("Genes that show upregulated expression 
          in Decitabine-exposed MEFs")


#write csvs

write.csv(go_results_pos_CON_DEC, "MEF_SD_DEC/mRNAseq/data/csvs/go_results_pos_CON_DEC.csv", row.names = FALSE)

write.csv(go_results_neg_CON_DEC, "MEF_SD_DEC/mRNAseq/data/csvs/go_results_neg_CON_DEC.csv", row.names = FALSE)

###CON_0 and DEC


design_table_CON_0_DEC <- design_table %>% filter(Treatment %in% c('CON_0', 'DEC'))

design_table_CON_0_DEC 

sample_info_CON_0_DEC <- design_table_CON_0_DEC %>%
  arrange(sample)
Treatment <- as.factor(sample_info_CON_0_DEC$Treatment)
Replicate <- as.factor(sample_info_CON_0_DEC$Replicate)
design_CON_0_DEC <- model.matrix(~Treatment + ~Replicate)

reads_per_gene_mRNAseq_CON_0_DEC <- mrna_seq_treatment_replicate %>% filter(Treatment %in% c('CON_0', 'DEC')) %>% select(gene_ID, second_read_strand, sample)

read_counts_CON_0_DEC <- reads_per_gene_mRNAseq_CON_0_DEC  %>%
  select(sample, gene_ID, count=second_read_strand) %>%
  spread(sample, count) %>%
  column_to_rownames('gene_ID')

dge_list <- calcNormFactors(DGEList(
  counts = as.matrix(read_counts_CON_0_DEC),
  group = Treatment))
dge_est <- estimateDisp(dge_list, design_CON_0_DEC)
fit <- glmQLFit(dge_est, design_CON_0_DEC)
qlf_CON_0_DEC <- glmQLFTest(fit, coef=2)
topTags(qlf_CON_0_DEC)

#add gene_id to column
CON_0_DEC_qlf_gene_ID <- topTags(qlf_CON_0_DEC, 1e9)[[1]] %>%
  rownames_to_column('gene_id')  %>%
  as_tibble()

CON_0_DEC_top10 <- topTags(qlf_CON_0_DEC)

#remove version to enable GO recognition - make top 10 plot

CON_0_DEC_qlf_gene_ID_ <- separate(CON_0_DEC_qlf_gene_ID, col=gene_id, into=c('gene_id'), extra='drop')

CON_0_DEC_top10 <- CON_0_DEC_qlf_gene_ID_ %>% head(10)


write.csv(CON_0_DEC_qlf_gene_ID_, "MEF_SD_DEC/mRNAseq/data/csvs/CON_0_DEC_qlf_gene_ID_.csv", row.names = FALSE)

write.csv(CON_0_DEC_top10, "MEF_SD_DEC/mRNAseq/data/csvs/CON_0_DEC_top10.csv", row.names = FALSE)

#FDR < 0.01

CON_0_DEC_FDR <- CON_0_DEC_qlf_gene_ID_ %>% filter(FDR < 0.01)


write.csv(CON_0_DEC_FDR, "MEF_SD_DEC/mRNAseq/data/csvs/CON_0_DEC_FDR.csv", row.names = FALSE)

CON_0_DEC_FDR_pos <- CON_0_DEC_FDR %>% filter(logFC > 0)

write.csv(CON_0_DEC_FDR, "MEF_SD_DEC/mRNAseq/data/csvs/CON_0_DEC_FDR_pos.csv", row.names = FALSE)

CON_0_DEC_FDR_neg <- CON_0_DEC_FDR %>% filter(logFC < 0)


write.csv(CON_0_DEC_FDR, "MEF_SD_DEC/mRNAseq/data/csvs/CON_0_DEC_FDR_neg.csv", row.names = FALSE)


#upreg and downreg tables for plot

CON_0_DEC_pos <- CON_0_DEC_qlf_gene_ID_ %>% filter(logFC > 0)

CON_0_DEC_neg <- CON_0_DEC_qlf_gene_ID_ %>% filter(logFC < 0)

write.csv(CON_0_DEC_pos, "MEF_SD_DEC/mRNAseq/data/csvs/CON_0_DEC_pos.csv", row.names = FALSE)

write.csv(CON_0_DEC_neg, "MEF_SD_DEC/mRNAseq/data/csvs/CON_0_DEC_neg.csv", row.names = FALSE)

#arrange by FDR

CON_0_DEC_pos <- CON_0_DEC_pos %>% arrange(FDR)

CON_0_DEC_neg <- CON_0_DEC_neg  %>% arrange(FDR)

#top 10 pos and neg

CON_0_DEC_pos_10 <- CON_0_DEC_pos %>% head(10)

CON_0_DEC_pos_10 %>% select(gene_id)

CON_0_DEC_neg_10 <- CON_0_DEC_neg %>% head(10)



##GO

#gene lists

gene_list_neg_con_0_dec_10 <- CON_0_DEC_neg_10$gene_id

gene_list_pos_con_0_dec_10 <- CON_0_DEC_pos_10$gene_id

#run for CON_0_DEC neg

ensembl_mm10 <- biomaRt::useEnsembl(biomart = "ensembl",
                                    dataset = 'mmusculus_gene_ensembl', version=106)

go_results_neg_CON_0_DEC <- biomaRt::getBM(
  attributes=c("ensembl_gene_id", "external_gene_name",
               "go_id", "name_1006", "namespace_1003"),
  filters="ensembl_gene_id",
  values=gene_list_neg_con_0_dec_10,
  mart=ensembl_mm10) %>%
  as_tibble()

#run for Con_0_dec pos

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




go_results_neg_CON_0_DEC %>% group_by(name_1006) %>% filter(name_1006!="")  %>% filter(n()>1)  %>% ggplot() + geom_bar(aes(y=forcats::fct_infreq(name_1006))) +
  facet_grid(rows=vars(namespace_1003), scales='free_y') +
  labs(x='Number of genes', y='GO Term') +
  theme_bw() +
  ggtitle("Genes that show downregulated expression in Decitabine-exposed MEFs relative to 
          Day 0 Controls")




go_results_pos_CON_0_DEC %>% group_by(name_1006) %>% filter(name_1006!="")  %>% filter(n()>1)  %>% ggplot() + geom_bar(aes(y=forcats::fct_infreq(name_1006))) +
  facet_grid(rows=vars(namespace_1003), scales='free_y') +
  labs(x='Number of genes', y='GO Term') +
  theme_bw() +
  ggtitle("Genes that show upregulated expression in Decitabine-exposed MEFs relative to Day 0 Controls")


#write csvs

write.csv(go_results_pos_CON_0_DEC, "MEF_SD_DEC/mRNAseq/data/csvs/go_results_pos_CON_0_DEC.csv", row.names = FALSE)

write.csv(go_results_neg_CON_0_DEC, "MEF_SD_DEC/mRNAseq/data/csvs/go_results_neg_CON_0_DEC.csv", row.names = FALSE)

