##rRNA coverage
##rRNA snp_extract
.libPaths('/data/Blizard-Rakyan/Fran/R/packages_4.0.2')
library(tidyverse)

##rRNA coverage (samtools depth)

rrna_fd <- '/data/Blizard-Rakyan/Selin/projects/MEF_rRNA/RNA-Seq/data/covs_rrna/'




list.files(rrna_fd)
#correct 
library(dplyr)
rrna_depth <- bind_rows(lapply(X = list.files(rrna_fd),
                               FUN=function(x) {
                                 read_tsv(paste0(rrna_fd, x), col_names=c('name', 'pos', 'depth', "sample")) %>% mutate(sample = x)
                               }))



rrna_depth_table <- rrna_depth %>% gsub(".cov", "", rrna_depth)
rrna_depth_table

rrna_depth %>% filter(sample == 'Ohr_B.cov')

rrna_depth_table <- rrna_depth %>% mutate(sample =gsub(".cov", "", sample))
rrna_depth_table


#3008 is -1, 3009 is 1. if position is below or equal to 3008, subtract 3009, if it is above 3008 subtract 3008.
library(dplyr)
rrna_depth_table_adj <- rrna_depth_table %>%
  mutate(adj_snp_pos=if_else(pos <= 3008, pos-3009, pos-3008))
rrna_depth_table_adj

rrna_depth_table_adj_hapsnps <- rrna_depth_table_adj %>% filter(adj_snp_pos %in% c(-104, 6832, 6007, 6777, 12736, 12857, 174, 743, 1091))

rrna_depth_table_adj_hapsnps$adj_snp_pos <- as.character(rrna_depth_table_adj_hapsnps$adj_snp_pos)

library(dplyr)
library(stringr)
rrna_depth_table_adj_hapsnps_haplotype <- rrna_depth_table_adj_hapsnps %>%
  mutate(haplotype = case_when(
    str_detect(adj_snp_pos, "-104") ~ "ATA/ATG",
    str_detect(adj_snp_pos, "1091") ~ "CTA",
    str_detect(adj_snp_pos, "12736") ~ "ATG",
    str_detect(adj_snp_pos, "12857") ~ "CTA",
    str_detect(adj_snp_pos, "174") ~ "CCA",
    str_detect(adj_snp_pos, "6007") ~ "ATA",
    str_detect(adj_snp_pos, "6777") ~ "ATG",
    str_detect(adj_snp_pos, "6832") ~ "ATA",
    str_detect(adj_snp_pos, "743") ~ "CTA"))
rrna_depth_table_adj_hapsnps_haplotype


library(ggplot2)
rrna_depth_graph <- rrna_depth_table_adj_hapsnps_haplotype %>% ggplot() +
  aes(x = adj_snp_pos, y = depth) + 
  geom_point(aes(color = sample, shape = haplotype), size = 2, alpha = 1, position = position_dodge(width = 0.35)) +
  theme_bw()+
  labs(x = "Haplotype-defining snp", y = "Depth") +
  ggtitle("rRNA-seq depth at Haplotype-defining snps")
rrna_depth_graph 


write.csv(rrna_depth_table_adj_hapsnps_haplotype, "/data/Blizard-Rakyan/Selin/projects/MEF_rRNA/RNA-Seq/data/csv_tables/rrna_depth_table_adj_hapsnps_haplotype.csv", row.names = FALSE)


rrna_depth_table_adj_sample %>% filter(sample == 'D3ConA')


##rRNA snp_extract

#have a look at onee sample
library(readr)
D3_Con_B_alleles <- read_csv("MEF_rRNA/RNA-Seq/data/snp_extract/D3_Con_B/D3_Con_B_alleles.csv")
D3_Con_B_alleles

#bind tables 
common_path_rRNA = "/data/Blizard-Rakyan/Selin/projects/MEF_rRNA/RNA-Seq/data/snp_extract/"
sub_dirs_rRNA = list.files(common_path_rRNA)
sub_dirs_rRNA

snp_extract_rRNA_table <- bind_rows(lapply(
  X=sub_dirs_rRNA,
  FUN=function(x){read_csv(
    paste0(common_path_rRNA, x, '/', x, '_alleles.csv'),
    col_types=cols()) %>%
      mutate(sample=x)}))
snp_extract_rRNA_table

#summarise valuees for number of reads column
snp_extract_rRNA_table_num_reads<- snp_extract_rRNA_table %>%
  group_by(sample, snp_allele, snp_pos) %>%
  summarise(num_of_reads = n_distinct(read_id))
snp_extract_rRNA_table_num_reads  

#3008 is -1, 3009 is 1. if position is below or equal to 3008, subtract 3009, if it is above 3008 subtract 3008.
library(dplyr)
snp_extract_rRNA_table_num_reads_adj <- snp_extract_rRNA_table_num_reads  %>%
  mutate(adj_snp_pos=if_else(snp_pos <= 3008, snp_pos-3009, snp_pos-3008))
snp_extract_rRNA_table_num_reads_adj

snp_extract_rRNA_6832 <- snp_extract_rRNA_table_num_reads_adj %>% filter(adj_snp_pos %in% 6832)

snp_extract_rRNA_table_num_reads_adj %>% filter(adj_snp_pos %in% 6007)

write.csv(snp_extract_rRNA_6832, "MEF_rRNA/RNA-Seq/data/csv_tables/snp_extract_rRNAseq_6832.csv", row.names = FALSE)

write.csv(snp_extract_rRNA_table_num_reads_adj, "MEF_rRNA/RNA-Seq/data/csv_tables/snp_extract_rRNAseq_table.csv", row.names = FALSE)
