.libPaths('/data/Blizard-Rakyan/Fran/R/packages_4.1.1/')



read_tsv("/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/insert_size_filtered/tsv/N1_0_depth.tsv", col_names=c("ref_seq", "base_index", "depth"))
N1_0_depth <- read_tsv("/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/insert_size_filtered/tsv/N1_0_depth.tsv", col_names=c("ref_seq", "base_index", "depth"))

fd = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/insert_size_filtered/tsv/"
library(dplyr)
read_depth <- bind_rows(lapply(X = list.files(fd),
                 FUN=function(x) {
                   read_tsv(paste0(fd,x),col_names=c("ref_seq", "base_index", "depth")) %>% mutate(sample=x) 
                  }))

read_depth_table <- read_depth %>% gsub("_depth.tsv","", read_depth)
read_depth_table

write.csv(read_depth, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/insert_size_filtered/tsv/read_depth.csv", row.names = FALSE)

library(dplyr)
library(ggplot2)
N1_0_depth_graph <- N1_0_depth %>% ggplot() +
  aes(x = ref_seq, y = depth) + 
  geom_point() + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "sequence name", y = "read depth") +
  ggtitle("Read depth - insert size 150bp for N1_0")+
  theme(legend.position='none')
N1_0_depth_graph

N1_0_rDNA_depth_graph <- N1_0_depth %>% filter(ref_seq == "BK000964.3_looped_3008") %>% ggplot() +
  aes(x=base_index, y = depth) +
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "rDNA base index", y = "read depth") +
  ggtitle("Read depth - rDNA - insert size 150bp for N1_0")+
  theme(legend.position='none')
N1_0_rDNA_depth_graph

##read tsv for mapped and filtered reads (150bp)

#idxstats

fd_idxstats = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/insert_size_filtered_masked_mouse/tsv/"

library(dplyr)
masked_filtered_idxstats <- bind_rows(lapply(X = list.files(fd_idxstats),
                               FUN=function(x) {
                                 read_tsv(paste0(fd_idxstats,x),col_names=c("ref_seq", "seq_length", "mapped_read_segments", "unmapped_read_segments")) %>% mutate(sample=x) 
                               }))

masked_filtered_idxstats_table <- masked_filtered_idxstats %>% group_by(sample) %>% mutate(sample=gsub("_idxstats.tsv","", sample))
masked_filtered_idxstats_table

write.csv(masked_filtered_idxstats_table, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/mapped_filtered_inserts_150/idxstats_150.csv", row.names = FALSE)

#depth

fd_depth = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/insert_size_filtered_masked_mouse/tsv_depth/"
library(dplyr)
masked_filtered_read_depth <- bind_rows(lapply(X = list.files(fd_depth),
                                               FUN=function(x) {
                                                 read_tsv(paste0(fd_depth,x),col_names=c("ref_seq", "base_index", "depth")) %>% mutate(sample=x) 
                                               }))

masked_filtered_read_depth_table <- masked_filtered_read_depth %>% group_by(sample) %>% mutate(sample=gsub("_depth.tsv","", sample))
masked_filtered_read_depth_table

write.csv(masked_filtered_read_depth_table, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/mapped_filtered_inserts_150/depth_150.csv", row.names = FALSE)


#calculate ratio betweeen rDNA and total 

mf_idxstats_ratio <- masked_filtered_idxstats_table %>% group_by(sample) %>% mutate(total_reads = sum(mapped_read_segments)) %>% mutate(rdna_reads = sum(mapped_read_segments[ref_seq == "BK000964.3_looped_3008"])) %>% mutate(rdna_ratio = rdna_reads/total_reads)

mf_idxstats_ratio

write.csv(mf_idxstats_ratio, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/mapped_filtered_inserts_150/rdna_ratio_idxstats_150.csv", row.names = FALSE)

#Treatment condition and replicate

mf_idxstats_ratio_treatment <- mf_idxstats_ratio %>%
  mutate(Treatment = case_when(
    endsWith(sample, "0") ~ "CON_0",
    endsWith(sample, "CON") ~ "CON",
    endsWith(sample, "VEH") ~ "VEH",
    endsWith(sample, "DEC") ~ "DEC",
    endsWith(sample, "S") ~ "SDEP",
    endsWith(sample, "S_D") ~ "SDEP_DEC"
  ))

mf_idxstats_ratio_treatment

mf_idxstats_ratio_treatment_replicate <- mf_idxstats_ratio_treatment %>%
  mutate(Replicate = case_when(
    startsWith(sample, "N1") ~ "1",
    startsWith(sample, "N2") ~ "2",
    startsWith(sample, "N3") ~ "3",
    startsWith(sample, "N4") ~ "4"
  ))
mf_idxstats_ratio_treatment_replicate


library(forcats)
library(dplyr)
library(ggplot2)
idxstats_rdna_ratio_graph <- mf_idxstats_ratio_treatment_replicate %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% ggplot() +
  aes(x = Treatment, y = rdna_ratio) + 
  geom_point(aes(colour = Replicate), size = 2, alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Treatment condition", y = "rdna_ratio") +
  theme_bw()+
  ggtitle("rDNA ratio - insert size 150bp")
idxstats_rdna_ratio_graph


library(forcats)
library(ggplot2)
mf_idxstats_ratio_treatment_replicate %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% ggplot() + aes(x= Treatment , y= rdna_ratio) +
  geom_jitter(aes(color=Replicate), width=0.1, alpha=0.5) +
  ggpubr::stat_compare_means(paired = TRUE, aes(label = paste0("p = ", ..p.format..)), ref.group = 'CON_0',
                             hide.ns=T,
                             bracket.size=1, label.y=0.002) +
  theme_bw() +
  labs(x="Treatment Condition", y="rDNA ratio (rDNA/total reads)") +
  coord_cartesian(ylim=c(0, 0.002)) +
  ggtitle("rDNA ratio - insert size 150bp")


#add Treatment condition


mf_idxstats_ratio_treatment_replicate



library(dplyr)


detach("package:methylKit")

mf_idxstats_reads_ratio <- mf_idxstats_ratio_treatment_replicate %>% select(sample, total_reads, rdna_reads, rdna_ratio) %>% group_by(sample)

mf_idxstats_reads_ratio <- unique(mf_idxstats_reads_ratio)


write.csv(mf_idxstats_reads_ratio, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/mapped_filtered_inserts_150/mf_idxstats_reads_ratio_150.csv", row.names = FALSE)


##read tsv for 140_160 depth table

fd_range_rdna = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/bam_rdna_140_160/tsv_depth/"
library(dplyr)
read_depth_rdna_140_160 <- bind_rows(lapply(X = list.files(fd_range_rdna),
                               FUN=function(x) {
                                 read_tsv(paste0(fd_range_rdna,x),col_names=c("ref_seq", "base_index", "depth")) %>% mutate(sample=x) 
                               }))

read_depth_rdna_140_160_table <- read_depth_rdna_140_160 %>% group_by(sample) %>% mutate(sample=gsub("_depth.tsv","", sample))
read_depth_rdna_140_160_table

read_depth_rdna_140_160_table %>% mutate(base_index = base_index - 3008)


#Treatment condition and replicate

read_depth_rdna_range_treatment <- read_depth_rdna_140_160_table %>%
  mutate(Treatment = case_when(
    endsWith(sample, "0") ~ "CON_0",
    endsWith(sample, "CON") ~ "CON",
    endsWith(sample, "VEH") ~ "VEH",
    endsWith(sample, "DEC") ~ "DEC",
    endsWith(sample, "S") ~ "SDEP",
    endsWith(sample, "S_D") ~ "SDEP_DEC"
  ))

read_depth_rdna_range_treatment

read_depth_rdna_range_treatment_replicate <- read_depth_rdna_range_treatment%>%
  mutate(Replicate = case_when(
    startsWith(sample, "N1") ~ "1",
    startsWith(sample, "N2") ~ "2",
    startsWith(sample, "N3") ~ "3",
    startsWith(sample, "N4") ~ "4"
  ))
read_depth_rdna_range_treatment_replicate

read_depth_rdna_range_minus3008 <- read_depth_rdna_range_treatment_replicate %>% mutate(base_index = base_index - 3008)
read_depth_rdna_range_minus3008 

library(dplyr)
library(ggplot2)
rDNA_all_range_depth_graph <- read_depth_rdna_range_minus3008  %>% ggplot() +
  aes(x=base_index, y = depth) +
  geom_rect(data = read_depth_rdna_range_minus3008, aes(xmin=4008, xmax=5877, ymin = -Inf, ymax=Inf), fill= "grey84", alpha = 0.05)+
  geom_rect(data = read_depth_rdna_range_minus3008, aes(xmin=8123, xmax=12852, ymin = -Inf, ymax=Inf), fill= "grey84", alpha = 0.05)+
  geom_rect(data = read_depth_rdna_range_minus3008, aes(xmin=6878, xmax=7034, ymin = -Inf, ymax=Inf), fill= "grey84", alpha = 0.05)+
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "rDNA base index", y = "read depth") +
  ggtitle("Read depth - 140-160 bp range")+
  facet_grid(vars(Treatment), vars(Replicate)) +
  theme(strip.text.y = element_text(size = 5)) +
  theme(aspect.ratio=1) +
  theme(legend.position='none')
rDNA_all_range_depth_graph


library(dplyr)
library(ggplot2)
rDNA_all_range_depth_graph <- read_depth_rdna_range_minus3008  %>% ggplot() +
  aes(x=base_index, y = depth) +
  geom_rect(data = read_depth_rdna_range_minus3008, aes(xmin=4008, xmax=5877, ymin = -Inf, ymax=Inf), fill= "grey84", alpha = 0.05)+
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "rDNA base index", y = "read depth") +
  ggtitle("Read depth - 140-160 bp range")+
  facet_grid(vars(Treatment), vars(Replicate)) +
  theme(strip.text.y = element_text(size = 5)) +
  theme(aspect.ratio=1) +
  theme(legend.position='none')
rDNA_all_range_depth_graph


rDNA_range_depth_graph <- read_depth_rdna_range_treatment_replicate %>% ggplot() +
  aes(x=base_index, y = depth) +
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "rDNA base index", y = "read depth (log10 scale)") +
  ggtitle("Read depth - rDNA - insert range 140-160bp")+
  theme(legend.position='none') +
  facet_grid(rows = vars(Treatment), cols = vars(Replicate)) +
  theme(aspect.ratio = 1) + 
scale_y_log10()
rDNA_range_depth_graph

#unfiltered depth

fd_unfil_rdna = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/bam_rdna/tsv_depth/"
library(dplyr)
read_depth_rdna_unfil <- bind_rows(lapply(X = list.files(fd_unfil_rdna),
                                            FUN=function(x) {
                                              read_tsv(paste0(fd_unfil_rdna,x),col_names=c("ref_seq", "base_index", "depth")) %>% mutate(sample=x) 
                                            }))

read_depth_rdna_unfil_table <- read_depth_rdna_unfil %>% group_by(sample) %>% mutate(sample=gsub("_depth.tsv","", sample))
read_depth_rdna_unfil_table

read_depth_rdna_unfil_table %>% mutate(base_index = base_index - 3008)

#Treatment condition and replicate

read_depth_rdna_unfil_treatment <- read_depth_rdna_unfil_table %>%
  mutate(Treatment = case_when(
    endsWith(sample, "0") ~ "CON_0",
    endsWith(sample, "CON") ~ "CON",
    endsWith(sample, "VEH") ~ "VEH",
    endsWith(sample, "DEC") ~ "DEC",
    endsWith(sample, "S") ~ "SDEP",
    endsWith(sample, "S_D") ~ "SDEP_DEC"
  ))

read_depth_rdna_unfil_treatment

read_depth_rdna_unfil_treatment_replicate <- read_depth_rdna_unfil_treatment%>%
  mutate(Replicate = case_when(
    startsWith(sample, "N1") ~ "1",
    startsWith(sample, "N2") ~ "2",
    startsWith(sample, "N3") ~ "3",
    startsWith(sample, "N4") ~ "4"
  ))
read_depth_rdna_unfil_treatment_replicate

read_depth_rdna_unfil_minus3008 <- read_depth_rdna_unfil_treatment_replicate %>% mutate(base_index = base_index - 3008)
read_depth_rdna_unfil_minus3008 


library(dplyr)
library(ggplot2)
rDNA_all_unfil_depth_graph <- read_depth_rdna_unfil_minus3008  %>% ggplot() +
  aes(x=base_index, y = depth) +
  geom_rect(data = read_depth_rdna_unfil_minus3008, aes(xmin=4008, xmax=5877, ymin = -Inf, ymax=Inf), fill= "grey84", alpha = 0.05)+
  geom_rect(data = read_depth_rdna_unfil_minus3008, aes(xmin=8123, xmax=12852, ymin = -Inf, ymax=Inf), fill= "grey84", alpha = 0.05)+
  geom_rect(data = read_depth_rdna_unfil_minus3008, aes(xmin=6878, xmax=7034, ymin = -Inf, ymax=Inf), fill= "grey84", alpha = 0.05)+
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "rDNA base index", y = "read depth") +
  ggtitle("Read depth - unfiltered inserts")+
  facet_grid(vars(Treatment), vars(Replicate)) +
  theme(strip.text.y = element_text(size = 5)) +
  theme(aspect.ratio=1) +
  theme(legend.position='none')
rDNA_all_unfil_depth_graph


library(dplyr)
library(ggplot2)
rDNA_all_unfil_depth_graph <- read_depth_rdna_unfil_minus3008  %>% ggplot() +
  aes(x=base_index, y = depth_unfil) +
  geom_rect(data = read_depth_rdna_unfil_minus3008, aes(xmin=4008, xmax=5877, ymin = -Inf, ymax=Inf), fill= "grey84", alpha = 0.05)+
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "rDNA base index", y = "read depth") +
  ggtitle("Read depth - unfiltered inserts")+
  facet_grid(vars(Treatment), vars(Replicate)) +
  theme(strip.text.y = element_text(size = 5)) +
  theme(aspect.ratio=1) +
  theme(legend.position='none')
rDNA_all_unfil_depth_graph





#30_51 depth

fd_30_51_rdna = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/bam_rdna_30_51/tsv_depth/"
library(dplyr)
read_depth_rdna_30_51 <- bind_rows(lapply(X = list.files(fd_30_51_rdna),
                                          FUN=function(x) {
                                            read_tsv(paste0(fd_30_51_rdna,x),col_names=c("ref_seq", "base_index", "depth")) %>% mutate(sample=x) 
                                          }))

read_depth_rdna_30_51_table <- read_depth_rdna_30_51 %>% group_by(sample) %>% mutate(sample=gsub("_depth.tsv","", sample))
read_depth_rdna_30_51_table

read_depth_rdna_30_51_table %>% mutate(base_index = base_index - 3008)


#Treatment condition and replicate

read_depth_rdna_30_51_treatment <- read_depth_rdna_30_51_table %>%
  mutate(Treatment = case_when(
    endsWith(sample, "0") ~ "CON_0",
    endsWith(sample, "CON") ~ "CON",
    endsWith(sample, "VEH") ~ "VEH",
    endsWith(sample, "DEC") ~ "DEC",
    endsWith(sample, "S") ~ "SDEP",
    endsWith(sample, "S_D") ~ "SDEP_DEC"
  ))

read_depth_rdna_30_51_treatment

read_depth_rdna_30_51_treatment_replicate <- read_depth_rdna_30_51_treatment%>%
  mutate(Replicate = case_when(
    startsWith(sample, "N1") ~ "1",
    startsWith(sample, "N2") ~ "2",
    startsWith(sample, "N3") ~ "3",
    startsWith(sample, "N4") ~ "4"
  ))
read_depth_rdna_30_51_treatment_replicate

read_depth_rdna_30_51_minus3008 <- read_depth_rdna_30_51_treatment_replicate %>% mutate(base_index = base_index - 3008)
read_depth_rdna_30_51_minus3008 


library(dplyr)
library(ggplot2)
rDNA_all_30_51_depth_graph <- read_depth_rdna_30_51_minus3008  %>% ggplot() +
  aes(x=base_index, y = depth) +
  geom_rect(data = read_depth_rdna_30_51_minus3008, aes(xmin=4008, xmax=5877, ymin = -Inf, ymax=Inf), fill= "grey84", alpha = 0.05)+
  geom_rect(data = read_depth_rdna_30_51_minus3008, aes(xmin=8123, xmax=12852, ymin = -Inf, ymax=Inf), fill= "grey84", alpha = 0.05)+
  geom_rect(data = read_depth_rdna_30_51_minus3008, aes(xmin=6878, xmax=7034, ymin = -Inf, ymax=Inf), fill= "grey84", alpha = 0.05)+
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "rDNA base index", y = "read depth") +
  ggtitle("Read depth - 30-51")+
  facet_grid(vars(Treatment), vars(Replicate)) +
  theme(strip.text.y = element_text(size = 5)) +
  theme(aspect.ratio=1) +
  theme(legend.position='none')
rDNA_all_30_51_depth_graph

library(dplyr)
library(ggplot2)
rDNA_all_30_51_depth_graph <- read_depth_rdna_30_51_minus3008  %>% ggplot() +
  aes(x=base_index, y = depth_30_51) +
  geom_rect(data = read_depth_rdna_30_51_minus3008, aes(xmin=4008, xmax=5877, ymin = -Inf, ymax=Inf), fill= "grey84", alpha = 0.05)+
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "rDNA base index", y = "read depth") +
  ggtitle("Read depth - 30-51")+
  facet_grid(vars(Treatment), vars(Replicate)) +
  theme(strip.text.y = element_text(size = 5)) +
  theme(aspect.ratio=1) +
  theme(legend.position='none')
rDNA_all_30_51_depth_graph


#bind 140_160_30_51_unfil

allranges_depth <- bind_rows(
  
  read_depth_rdna_unfil_minus3008 %>%
    mutate(mode = 'unfil'),
  read_depth_rdna_30_51_minus3008 %>%
    mutate(mode = '30-51bp'),
  read_depth_rdna_range_minus3008 %>%
    mutate(mode = '140-160bp'))


filter(meth + unmeth > 100, start >= 4008+3009, end <=5877+3009)

library(dplyr)
library(ggplot2)
rDNA_all_depth_graph <- allranges_depth %>% filter(base_index >= 4008, base_index <=5877) %>% ggplot() +
  aes(x=base_index, y = depth) +
  geom_rect(data = allranges_depth, aes(xmin=4008, xmax=5877, ymin = -Inf, ymax=Inf), fill= "grey84", alpha = 0.05)+
  geom_line(aes(colour = mode))+
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "rDNA base index", y = "read depth") +
  ggtitle("Read depth : filtered and unfiltered")+
  facet_grid(vars(Treatment), vars(Replicate)) +
  theme(strip.text.y = element_text(size = 5)) +
  theme(aspect.ratio=1)+
  scale_x_continuous(limits = c(4005, 5880)) +
  theme_bw()
rDNA_all_depth_graph

library(dplyr)
library(ggplot2)
rDNA_all_depth_graph <- allranges_depth %>% filter(base_index >= 4008, base_index <=5877) %>% ggplot() +
  aes(x=base_index, y = depth) +
  geom_line(aes(colour = mode))+
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "rDNA base index", y = "read depth") +
  ggtitle("Read depth at 18S : filtered and unfiltered")+
  facet_grid(vars(Treatment), vars(Replicate)) +
  theme(strip.text.y = element_text(size = 5)) +
  theme(aspect.ratio=1)+
  scale_x_continuous(limits = c(4008, 5877)) +
  theme_bw()
rDNA_all_depth_graph



#range 140-160 rdna idxstats

fd_rdna_range_idxstats = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/bam_rdna_140_160/tsv_idxstats/"

library(dplyr)
rdna_range_filtered_idxstats <- bind_rows(lapply(X = list.files(fd_rdna_range_idxstats),
                                             FUN=function(x) {
                                               read_tsv(paste0(fd_rdna_range_idxstats,x),col_names=c("ref_seq", "seq_length", "mapped_read_segments", "unmapped_read_segments")) %>% mutate(sample=x) 
                                             }))

rdna_range_filtered_idxstats_table <- rdna_range_filtered_idxstats %>% group_by(sample) %>% mutate(sample=gsub("_idxstats.tsv","", sample))

rdna_range_filtered_idxstats_table

rdna_reads_number <- rdna_range_filtered_idxstats_table %>% group_by(sample) %>% mutate(total_reads = sum(mapped_read_segments)) %>% mutate(rdna_reads = sum(mapped_read_segments[ref_seq == "BK000964.3_looped_3008"])) 

rdna_reads_num_table <- rdna_reads_number %>% filter(ref_seq == "BK000964.3_looped_3008") %>% select(sample, rdna_reads)

#whole genome 140_160

fd_wg_range_idxstats = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/bam_mask_mouse_140_160/tsv_idxstats/"

library(dplyr)
wg_range_filtered_idxstats <- bind_rows(lapply(X = list.files(fd_wg_range_idxstats),
                                                 FUN=function(x) {
                                                   read_tsv(paste0(fd_wg_range_idxstats,x),col_names=c("ref_seq", "seq_length", "mapped_read_segments", "unmapped_read_segments")) %>% mutate(sample=x) 
                                                 }))

wg_range_filtered_idxstats_table <- wg_range_filtered_idxstats %>% group_by(sample) %>% mutate(sample=gsub("_idxstats.tsv","", sample))
#view below:
wg_range_filtered_idxstats_table

wg_range_idxstats_ratio <- wg_range_filtered_idxstats_table %>% group_by(sample) %>% mutate(total_reads = sum(mapped_read_segments)) %>% mutate(rdna_reads = sum(mapped_read_segments[ref_seq == "BK000964.3_looped_3008"])) %>% mutate(rdna_ratio = rdna_reads/total_reads)

wg_range_idxstats_ratio

wg_range_idxstats_ratio_select <- wg_range_idxstats_ratio %>% select(sample, total_reads, rdna_reads, rdna_ratio) 

wg_140_160_ratio <- unique(wg_range_idxstats_ratio_select)

write.csv(wg_140_160_ratio, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/140_160/140_160_ratio.csv", row.names = FALSE)

##range coding region

fd_rdna_cr = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/bam_rdna_140_160/tsv_coding_region/"

library(dplyr)
rdna_codingregion_count <- bind_rows(lapply(X = list.files(fd_rdna_cr),
                                               FUN=function(x) {
                                                 read_tsv(paste0(fd_rdna_cr,x),col_names=c("rdna_coding_region_reads")) %>% mutate(sample=x) 
                                               }))

rdna_codingregion_count_table <- rdna_codingregion_count %>% group_by(sample) %>% mutate(sample=gsub("_codingregion.tsv","", sample))
rdna_codingregion_count_table 

wg_rdna_cr_table <- inner_join(rdna_codingregion_count_table, wg_140_160_ratio, by = "sample")
wg_rdna_cr_table

wg_rdna_cr_ratio_table <- wg_rdna_cr_table %>% mutate(rdna_cr_ratio = rdna_coding_region_reads/total_reads)

col_order <- c("sample", "total_reads", "rdna_reads", "rdna_coding_region_reads", "rdna_ratio", "rdna_cr_ratio")
wg_rdna_cr_ratio_table_orderattempt <- wg_rdna_cr_ratio_table[, col_order]
wg_rdna_cr_ratio_table_orderattempt

wg_rdna_cr_ratio_table_order <- wg_rdna_cr_ratio_table_orderattempt

write.csv(wg_rdna_cr_ratio_table_order, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/140_160/wg_rdna_cr_ratio_140_160.csv", row.names = FALSE)

#Treatment condition and replicate

wg_rdna_cr_ratio_table_order_treatment <- wg_rdna_cr_ratio_table_order %>%
  mutate(Treatment = case_when(
    endsWith(sample, "0") ~ "CON_0",
    endsWith(sample, "CON") ~ "CON",
    endsWith(sample, "VEH") ~ "VEH",
    endsWith(sample, "DEC") ~ "DEC",
    endsWith(sample, "S") ~ "SDEP",
    endsWith(sample, "S_D") ~ "SDEP_DEC"
  ))

wg_rdna_cr_ratio_table_order_treatment

wg_rdna_cr_ratio_table_order_treatment_replicate <- wg_rdna_cr_ratio_table_order_treatment %>%
  mutate(Replicate = case_when(
    startsWith(sample, "N1") ~ "1",
    startsWith(sample, "N2") ~ "2",
    startsWith(sample, "N3") ~ "3",
    startsWith(sample, "N4") ~ "4"
  ))
wg_rdna_cr_ratio_table_order_treatment_replicate

library(forcats)
library(dplyr)
library(ggplot2)
rdna_ratio_140_160_graph <- wg_rdna_cr_ratio_table_order_treatment_replicate %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% ggplot() +
  aes(x = Treatment, y = rdna_ratio) + 
  geom_point(aes(colour = Replicate), size = 2, alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Treatment condition", y = "rdna ratio (rdna reads/total reads") +
  theme_bw()+
  ggtitle("rDNA ratio - insert size range 140-160bp")
rdna_ratio_140_160_graph 

library(forcats)
library(dplyr)
library(ggplot2)
rdna_cr_ratio_140_160_graph <- wg_rdna_cr_ratio_table_order_treatment_replicate %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% ggplot() +
  aes(x = Treatment, y = rdna_cr_ratio) + 
  geom_point(aes(colour = Replicate), size = 2, alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Treatment condition", y = "rdna coding region ratio (rdna coding region reads/total reads") +
  theme_bw()+
  ggtitle("rDNA coding region ratio - insert size range 140-160bp")
rdna_cr_ratio_140_160_graph 

##human

#whole genome 140_160

fd_human_140_160_idxstats = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/Human/bam_mask_140_160/tsv_idxstats/"

library(dplyr)
human_140_160_idxstats <- bind_rows(lapply(X = list.files(fd_human_140_160_idxstats),
                                               FUN=function(x) {
                                                 read_tsv(paste0(fd_human_140_160_idxstats,x),col_names=c("ref_seq", "seq_length", "mapped_read_segments", "unmapped_read_segments")) %>% mutate(sample=x) 
                                               }))

human_140_160_idxstats_table <- human_140_160_idxstats %>% group_by(sample) %>% mutate(sample=gsub("_idxstats.tsv","", sample))
#view below:
human_140_160_idxstats_table

human_140_160_idxstats_ratio <- human_140_160_idxstats_table %>% group_by(sample) %>% mutate(total_reads = sum(mapped_read_segments)) %>% mutate(rdna_reads = sum(mapped_read_segments[ref_seq == "KY962518.1_looped_2120"])) %>% mutate(rdna_ratio = rdna_reads/total_reads)

human_140_160_idxstats_ratio

human_140_160_idxstats_ratio_select <- human_140_160_idxstats_ratio %>% select(sample, total_reads, rdna_reads, rdna_ratio) 
human_140_160_idxstats_ratio_select

human_140_160_ratio <- unique(human_140_160_idxstats_ratio_select)
human_140_160_ratio

write.csv(human_140_160_ratio, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/140_160/human_140_160_ratio.csv", row.names = FALSE)

#whole genome 130_170

fd_human_130_170_idxstats = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/Human/bam_mask_130_170/tsv_idxstats/"

library(dplyr)
human_130_170_idxstats <- bind_rows(lapply(X = list.files(fd_human_130_170_idxstats),
                                           FUN=function(x) {
                                             read_tsv(paste0(fd_human_130_170_idxstats,x),col_names=c("ref_seq", "seq_length", "mapped_read_segments", "unmapped_read_segments")) %>% mutate(sample=x) 
                                           }))

human_130_170_idxstats_table <- human_130_170_idxstats %>% group_by(sample) %>% mutate(sample=gsub("_idxstats.tsv","", sample))
#view below:
human_130_170_idxstats_table

human_130_170_idxstats_ratio <- human_130_170_idxstats_table %>% group_by(sample) %>% mutate(total_reads = sum(mapped_read_segments)) %>% mutate(rdna_reads = sum(mapped_read_segments[ref_seq == "KY962518.1_looped_2120"])) %>% mutate(rdna_ratio = rdna_reads/total_reads)

human_130_170_idxstats_ratio

human_130_170_idxstats_ratio_select <- human_130_170_idxstats_ratio %>% select(sample, total_reads, rdna_reads, rdna_ratio) 
human_130_170_idxstats_ratio_select

human_130_170_ratio <- unique(human_130_170_idxstats_ratio_select)
human_130_170_ratio

write.csv(human_130_170_ratio, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/130_170_human/human_130_170_ratio.csv", row.names = FALSE)

#coding region 140_160

fd_human_rdna_cr_140_160 = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/Human/bam_mask_140_160/tsv_coding_region/"

library(dplyr)
human_140_160_rdna_codingregion <- bind_rows(lapply(X = list.files(fd_human_rdna_cr_140_160),
                                            FUN=function(x) {
                                              read_tsv(paste0(fd_human_rdna_cr_140_160,x),col_names=c("rdna_coding_region_reads")) %>% mutate(sample=x) 
                                            }))

human_140_160_rdna_codingregion_table <- human_140_160_rdna_codingregion %>% group_by(sample) %>% mutate(sample=gsub("_codingregion.tsv","", sample))
human_140_160_rdna_codingregion_table 

human_140_160_rdna_cr_table <- inner_join(human_140_160_rdna_codingregion_table, human_140_160_ratio, by = "sample")
human_140_160_rdna_cr_table

human_140_160_rdna_cr_ratio_table <- human_140_160_rdna_cr_table %>% mutate(rdna_cr_ratio = rdna_coding_region_reads/total_reads)
human_140_160_rdna_cr_ratio_table


col_order <- c("sample", "total_reads", "rdna_reads", "rdna_coding_region_reads", "rdna_ratio", "rdna_cr_ratio")
human_140_160_rdna_cr_ratio_table_orderattempt <- human_140_160_rdna_cr_ratio_table[, col_order]
human_140_160_rdna_cr_ratio_table_orderattempt

human_140_160_rdna_cr_ratio_table_order <- human_140_160_rdna_cr_ratio_table_orderattempt

write.csv(human_140_160_rdna_cr_ratio_table_order, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/140_160/human_rdna_cr_ratio_140_160.csv", row.names = FALSE)

human_140_160_rdna_cr_ratio_table_order_treatment <- human_140_160_rdna_cr_ratio_table_order %>%
  mutate(Treatment = case_when(
    endsWith(sample, "Dep") ~ "SDEP",
    endsWith(sample, "Con") ~ "CON"
  ))
human_140_160_rdna_cr_ratio_table_order_treatment



#coding region 130_170

fd_human_rdna_cr_130_170 = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/Human/bam_mask_130_170/tsv_coding_region/"

library(dplyr)
human_130_170_rdna_codingregion <- bind_rows(lapply(X = list.files(fd_human_rdna_cr_130_170),
                                                    FUN=function(x) {
                                                      read_tsv(paste0(fd_human_rdna_cr_130_170,x),col_names=c("rdna_coding_region_reads")) %>% mutate(sample=x) 
                                                    }))

human_130_170_rdna_codingregion_table <- human_130_170_rdna_codingregion %>% group_by(sample) %>% mutate(sample=gsub("_codingregion.tsv","", sample))
human_130_170_rdna_codingregion_table 

human_130_170_rdna_cr_table <- inner_join(human_130_170_rdna_codingregion_table, human_130_170_ratio, by = "sample")
human_130_170_rdna_cr_table

human_130_170_rdna_cr_ratio_table <- human_130_170_rdna_cr_table %>% mutate(rdna_cr_ratio = rdna_coding_region_reads/total_reads)
human_130_170_rdna_cr_ratio_table


col_order <- c("sample", "total_reads", "rdna_reads", "rdna_coding_region_reads", "rdna_ratio", "rdna_cr_ratio")
human_130_170_rdna_cr_ratio_table_orderattempt <- human_130_170_rdna_cr_ratio_table[, col_order]
human_130_170_rdna_cr_ratio_table_orderattempt

human_130_170_rdna_cr_ratio_table_order <- human_130_170_rdna_cr_ratio_table_orderattempt

write.csv(human_130_170_rdna_cr_ratio_table_order, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/130_170_human/human_rdna_cr_ratio_130_170.csv", row.names = FALSE)

human_130_170_rdna_cr_ratio_table_order_treatment <- human_130_170_rdna_cr_ratio_table_order %>%
  mutate(Treatment = case_when(
    endsWith(sample, "Dep") ~ "SDEP",
    endsWith(sample, "Con") ~ "CON"
  ))
human_130_170_rdna_cr_ratio_table_order_treatment

library(forcats)
library(dplyr)
library(ggplot2)
human_130_170_graph <- human_130_170_rdna_cr_ratio_table_order_treatment %>% mutate(Treatment=fct_relevel(Treatment, 'CON')) %>% ggplot() +
  aes(x = Treatment, y = rdna_cr_ratio) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Treatment condition", y = "rdna coding region ratio (rdna coding region reads/total reads") +
  theme_bw()+
  ggtitle("Human rDNA coding region ratio - insert size range 130-170bp")
human_130_170_graph

library(forcats)
library(dplyr)
library(ggplot2)
human_140_160_graph <- human_140_160_rdna_cr_ratio_table_order_treatment %>% mutate(Treatment=fct_relevel(Treatment, 'CON')) %>% ggplot() +
  aes(x = Treatment, y = rdna_cr_ratio) + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Treatment condition", y = "rdna coding region ratio (rdna coding region reads/total reads") +
  theme_bw()+
  ggtitle("Human rDNA coding region ratio - insert size range 140-160bp")
human_140_160_graph


##reads_per_gene 

#mRNAseq

fd_star_mRNAseq = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/mRNAseq/data/star/"

library(dplyr)
reads_per_gene_mRNAseq <- bind_rows(lapply(X = list.files(fd_star_mRNAseq), 
                                             FUN=function(x){ fn <- list.files(paste0(fd_star_mRNAseq,x), pattern = "ReadsPerGene.out.tab")
                                               read_tsv(paste0(fd_star_mRNAseq,x,'/',fn),col_names=c("gene_ID", "unstranded", "first_read_strand", "second_read_strand")) %>% mutate(sample=x) 
                                             }))

write.csv(reads_per_gene_mRNAseq, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/mRNAseq/data/star/reads_per_gene_mRNAseq.csv", row.names = FALSE)

reads_per_gene_mRNAseq %>% filter(gene_ID == "rDNA")

#rRNAseq

fd_star_rRNAseq = "/data/Blizard-Rakyan/Selin/projects/MEF_rRNA/RNA-Seq/data/star/"

library(dplyr)
reads_per_gene_rRNAseq <- bind_rows(lapply(X = list.files(fd_star_rRNAseq), 
                                           FUN=function(x){ fn <- list.files(paste0(fd_star_rRNAseq,x), pattern = "ReadsPerGene.out.tab")
                                           read_tsv(paste0(fd_star_rRNAseq,x,'/',fn),col_names=c("gene_ID", "unstranded", "first_read_strand", "second_read_strand")) %>% mutate(sample=x) 
                                           }))
rDNA_reads_rRNAseq <- reads_per_gene_rRNAseq %>% filter(gene_ID == "rDNA")

reads_per_gene_mRNAseq %>% filter(gene_ID == "rDNA")

###unfiltered reads

##read tsv

#idxstats

fd_idxstats = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/bam_mask_mouse/tsv_idxstats/"

library(dplyr)
mask_idxstats <- bind_rows(lapply(X = list.files(fd_idxstats),
                                             FUN=function(x) {
                                               read_tsv(paste0(fd_idxstats,x),col_names=c("ref_seq", "seq_length", "mapped_read_segments", "unmapped_read_segments")) %>% mutate(sample=x) 
                                             }))

mask_idxstats_table <- mask_idxstats %>% group_by(sample) %>% mutate(sample=gsub("_idxstats.tsv","", sample))
mask_idxstats_table

write.csv(mask_idxstats_table, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/mask_not_filtered_rDNAratio/idxstats_nf.csv", row.names = FALSE)


#calculate ratio betweeen rDNA and total 

mask_idxstats_ratio <- mask_idxstats_table %>% group_by(sample) %>% mutate(total_reads = sum(mapped_read_segments)) %>% mutate(rdna_reads = sum(mapped_read_segments[ref_seq == "BK000964.3_looped_3008"])) %>% mutate(rdna_ratio = rdna_reads/total_reads)

mask_idxstats_ratio

write.csv(mask_idxstats_ratio, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/mask_not_filtered_rDNAratio/rdna_ratio_idxstats_nf.csv", row.names = FALSE)

#Treatment condition and replicate

mask_idxstats_ratio_treatment <- mask_idxstats_ratio %>%
  mutate(Treatment = case_when(
    endsWith(sample, "0") ~ "CON_0",
    endsWith(sample, "CON") ~ "CON",
    endsWith(sample, "VEH") ~ "VEH",
    endsWith(sample, "DEC") ~ "DEC",
    endsWith(sample, "S") ~ "SDEP",
    endsWith(sample, "S_D") ~ "SDEP_DEC"
  ))

mask_idxstats_ratio_treatment

mask_idxstats_ratio_treatment_replicate <- mask_idxstats_ratio_treatment %>%
  mutate(Replicate = case_when(
    startsWith(sample, "N1") ~ "1",
    startsWith(sample, "N2") ~ "2",
    startsWith(sample, "N3") ~ "3",
    startsWith(sample, "N4") ~ "4"
  ))
mask_idxstats_ratio_treatment_replicate


library(forcats)
library(dplyr)
library(ggplot2)
idxstats_rdna_ratio_graph <- mask_idxstats_ratio_treatment_replicate %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% ggplot() +
  aes(x = Treatment, y = rdna_ratio) + 
  geom_point(aes(colour = Replicate), size = 2, alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Treatment condition", y = "rdna_ratio") +
  theme_bw()+
  ggtitle("rDNA ratio - not filtered")
idxstats_rdna_ratio_graph

library(dplyr)
library(forcats)
library(ggplot2)
mask_idxstats_ratio_treatment_replicate %>% mutate(Treatment=fct_relevel(Treatment, 'CON_0')) %>% distinct() %>% ggplot() + aes(x= Treatment , y= rdna_ratio) +
  geom_jitter(aes(color=Replicate), width=0.1, alpha=0.5) +
  ggpubr::stat_compare_means(paired = TRUE, aes(label = paste0("p = ", ..p.format..)), ref.group = 'CON_0',
                             hide.ns=T,
                             bracket.size=1, label.y=0.002) +
  theme_bw() +
  labs(x="Treatment Condition", y="rDNA ratio (rDNA/total reads)") +
  coord_cartesian(ylim=c(0.05, 0.07)) +
  ggtitle("rDNA ratio - not filtered")


#add Treatment condition


mask_idxstats_ratio_treatment_replicate



library(dplyr)


detach("package:methylKit")

mask_idxstats_reads_ratio <- mask_idxstats_ratio_treatment_replicate %>% select(sample, rdna_ratio, Treatment, Replicate) %>% group_by(sample)

mask_idxstats_reads_ratio_unique <- distinct(mask_idxstats_reads_ratio)





mask_idxstats_ratio_unique_treatment <- mask_idxstats_reads_ratio_unique %>%
  mutate(Treatment = case_when(
    endsWith(sample, "0") ~ "CON_0",
    endsWith(sample, "CON") ~ "CON",
    endsWith(sample, "VEH") ~ "VEH",
    endsWith(sample, "DEC") ~ "DEC",
    endsWith(sample, "S") ~ "SDEP",
    endsWith(sample, "S_D") ~ "SDEP_DEC"
  ))

mask_idxstats_ratio_unique_treatment

mask_idxstats_ratio_unique_treatment_replicate <- mask_idxstats_ratio_unique_treatment %>%
  mutate(Replicate = case_when(
    startsWith(sample, "N1") ~ "1",
    startsWith(sample, "N2") ~ "2",
    startsWith(sample, "N3") ~ "3",
    startsWith(sample, "N4") ~ "4"
  ))
mask_idxstats_ratio_unique_treatment_replicate




write.csv(mask_idxstats_ratio_unique_treatment_replicate, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/mask_not_filtered_rDNAratio/mask_idxstats_reads_ratio_nf_tr.csv", row.names = FALSE)


###30_51

##read tsv

#idxstats

fd_idxstats = "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/bam_mask_mouse_30_51/tsv_idxstats/"

library(dplyr)
mask_30_51_idxstats <- bind_rows(lapply(X = list.files(fd_idxstats),
                                  FUN=function(x) {
                                    read_tsv(paste0(fd_idxstats,x),col_names=c("ref_seq", "seq_length", "mapped_read_segments", "unmapped_read_segments")) %>% mutate(sample=x) 
                                  }))

mask_30_51_idxstats_table <- mask_30_51_idxstats %>% group_by(sample) %>% mutate(sample=gsub("_idxstats.tsv","", sample))
mask_30_51_idxstats_table

write.csv(mask_30_51_idxstats_table, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/mask_30_51_rDNAratio/idxstats_30_51.csv", row.names = FALSE)


#calculate ratio betweeen rDNA and total 

mask_30_51_idxstats_ratio <- mask_30_51_idxstats_table %>% group_by(sample) %>% mutate(total_reads = sum(mapped_read_segments)) %>% mutate(rdna_reads = sum(mapped_read_segments[ref_seq == "BK000964.3_looped_3008"])) %>% mutate(rdna_ratio = rdna_reads/total_reads)

mask_30_51_idxstats_ratio

write.csv(mask_30_51_idxstats_ratio, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/mask_30_51_rDNAratio/rdna_ratio_idxstats_30_51.csv", row.names = FALSE)

#Treatment condition and replicate

mask_30_51_idxstats_ratio_treatment <- mask_30_51_idxstats_ratio %>%
  mutate(Treatment = case_when(
    endsWith(sample, "0") ~ "CON_0",
    endsWith(sample, "CON") ~ "CON",
    endsWith(sample, "VEH") ~ "VEH",
    endsWith(sample, "DEC") ~ "DEC",
    endsWith(sample, "S") ~ "SDEP",
    endsWith(sample, "S_D") ~ "SDEP_DEC"
  ))

mask_30_51_idxstats_ratio_treatment

mask_30_51_idxstats_ratio_treatment_replicate <- mask_30_51_idxstats_ratio_treatment %>%
  mutate(Replicate = case_when(
    startsWith(sample, "N1") ~ "1",
    startsWith(sample, "N2") ~ "2",
    startsWith(sample, "N3") ~ "3",
    startsWith(sample, "N4") ~ "4"
  ))
mask_30_51_idxstats_ratio_treatment_replicate


#make table simpler
mask_30_51_ratio <- mask_30_51_idxstats_ratio_treatment_replicate %>% select(sample, rdna_ratio, Treatment, Replicate) %>% group_by(sample)

mask_30_51_ratio_unique <- distinct(mask_30_51_ratio)


write.csv(mask_30_51_ratio_unique, "/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/mask_30_51_ratio_unique.csv", row.names = FALSE)

