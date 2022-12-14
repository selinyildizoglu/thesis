.libPaths('/data/Blizard-Rakyan/Fran/R/packages_4.1.1/')

if(any(grepl("package:methylKit", search()))) 
  detach("package:methylKit") else message("methylKit not loaded")

detach("package:biomaRt")

library(tidyverse)
fd <- "/data/Blizard-Rakyan/Selin/projects/MEF_rRNA/RNA-Seq/data/lofreq/"
read_snps_lofreq_rRNA <- bind_rows(lapply(X=list.files(fd), FUN = function(x) { 
  read_tsv(paste0(fd, x, "/", x, ".snps.vcf"), 
           comment = "#", col_names = F, col_types=cols()) %>% mutate(sample=x)})) %>% select(sample, pos=X2, ref=X4, alt=X5, info=X8)%>%
  separate(info, into=c('dp', 'af'), sep=';', extra='drop') %>%
  mutate(af=as.numeric(gsub('AF=','',af)), dp = as.numeric(gsub('DP=','',dp))) 
  
read_snps_lofreq_rRNA

write.csv(read_snps_lofreq_rRNA, "/data/Blizard-Rakyan/Selin/projects/MEF_rRNA/RNA-Seq/data/lofreq/read_snps_lofreq_rRNA.csv")


##mouse strains rRNA-seq

library(tidyverse)
fd <- "/data/Blizard-Rakyan/Selin/projects/mouse_strains/RNA-SeqSel/data/lofreq/"
strain_snps_lofreq_rRNA <- bind_rows(lapply(X=list.files(fd), FUN = function(x) { 
  read_tsv(paste0(fd, x, "/", x, ".snps.vcf"), 
           comment = "#", col_names = F, col_types=cols()) %>% mutate(sample=x)})) %>% select(sample, pos=X2, ref=X4, alt=X5, info=X8)%>%
  separate(info, into=c('dp', 'af'), sep=';', extra='drop') %>%
  mutate(af=as.numeric(gsub('AF=','',af)), dp = as.numeric(gsub('DP=','',dp))) 

strain_snps_lofreq_rRNA

write.csv(strain_snps_lofreq_rRNA, "/data/Blizard-Rakyan/Selin/projects/mouse_strains/RNA-SeqSel/csv_tables/strain_snps_lofreq_rRNA.csv", row.names = FALSE)

#many - rRNA-seq

library(tidyverse)
fd <- "/data/Blizard-Rakyan/Selin/projects/mouse_strains/RNA-SeqSel/data/lofreq/"
strain_many_lofreq_rRNA <- bind_rows(lapply(X=list.files(fd), FUN = function(x) { 
  read_tsv(paste0(fd, x, "/", x, ".many.vcf"), 
           comment = "#", col_names = F, col_types=cols()) %>% mutate(sample=x)})) %>% select(sample, pos=X2, ref=X4, alt=X5, info=X8)%>%
  separate(info, into=c('dp', 'af'), sep=';', extra='drop') %>%
  mutate(af=as.numeric(gsub('AF=','',af)), dp = as.numeric(gsub('DP=','',dp))) 

strain_many_lofreq_rRNA

write.csv(strain_many_lofreq_rRNA, "/data/Blizard-Rakyan/Selin/projects/mouse_strains/RNA-SeqSel/csv_tables/strain_many_lofreq_rRNA.csv", row.names = FALSE)




##Cut & Tag

library(tidyverse)
fd <- "/data/Blizard-Rakyan/Selin/projects/mouse_strains/Cut_and_Tag_thesis/data/lofreq_indels/"
cuttag_snps_lofreq <- bind_rows(lapply(X=list.files(fd), FUN = function(x) { 
  read_tsv(paste0(fd, x, "/", x, ".snps.vcf"), 
           comment = "#", col_names = F, col_types=cols()) %>% mutate(sample=x)})) %>% select(sample, pos=X2, ref=X4, alt=X5, info=X8)%>%
  separate(info, into=c('dp', 'af'), sep=';', extra='drop') %>%
  mutate(af=as.numeric(gsub('AF=','',af)), dp = as.numeric(gsub('DP=','',dp))) 

cuttag_snps_lofreq

write.csv(cuttag_snps_lofreq, "/data/Blizard-Rakyan/Selin/projects/mouse_strains/Cut_and_Tag_thesis/csv_tables/cuttag_snps_lofreq.csv", row.names = FALSE)


