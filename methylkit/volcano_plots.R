#volcano plots

#p + theme(text=element_text(size=20), #change font size of all text
#axis.text=element_text(size=20), #change font size of axis text
#axis.title=element_text(size=20), #change font size of axis titles
#plot.title=element_text(size=20), #change font size of plot title
#legend.text=element_text(size=20), #change font size of legend text
#legend.title=element_text(size=20)) #change font size of legend title 


#... %>% mutate(fdr=p.adjust(pvalue, method='fdr'))


library(readr)
c_c_0_mask_no_pc2 <- read_csv("MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/dmr_mask/c_c_0_mask_no_pc2.csv")
c_c_0_mask_no_pc2

rDNA <- c_c_0_mask_no_pc2 %>% filter(seqnames == 'BK000964.3_looped_3008')

#control day 5 - control day 0 
library(ggplot2)
con5_0_no_pc2_volcano_rdna <-c_c_0_mask_no_pc2 %>% filter(meth.diff > 5 | meth.diff < -5) %>% ggplot() +
  aes(x = meth.diff, y=-log10(pvalue))+
  geom_point(data = rDNA, aes(x = meth.diff, y=-log10(pvalue)), colour = 'red', alpha = 1) +
  theme_bw()+
  labs(x = "Methylation difference", y = "-log10(pvalue)") +
  theme(legend.position='none') +
  ggtitle("DMRs:Control day 5 - Control Day 0")
con5_0_no_pc2_volcano_rdna

library(ggplot2)
con5_0_no_pc2_volcano <-c_c_0_mask_no_pc2 %>% filter(meth.diff > 5 | meth.diff < -5) %>% ggplot() +
  aes(x = meth.diff, y=-log10(pvalue))+
  geom_point(aes(color = seqnames), alpha = 1)  +
  theme_bw()+
  labs(x = "Methylation difference", y = "-log10(pvalue)") +
  ggtitle("DMRs:Control day 5 - Control Day 0")
con5_0_no_pc2_volcano


#13.7.22

c_c_0_mask_no_pc2 %>% mutate(fdr=p.adjust(pvalue, method='fdr'))

#geom_hline(yintercept=-log10(0.01), col="red")


#mutate(fdr=...,sig=(fdr < 0.01 & meth.diff > 5)) %>%ggplot() + aes(x = ..., y = ..., colour = sig) + ...scale_colour_manual(values=c('FALSE'='black', 'TRUE'='red') +

library(ggplot2)
con5_0_no_pc2_volcano <-c_c_0_mask_no_pc2 %>% mutate(fdr=p.adjust(pvalue, method='fdr'), sig=(fdr < 0.01 & meth.diff > 5 & meth.diff < -5)) %>% ggplot() +
  aes(x = meth.diff, y= -log10(fdr), col = sig) +
  scale_colour_manual(values=c('FALSE'='black', 'TRUE'='red')) + 
  geom_point(alpha = 1)  +
  theme_bw()+
  geom_vline(xintercept=c(-5, 5), col="grey") +
  geom_hline(yintercept=-log10(0.01), col="grey") +
  labs(x = "Methylation difference", y = "-log10 (FDR)") +
  ggtitle("DMRs:Control day 5 - Control Day 0")
con5_0_no_pc2_volcano


library(readr)
s_c_0_mask_no_pc2 <- read_csv("MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/dmr_mask/s_c_0_mask_no_pc2.csv")
s_c_0_mask_no_pc2

#serum deprivation - control day 0
library(ggplot2)
s5_0_no_pc2_volcano <-s_c_0_mask_no_pc2 %>% filter(meth.diff > 5 | meth.diff < -5) %>% ggplot() +
  aes(x = meth.diff, y=-log10(pvalue))+
  geom_point(aes(color = seqnames), alpha = 1)+
  theme_bw()+
  labs(x = "Methylation difference", y = "-log10(pvalue)") +
  ggtitle("DMRs: Serum deprivation - Control Day 0")
s5_0_no_pc2_volcano

s_c_0_rDNA <- s_c_0_mask_no_pc2 %>% filter(seqnames == 'BK000964.3_looped_3008')

write.csv(s_c_0_rDNA, "MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/DMR/s_c_0_rDNA.csv", row.names = FALSE)




library(readr)
s_c_mask_no_pc2 <- read_csv("MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/dmr_mask/s_c_mask_no_pc2.csv")
s_c_mask_no_pc2


rDNA <- s_c_mask_no_pc2 %>% filter(meth.diff > 5 | meth.diff < -5) %>% filter(seqnames == 'BK000964.3_looped_3008')


#serum dep - control 

write.csv(s_c_0_rDNA, "MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/DMR/s_c_0_rDNA.csv", row.names = FALSE)

library(ggplot2)
s_c_volcano_no_pc2 <- s_c_mask_no_pc2 %>%  filter(meth.diff > 5 | meth.diff < -5) %>% ggplot() +
aes(x = meth.diff, y=-log10(pvalue))+
geom_point()+
theme_bw()+
labs(x = "Methylation difference", y = "-log10(pvalue)") +
ggtitle("DMRs: Serum deprivation - Control day 5")
s_c_volcano_no_pc2

coloured_coded_s_c <- s_c_volcano_no_pc2 + geom_point(data = rDNA, aes(x = meth.diff, y=-log10(pvalue)), colour = 'red', alpha = 1)
coloured_coded_s_c

s_c_rDNA <- s_c_mask_no_pc2 %>% filter(meth.diff > 5 | meth.diff < -5) %>% filter(seqnames == 'BK000964.3_looped_3008')

write.csv(s_c_rDNA, "MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/csv_tables/DMR/s_c_rDNA.csv", row.names = FALSE)

#vehicle - control

library(readr)
v_c_mask_no_pc2 <- read_csv("MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/dmr_mask/v_c_mask_no_pc2.csv")
v_c_mask_no_pc2

rDNA <- v_c_mask_no_pc2 %>% filter(meth.diff > 5 | meth.diff < -5) %>% filter(seqnames == 'BK000964.3_looped_3008')

library(ggplot2)
v_c5_no_pc2_volcano <- v_c_mask_no_pc2 %>% filter(meth.diff > 5 | meth.diff < -5) %>% ggplot() +
  aes(x = meth.diff, y=-log10(pvalue))+
  geom_point(aes(color = seqnames), alpha = 1)+
  theme_bw()+
  labs(x = "Methylation difference", y = "-log10(pvalue)") +
  ggtitle("DMRs: Vehicle - Control day 5")
v_c5_no_pc2_volcano

coloured_coded_v_c <- v_c5_no_pc2_volcano + geom_point(data = rDNA, aes(x = meth.diff, y=-log10(pvalue)), colour = 'red', alpha = 1)

#decitabine - vehicle 

#rDNA <- v_c_mask_no_pc2 %>% filter(meth.diff > 5 | meth.diff < -5) %>% filter(seqnames == 'BK000964.3_looped_3008')

rDNA <- d_v_mask_no_pc2 %>% filter(meth.diff > 5 | meth.diff < -5) %>% filter(seqnames == 'BK000964.3_looped_3008')

library(readr)
d_v_mask_no_pc2 <- read_csv("MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/dmr_mask/d_v_mask_no_pc2.csv")
d_v_mask_no_pc2

library(ggplot2)
d_v_no_pc2_volcano <- d_v_mask_no_pc2 %>% filter(meth.diff > 5 | meth.diff < -5) %>% ggplot() +
  aes(x = meth.diff, y=-log10(pvalue))+
  geom_point()+
  theme_bw()+
  labs(x = "Methylation difference", y = "-log10(pvalue)") +
  ggtitle("DMRs: Decitabine - Vehicle")
d_v_no_pc2_volcano 

dv_rdna <- d_v_no_pc2_volcano + geom_point(data = rDNA, aes(x = meth.diff, y=-log10(pvalue)), colour = 'red', alpha = 1)
dv_rdna

#combined serum deprivation and decitabine + control

library(readr)
s_d_c_mask_no_pc2 <- read_csv("MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/dmr_mask/s_d_c_mask_no_pc2.csv")
s_d_c_mask_no_pc2

rDNA <- s_d_c_mask_no_pc2 %>% filter(meth.diff > 5 | meth.diff < -5) %>% filter(seqnames == 'BK000964.3_looped_3008')

library(ggplot2)
sd_c_no_pc2_volcano <- s_d_c_mask_no_pc2 %>% filter(meth.diff > 5 | meth.diff < -5) %>% ggplot() +
  aes(x = meth.diff, y=-log10(pvalue))+
  geom_point(aes(color = seqnames), alpha = 1)+
  theme_bw()+
  labs(x = "Methylation difference", y = "-log10(pvalue)") +
  ylim(0,1000) +
  ggtitle("DMRs: combined serum deprivation + decitabine - Control day 5")
sd_c_no_pc2_volcano

sdc_rdna <- sd_c_no_pc2_volcano + geom_point(data = rDNA, aes(x = meth.diff, y=-log10(pvalue)), colour = 'red', alpha = 1)
sdc_rdna

#combined S_D - D

#filter(meth.diff > 5 | meth.diff < -5) %>% 

library(readr)
s_d_d_mask_no_pc2 <- read_csv("MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/dmr_mask/s_d_d_mask_no_pc2.csv")
s_d_d_mask_no_pc2

rDNA <- s_d_d_mask_no_pc2 %>% filter(meth.diff > 5 | meth.diff < -5) %>% filter(seqnames == 'BK000964.3_looped_3008')


library(ggplot2)
s_d_d_no_pc2_volcano <- s_d_d_mask_no_pc2 %>% filter(meth.diff > 5 | meth.diff < -5) %>% ggplot() +
  aes(x = meth.diff, y=-log10(pvalue))+
  geom_point()+
  theme_bw()+
  labs(x = "Methylation difference", y = "-log10(pvalue)") +
  ggtitle("DMRs: combined serum deprivation & decitabine - decitabine only") +
  theme(legend.position='none')
s_d_d_no_pc2_volcano

sdd_rdna <- s_d_d_no_pc2_volcano + geom_point(data = rDNA, aes(x = meth.diff, y=-log10(pvalue)), colour = 'red', alpha = 1)
sdd_rdna

#decitabine and day 0

library(readr)
d_c_0_mask_no_pc2 <- read_csv("MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/dmr_mask/d_c_0_mask_no_pc2.csv")
d_c_0_mask_no_pc2

rDNA <- d_c_0_mask_no_pc2 %>% filter(meth.diff > 5 | meth.diff < -5) %>% filter(seqnames == 'BK000964.3_looped_3008')

library(ggplot2)
d_c_0_no_pc2_volcano <- d_c_0_mask_no_pc2 %>% filter(meth.diff > 5 | meth.diff < -5) %>% ggplot() +
  aes(x = meth.diff, y=-log10(pvalue))+
  geom_point()+
  theme_bw()+
  labs(x = "Methylation difference", y = "-log10(pvalue)") +
  ggtitle("DMRs: decitabine - control day 0") +
  theme(legend.position='none')
d_c_0_no_pc2_volcano

dc0_rdna <- d_c_0_no_pc2_volcano + geom_point(data = rDNA, aes(x = meth.diff, y=-log10(pvalue)), colour = 'red', alpha = 1)
dc0_rdna

#no filter


rDNA <- d_c_0_mask_no_pc2 %>% filter(seqnames == 'BK000964.3_looped_3008')

library(ggplot2)
d_c_0_no_pc2_volcano <- d_c_0_mask_no_pc2 %>%  ggplot() +
  aes(x = meth.diff, y=-log10(pvalue))+
  geom_point()+
  theme_bw()+
  labs(x = "Methylation difference", y = "-log10(pvalue)") +
  ggtitle("DMRs: decitabine - control day 0") +
  theme(legend.position='none')
d_c_0_no_pc2_volcano

dc0_rdna <- d_c_0_no_pc2_volcano + geom_point(data = rDNA, aes(x = meth.diff, y=-log10(pvalue)), colour = 'red', alpha = 1)
dc0_rdna


###################
###########13.7.22

library(ggplot2)
con5_0_no_pc2_volcano <-c_c_0_mask_no_pc2 %>% mutate(fdr=p.adjust(pvalue, method='fdr'), sig=(fdr < 0.01 & meth.diff > 5 & meth.diff < -5)) %>% ggplot() +
  aes(x = meth.diff, y= -log10(fdr), col = sig) +
  scale_colour_manual(values=c('FALSE'='black', 'TRUE'='red')) + 
  geom_point(alpha = 1)  +
  theme_bw()+
  geom_vline(xintercept=c(-5, 5), col="grey") +
  geom_hline(yintercept=-log10(0.01), col="grey") +
  labs(x = "Methylation difference", y = "-log10 (FDR)") +
  ggtitle("DMRs:Control day 5 - Control Day 0")
con5_0_no_pc2_volcano


#s_c_0

s_c_0_mask_no_pc2 %>% mutate(fdr=p.adjust(pvalue, method='fdr'), sig=(fdr < 0.01 & abs(meth.diff) > 5)) %>% filter(sig) %>% arrange(desc(fdr))


s_c0_rDNA <- s_c_0_mask_no_pc2 %>% mutate(fdr=p.adjust(pvalue, method='fdr')) %>% filter(seqnames == 'BK000964.3_looped_3008')

library(ggplot2)
s5_0_no_pc2_volcano <-s_c_0_mask_no_pc2 %>% mutate(fdr=p.adjust(pvalue, method='fdr'), sig=(fdr < 0.01 & abs(meth.diff) > 5)) %>% ggplot() +
  aes(x = meth.diff, y= -log10(fdr), colour = factor(sig)) +
  scale_colour_manual(values=c('FALSE'='black', 'TRUE'='red')) + 
  geom_point(alpha = 1, size = 1)+
  theme_bw()+
  geom_vline(xintercept=c(-5, 5), col="grey") +
  geom_hline(yintercept=-log10(0.01), col="grey") +
  labs(x = "Methylation difference", y = "-log10(FDR)") +
  ggtitle("DMRs: Serum deprivation - Control Day 0")
s5_0_no_pc2_volcano


coloured_coded_s_c0 <- s5_0_no_pc2_volcano + geom_point(data = s_c0_rDNA, aes(x = meth.diff, y=-log10(fdr)), colour = 'blue', alpha = 1)
coloured_coded_s_c0


mutate(type = case_when(seqnames == "BK..." ~ 'rDNA', fdr < 0.01 & abs(meth.diff) > 5 ~ 'Significant', TRUE ~ 'Non-Significant'))



##############13.7.22 improved 

#con 5 - con 0 

library(ggplot2)
con5_0_no_pc2_volcano <-c_c_0_mask_no_pc2 %>% mutate(fdr=p.adjust(pvalue, method='fdr')) %>% mutate(Type = case_when(seqnames == "BK000964.3_looped_3008" ~ 'rDNA', fdr < 0.01 & abs(meth.diff) > 5 ~ 'Significant', TRUE ~ 'Non-Significant')) %>% ggplot() +
  aes(x = meth.diff, y= -log10(fdr), colour = Type) +
  scale_colour_manual(values=c('rDNA'='blue', 'Significant'='red', 'Non-significant' = 'black')) + 
  geom_point(alpha = 1, size = 2)  +
  theme_bw()+
  geom_vline(xintercept=c(-5, 5), col="grey") +
  geom_hline(yintercept=-log10(0.01), col="grey") +
  labs(x = "Methylation difference", y = "-log10 (FDR)") +
  scale_x_continuous(limits=c(-40,40), breaks=seq(-40,40,10)) + 
  ggtitle("Control Day 5 - Control Day 0") +
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(), 
        plot.title=element_text(size=20, hjust = 0.5), axis.title=element_text(size=10))
con5_0_no_pc2_volcano


#S_DEP - CON_0 

library(ggplot2)
s5_0_no_pc2_volcano <-s_c_0_mask_no_pc2 %>% mutate(fdr=p.adjust(pvalue, method='fdr')) %>% mutate(Type = case_when(seqnames == "BK000964.3_looped_3008" ~ 'rDNA', fdr < 0.01 & abs(meth.diff) > 5 ~ 'Significant', TRUE ~ 'Non-Significant')) %>% 
  ggplot() +
  aes(x = meth.diff, y= -log10(fdr), colour = Type) +
  scale_colour_manual(values=c('rDNA'='blue', 'Significant'='red', 'Non-significant' = 'black')) + 
  geom_point(alpha = 1, size = 2)+
  theme_bw()+
  geom_vline(xintercept=c(-5, 5), col="grey") +
  geom_hline(yintercept=-log10(0.01), col="grey") +
  labs(x = "Methylation difference", y = "-log10(FDR)") +
  scale_x_continuous(limits=c(-40,40), breaks=seq(-40,40,10)) + 
  ggtitle("Serum deprivation - Control Day 0") +
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(), 
        plot.title=element_text(size=17, hjust = 0.5), axis.title=element_text(size=10))
s5_0_no_pc2_volcano


#S_DEP - CON

library(ggplot2)
s_c_volcano_no_pc2 <- s_c_mask_no_pc2 %>%  mutate(fdr=p.adjust(pvalue, method='fdr')) %>% 
  mutate(Type = case_when(seqnames == "BK000964.3_looped_3008" ~ 'rDNA', fdr < 0.01 & abs(meth.diff) > 5 ~ 'Significant', TRUE ~ 'Non-Significant')) %>% 
  ggplot() +
  aes(x = meth.diff, y=-log10(fdr), colour = Type)+
  scale_colour_manual(values=c('rDNA'='blue', 'Significant'='red', 'Non-significant' = 'black')) + 
  geom_point(alpha = 1, size = 2)+
  theme_bw()+
  geom_vline(xintercept=c(-5, 5), col="grey") +
  geom_hline(yintercept=-log10(0.01), col="grey") +
  labs(x = "Methylation difference", y = "-log10(FDR)") +
  scale_x_continuous(limits=c(-40,40), breaks=seq(-40,40,10)) + 
  ggtitle("Serum deprivation - Control Day 5") +
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(), 
        plot.title=element_text(size=17, hjust = 0.5), axis.title=element_text(size=10))
s_c_volcano_no_pc2

# vehicle - control

library(ggplot2)
v_c5_no_pc2_volcano <- v_c_mask_no_pc2 %>% mutate(fdr=p.adjust(pvalue, method='fdr')) %>% 
  mutate(Type = case_when(seqnames == "BK000964.3_looped_3008" ~ 'rDNA', fdr < 0.01 & abs(meth.diff) > 5 ~ 'Significant', TRUE ~ 'Non-Significant')) %>%
  ggplot() +
  aes(x = meth.diff, y=-log10(fdr), colour = Type) +
  scale_colour_manual(values=c('rDNA'='blue', 'Significant'='red', 'Non-significant' = 'black')) + 
  geom_point(alpha = 1, size = 2)+
  theme_bw()+
  geom_vline(xintercept=c(-5, 5), col="grey") +
  geom_hline(yintercept=-log10(0.01), col="grey") +
  labs(x = "Methylation difference", y = "-log10(FDR)") +
  scale_x_continuous(limits=c(-40,40), breaks=seq(-40,40,10)) + 
  scale_y_continuous(limits=c(0,300), breaks=seq(0,300,100)) + 
  ggtitle("Vehicle - Control Day 5") +
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(), 
        plot.title=element_text(size=20, hjust = 0.5), axis.title=element_text(size=10))
v_c5_no_pc2_volcano


# decitabine - vehicle

library(ggplot2)
d_v_no_pc2_volcano <- d_v_mask_no_pc2 %>% mutate(fdr=p.adjust(pvalue, method='fdr')) %>% 
  mutate(Type = case_when(seqnames == "BK000964.3_looped_3008" ~ 'rDNA', fdr < 0.01 & abs(meth.diff) > 5 ~ 'Significant', TRUE ~ 'Non-Significant')) %>%
  ggplot() +
  aes(x = meth.diff, y=-log10(fdr), colour = Type) +
  scale_colour_manual(values=c('rDNA'='blue', 'Significant'='red', 'Non-significant' = 'black')) + 
  geom_point(alpha = 1, size = 2)+
  theme_bw()+
  geom_vline(xintercept=c(-5, 5), col="grey") +
  geom_hline(yintercept=-log10(0.01), col="grey") +
  labs(x = "Methylation difference", y = "-log10(FDR)") +
  scale_x_continuous(limits=c(-40,40), breaks=seq(-40,40,10)) + 
  ggtitle("Decitabine - Vehicle") +
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(), 
        plot.title=element_text(size=20, hjust = 0.5), axis.title=element_text(size=10))
d_v_no_pc2_volcano 

# decitabine - control

library(readr)
d_c_mask_no_pc2 <- read_csv("MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/dmr_mask/d_c_mask_no_pc2.csv")
d_c_mask_no_pc2


library(ggplot2)
d_c_no_pc2_volcano <- d_c_mask_no_pc2 %>%  mutate(fdr=p.adjust(pvalue, method='fdr')) %>% 
  mutate(Type = case_when(seqnames == "BK000964.3_looped_3008" ~ 'rDNA', fdr < 0.01 & abs(meth.diff) > 5 ~ 'Significant', TRUE ~ 'Non-Significant')) %>%
  ggplot() +
  aes(x = meth.diff, y=-log10(fdr), colour = Type) +
  scale_colour_manual(values=c('rDNA'='blue', 'Significant'='red', 'Non-significant' = 'black')) + 
  geom_point(alpha = 1, size = 2) +
  theme_bw()+
  geom_vline(xintercept=c(-5, 5), col="grey") +
  geom_hline(yintercept=-log10(0.01), col="grey") +
  labs(x = "Methylation difference", y = "-log10(FDR)") +
  scale_x_continuous(limits=c(-40,40), breaks=seq(-40,40,10)) + 
  ggtitle("Decitabine - Control Day 5") +
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(), 
        plot.title=element_text(size=20, hjust = 0.5), axis.title=element_text(size=10))
d_c_no_pc2_volcano




# decitabine - control day 0

library(ggplot2)
d_c_0_no_pc2_volcano <- d_c_0_mask_no_pc2 %>%  mutate(fdr=p.adjust(pvalue, method='fdr')) %>% 
  mutate(Type = case_when(seqnames == "BK000964.3_looped_3008" ~ 'rDNA', fdr < 0.01 & abs(meth.diff) > 5 ~ 'Significant', TRUE ~ 'Non-Significant')) %>%
  ggplot() +
  aes(x = meth.diff, y=-log10(fdr), colour = Type) +
  scale_colour_manual(values=c('rDNA'='blue', 'Significant'='red', 'Non-significant' = 'black')) + 
  geom_point(alpha = 1, size = 2)+
  theme_bw()+
  geom_vline(xintercept=c(-5, 5), col="grey") +
  geom_hline(yintercept=-log10(0.01), col="grey") +
  labs(x = "Methylation difference", y = "-log10(FDR)") +
  ggtitle("Decitabine - Control Day 0") +
  scale_x_continuous(limits=c(-40,40), breaks=seq(-40,40,10)) + 
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(), 
        plot.title=element_text(size=20, hjust = 0.5), axis.title=element_text(size=10))
d_c_0_no_pc2_volcano




# S_DEP_DEC - DEC

library(ggplot2)
s_d_d_no_pc2_volcano <- s_d_d_mask_no_pc2 %>%  mutate(fdr=p.adjust(pvalue, method='fdr')) %>% 
  mutate(Type = case_when(seqnames == "BK000964.3_looped_3008" ~ 'rDNA', fdr < 0.01 & abs(meth.diff) > 5 ~ 'Significant', TRUE ~ 'Non-Significant')) %>%
  ggplot() +
  aes(x = meth.diff, y=-log10(fdr), colour = Type) +
  scale_colour_manual(values=c('rDNA'='blue', 'Significant'='red', 'Non-significant' = 'black')) + 
  geom_point(alpha = 1, size = 2)+
  theme_bw()+
  geom_vline(xintercept=c(-5, 5), col="grey") +
  geom_hline(yintercept=-log10(0.01), col="grey") +
  labs(x = "Methylation difference", y = "-log10(FDR)") +
  scale_x_continuous(limits=c(-40,40), breaks=seq(-40,40,10)) + 
  ggtitle("Serum deprivation & Decitabine 
          - Decitabine") +
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(), 
        plot.title=element_text(size=19, hjust = 0.5), axis.title=element_text(size=10))
s_d_d_no_pc2_volcano

#S_DEP_DEC - CON

library(ggplot2)
sd_c_no_pc2_volcano <- s_d_c_mask_no_pc2 %>% mutate(fdr=p.adjust(pvalue, method='fdr')) %>% 
  mutate(Type = case_when(seqnames == "BK000964.3_looped_3008" ~ 'rDNA', fdr < 0.01 & abs(meth.diff) > 5 ~ 'Significant', TRUE ~ 'Non-Significant')) %>%
  ggplot() +
  aes(x = meth.diff, y=-log10(fdr), colour = Type) +
  scale_colour_manual(values=c('rDNA'='blue', 'Significant'='red', 'Non-significant' = 'black')) + 
  geom_point(alpha = 1, size = 2)+
  theme_bw()+
  geom_vline(xintercept=c(-5, 5), col="grey") +
  geom_hline(yintercept=-log10(0.01), col="grey") +
  labs(x = "Methylation difference", y = "-log10(FDR)") +
  scale_x_continuous(limits=c(-40,40), breaks=seq(-40,40,10)) + 
  ggtitle("Serum deprivation & Decitabine 
          - Control Day 5") +
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank(),
        plot.title=element_text(size=19, hjust = 0.5), axis.title=element_text(size=10))
sd_c_no_pc2_volcano



