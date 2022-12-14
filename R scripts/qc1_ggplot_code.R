#30.7.20

#17.4.22


library(tidyverse)

qc1 %>% ggplot() +
  geom_point(
    aes(x = `Mouse ID`,
        y = reads_mappingefficiency),
    size = 2, color = "blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="Mouse ID", y="Mapping Efficiency (%)") +
  ggtitle("Bismark Alignment Rate")

`Mouse ID`
#have to save the inner_join action so you update qc1. 
qc1 <- inner_join(qc1, sampleID_mouseID, by = "sample")
qc1 %>% ggplot() +
  geom_point(
    aes(x = `Mouse ID`,
                    y = reads_mappingefficiency),
         size = 2, color = "blue") +
     theme_bw()+
     theme(axis.text.x = element_text(angle=45, hjust=1)) +
     labs(x="Mouse ID", y="Mapping Efficiency (%)") +
     ggtitle("Bismark Alignment Rate")

qc1 <- inner_join(qc1, sampleID_mouseID, by = "sample")
qc1 %>% ggplot() +
   geom_point(
    aes(x = `Mouse ID`,
             y = reads_1),
    size = 2, color = "blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="Mouse ID", y="Number of Reads") +
  ggtitle("Number of reads per sample")

reads1 <- qc1$reads_1
mean(reads1)
206666847
sd(reads1)
 32658707

SD <- sd(reads1)
avg <- mean(reads1)
SD/avg

qc1$`Mouse ID`[which.max(qc1$reads_1)]
qc1$`Mouse ID`[which.min(qc1$reads_1)]

qc1 %>%
  filter(reads_1==max(reads_1) | reads_1==min(reads_1)) %>%
  select(`Mouse ID`, reads_1) %>%
  arrange(reads_1)

qc1 %>%
  filter(reads_1==max(reads_1) | reads_1==min(reads_1)) %>%
  select(`Mouse ID`, reads_1) %>%
  arrange(reads_1)

#31.7.20

qc <- inner_join(qc, sampleID_mouseID, by = "sample")

qc <- inner_join(qc1, sampleID_mouseID, by = "sample")
qc %>% ggplot() +
  geom_point(
    aes(x = `Mouse ID.x`,
        y = reads_1),
    size = 2, color = "blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="Mouse ID", y="Number of Reads") +
  ggtitle("Number of reads per sample")


qc <- inner_join(qc, sampleID_mouseID, by = "sample")
qc %>% ggplot() +
  geom_point(
    aes(x = `Mouse ID`,
        y = reads_1),
    size = 2, color = "blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="Mouse ID", y="Number of Reads") +
  ggtitle("Number of reads per sample")

qc <- inner_join(qc, sampleID_mouseID, by = "sample")
qc %>% ggplot() +
  geom_point(
    aes(x = `Mouse ID`,
        y = reads_mappingefficiency),
    size = 2, color = "blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="Mouse ID", y="Mapping Efficiency (%)") +
  ggtitle("Bismark Alignment Rate")

qc %>% ggplot() +
  geom_point(
    aes(x = mean_length_trimmed_1,
        y = mean_length_trimmed_2),
    size = 2, color = "blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="Mean Length Trimmed (1)", y="Mean Length Trimmed (2)") +
  ggtitle("Mean Length Trimmed")

qc %>% ggplot() +
  geom_point(
    aes(x = `Mouse ID`,
        y = mean_length_trimmed_1),
    size = 2, color = "blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="Mouse ID", y="Mean Length Trimmed") +
  ggtitle("Average Read Lengths")

qc %>% ggplot() +
  geom_point(
    aes(x = `Mouse ID`,
        y = total_number_of_Cs_analysed),
    size = 2, color = "blue") +
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="Mouse ID", y="Total number of C's analysed") +
  ggtitle("Number of C's analysed in Bismark")


C_methylated_in_CpG_percentage
"total_number_of_C's_analysed"
"total_methylated_C's_in_CpG"
"total_unmethylated_C's_in_CpG"
"C_methylated_in_CpG_percentage"

#need to remove apostrophe


298458/(298458+4484665)

4484665

#adding new column

#Simply create a variable name for the new column and pass in a calculation formula as its value if, for example, you want a new column that's the sum of two existing columns:
dataFrame$newColumn <- dataFrame$oldColumn1 + dataFrame$oldColumn2
#As you can probably guess, this creates a new column called "newColumn" with the sum of oldColumn1 + oldColumn2 in each row.
#For our sample data frame called data, we could add a column for profit margin by dividing profit by revenue and then multiplying by 100:
companiesData$margin <- (companiesData$profit / companiesData$revenue) * 100


qc$percentage_of_C_meth__in_unknown <- qc$Total_C_methylated_in_unknown/(qc$Total_C_methylated_in_unknown + qc$Total_C_unmethylated_in_unknown)*100

qc$bisulfite_conversion_efficiency <- 100 - qc$percentage_of_C_meth__in_unknown

qc %>% ggplot() +
  geom_point(
aes(x = `Mouse ID`,
y = bisulfite_conversion_efficiency),
size = 2, color = "blue") +
theme_bw()+
theme(axis.text.x = element_text(angle=45, hjust=1)) +
labs(x="Mouse ID", y="Estimate of Bisulfite Conversion Efficiency (%)") +
ggtitle("Estimate of Bisulfite Conversion Rate")


#6.8.20

#WGS data

qcWGS <- inner_join(qcWGS, sampleID_mouseID, by = "sample")


#12.8.20

#WGBS complete bismark

qccompletebismark <- inner_join(qccompletebismark, sampleID_mouseID, by = "sample")

qccompletebismark %>% ggplot() +
  geom_point(
aes(x = `Mouse ID.x`,
y = reads_mappingefficiency),
size = 2, color = "blue") +
theme_bw()+
theme(axis.text.x = element_text(angle=45, hjust=1)) +
labs(x="Mouse ID", y="Mapping Efficiency (%)") +
ggtitle("Bismark Alignment Rate")

qccompletebismark %>% ggplot() +
  geom_point(
aes(x = `Mouse ID.x`,
y = C_methylated_in_CpG_percentage),
size = 2, color = "blue") +
theme_bw()+
theme(axis.text.x = element_text(angle=45, hjust=1)) +
labs(x="Mouse ID", y="C methylated in CpG (%)") +
ggtitle("C methylated in CpG (%)")

qccompletebismark %>% ggplot() +
  geom_point(
aes(x = `Mouse ID.x`,
y = total_methylated_Cs_in_CpG),
size = 2, color = "blue") +
theme_bw()+
theme(axis.text.x = element_text(angle=45, hjust=1)) +
labs(x="Mouse ID", y="Total methylated C's in CpG") +
ggtitle("Total Methylated C's in CpG")
