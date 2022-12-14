#6.8.20

#WGS data

qcWGS <- inner_join(qcWGS, sampleID_mouseID, by = "sample")

qcWGS %>% filter(sample == "S1")

S1 <- qcWGS %>%  filter(sample == "S1")
AverageS1 <- mean(qcWGS  %>% slice(1:2))

