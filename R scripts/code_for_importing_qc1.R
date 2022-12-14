#30.7.20
#checked working direcotyr
getwd()
#changed workinign directory via  "session" in menu
setwd("~/Desktop/bismark")
setwd("~/Desktop/bismark/bismark")
#import your csv file to your Global Enviironment
qc <- read.csv("qc1.csv", header = TRUE, sep =",")

qc <- read.csv("qc1.csv", head = TRUE, sep =",")

#alt code, and also for viewing the file 
library(readr)
qc1 <- read_csv("qc1.csv")
View(qc1)

#to avoid annoying message
read_csv('qc1.csv', col_types=cols())

