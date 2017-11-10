# 20-11-2017 MRC-Epid JHZ

options(digits=3, scipen=20, width=200)

library(openxlsx)
p <- Sys.getenv("p")
xlsx <- paste0("finemap-top", p, ".xlsx")
wb <- createWorkbook(xlsx)
snp <- read.table(paste0("snp", p, ".dat"), as.is=TRUE, header=TRUE, sep=" ")
addWorksheet(wb, "snp")
writeDataTable(wb, "snp", snp)
config <- read.table("config", p, ".dat"), as.is=TRUE, header=TRUE, sep=" ")
addWorksheet(wb, "config")
writeDataTable(wb, "config", subset(config,config_prob>0.4))
saveWorkbook(wb, file=xlsx, overwrite=TRUE)
