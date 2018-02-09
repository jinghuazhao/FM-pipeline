# 9-2-2018 MRC-Epid JHZ

library(openxlsx)
xlsx <- "gcta-jam.xlsx"
ld <- read.table("ld.txt",as.is=TRUE,fill=TRUE,header=TRUE)
wb <- loadWorkbook(xlsx)
addWorksheet(wb,"LD")
writeDataTable(wb,"LD",ld)
saveWorkbook(wb,file=xlsx,overwrite=TRUE)
