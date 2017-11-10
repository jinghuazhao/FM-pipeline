# 10-11-2017 MRC-Epid JHZ

options(digits=3, scipen=20, width=200)

library(openxlsx)
bed <- read.table("st.bed", as.is=TRUE, header=TRUE)
for(r in 1:nrow(bed))
{
   rsid <- with(bed, rsid)[r]
   f <- paste(bed[r,1], bed[r,2],bed[r,3],sep="_")
   xlsx <- paste0("chr", f, "-", rsid, ".xlsx")
   unlink(xlsx, recursive = FALSE, force = TRUE)
   wb <- createWorkbook(xlsx)
   if(file.exists(paste0(rsid, "-000001.png")))
   {
     addWorksheet(wb, paste0(rsid, "_plot"))
     insertImage(wb, paste0(rsid, "_plot"), paste0(rsid, "-000001.png"), width=18, height=12)
   }
   snp <- read.table(paste0("chr", f, ".snp"), as.is=TRUE, header=TRUE)
   addWorksheet(wb, paste0(f, ".snp"))
   writeDataTable(wb, paste0(f, ".snp"), snp)
   config <- read.table(paste0("chr", f, ".config"), as.is=TRUE, header=TRUE)
   addWorksheet(wb, paste0(f, ".cfg"))
   writeDataTable(wb, paste0(f, ".cfg"), config)
   chk <- readLines(paste0("chr", f, ".chk"))
   addWorksheet(wb, paste0(f, ".chk"))
   writeDataTable(wb, paste0(f, ".chk", noquote(chk)))
   saveWorkbook(wb, file=xlsx, overwrite=TRUE)
}
