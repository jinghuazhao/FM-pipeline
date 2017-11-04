# 4-11-2017 MRC-Epid JHZ

dir <- Sys.getenv("dir")
setwd(dir)
options(digits=3, scipen=20, width=200)

library(openxlsx)
bed <- read.table("st.bed", as.is=TRUE, header=TRUE)
for(r in nrows(bed))
{
   rsid <- with(bed, rsid)[r]
   f <- paste(bed[r,1:3], sep="_")
   xlsx <- paste0(dir, "/", rsid, ".xlsx")
   unlink(xlsx, recursive = FALSE, force = FALSE)
   wb <- createWorkbook(xlsx)
   if(file.exists(paste0(dir, "/", rsid, ".pdf")))
   {
     addWorksheet(wb, paste0(rsid, "_plot"))
     insertImage(wb, paste0(rsid, "_plot"), paste0(dir, "/", rsid, "-000001.png"), width=18, height=12)
   }
   snp <- read.table(paste0(dir, "/chr", f, ".snp"), as.is=TRUE, header=TRUE)
   addWorksheet(wb, paste0(f, ".snp"))
   writeDataTable(wb, paste0(f, ".snp"), snp)
   config <- read.table(paste0(dir, "/chr", f, ".config"), as.is=TRUE, header=TRUE)
   addWorksheet(wb, paste0(f, ".cfg"))
   writeDataTable(wb, paste0(f, ".cfg"), config)
   saveWorkbook(wb, file=xlsx, overwrite=TRUE)
}
