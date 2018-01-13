# 13-1-2018 MRC-Epid JHZ

options(digits=3, scipen=20, width=200)

library(openxlsx)
bed <- read.table("st.bed", as.is=TRUE, header=TRUE)
xlsx <- paste0("lz.xlsx")
unlink(xlsx, recursive = FALSE, force = TRUE)
wb <- createWorkbook(xlsx)
for(r in 1:nrow(bed))
{
   rsid <- with(bed, rsid)[r]
   f <- paste(bed[r,1], bed[r,2],bed[r,3],sep="_")
   if(file.exists(paste0(rsid, "-000001.png")))
   {
     addWorksheet(wb, paste0(rsid, "_plot"))
     insertImage(wb, paste0(rsid, "_plot"), paste0(rsid, "-000001.png"), width=18, height=12)
   }
}
saveWorkbook(wb, file=xlsx, overwrite=TRUE)

