# 20-9-2017 MRC-Epid JHZ

options(digits=3, scipen=20, width=200)
setwd("/genetics/data/gwas/6-7-17")
library(openxlsx)

region <- read.table("/genetics/data/gwas/6-7-17/doc/region.dat",as.is=TRUE,col.names=c("rsid","chr","pos"),sep="\t")
bed <- read.table("/genetics/data/gwas/6-7-17/doc/st.bed",col.names=c("chr","start","end"))
for(r in 1:97)
{
   rsid <- with(region,rsid)[r]
   f <- paste(bed[r,1],bed[r,2],bed[r,3],sep="_")
   xlsx <- paste0("HRC/",rsid,".xlsx")
   unlink(xlsx, recursive = FALSE, force = FALSE)
   wb <- createWorkbook(xlsx)
   if(file.exists(paste0("HRC/",rsid,".pdf")))
   {
     addWorksheet(wb, paste0(rsid,"zp"))
     insertImage(wb, paste0(rsid,"zp"), paste0("HRC/",rsid,"-000001.png"),width=18,height=12)
   }
   snp <- read.table(paste0("HRC/chr",f,"r.snp"),as.is=TRUE,header=TRUE)
   addWorksheet(wb, paste0(f,".snp"))
   writeDataTable(wb, paste0(f,".snp"), snp)
   config <- read.table(paste0("HRC/chr",f,"r.config"),as.is=TRUE,header=TRUE)
   addWorksheet(wb, paste0(f,".config"))
   writeDataTable(wb, paste0(f,".config"), config)
   saveWorkbook(wb, file=xlsx, overwrite=TRUE)
}
