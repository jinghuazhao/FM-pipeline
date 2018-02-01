# 1-2-2018 MRC-Epid JHZ

options(scipen=20, width=200)

if(file.exists("gcta-slct.csv")&file.exists("jam.cs"))
{

  slct <- read.csv("gcta-slct.csv",as.is=TRUE)
  cs <- read.table("jam.cs",header=TRUE,as.is=TRUE)
  stbed <- read.table("st.bed",as.is=TRUE,header=TRUE)
  require(openxlsx)
  xlsx <- "gcta-jam-finemap.xlsx"
  wb <- createWorkbook(xlsx)
  addWorksheet(wb, "gcta")
  writeDataTable(wb, "gcta", slct)
  addWorksheet(wb, "jam")
  writeDataTable(wb, "jam", cs)
  for(i in seq(nrow(stbed)))
  {
    r <- paste0("chr",stbed[i,1],"_",stbed[i,2],"_",stbed[i,3])
    slct.r <- subset(slct,region=="r")
    cs.r <- subset(cs,region=="r")
    slct_cs.r <- cbind(slct.r[c("region","SNP","rsid","pJ")],cs.r[c("snpid","PostProb","BF")])
    slct_cs <- paste0(r,"slct_cs")
    addWorksheet(wb, slct_cs)
    writeDataTable(wb, slct_cs, slct_cs.r)
    snplist <- unique(sort(c(with(slct.r, SNP),with(cs.r, snpid))))
    fm <- paste0(r,c(".snp",".config",".ld",".z"))
    if(file.exists(fm[1])&file.exists(fm[2])&file.exists(fm[3])&file.exists(fm[4]))
    {
      snp <- read.table(fm[1],as.is=TRUE,header=TRUE)
      config <- read.table(fm[2],as.is=TRUE,header=TRUE)
      ld <- read.table(fm[3])
      z <- read.table(fm[4],col.names=c("snp","z"))
      index <- with(subset(snp,snp%in%snplist),index)
      sumstat <- subset(z,snp%in%snplist)
      ld[index,index][upper.tri(ld[index,index])] <- NA
      check <- data.frame(index=rownames(ld[index,index]),sumstat,ld[index,index])
      slct_cs_check <- paste0(r,".check")
      addWorksheet(wb, slct_cs_check)
      writeDataTable(wb, slct_cs_check, check)
    }
  }
  saveWorkbook(wb, file=xlsx, overwrite=TRUE)
}
