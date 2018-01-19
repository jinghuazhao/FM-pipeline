# 15-1-2018 MRC-Epid JHZ

options(scipen=20, width=200)

# collate all Credible sets
cs <- data.frame()
stbed <- read.table("st.bed",as.is=TRUE,header=TRUE)
for(i in seq(nrow(stbed)))
{
  f <- paste0("chr",stbed[i,1],"_",stbed[i,2],"_",stbed[i,3],"p.cs")
  if(file.exists(f)) cs <- rbind(cs,read.table(f,as.is=TRUE,header=TRUE))
}

# split snpid
cs <- within(cs, {
   sp <- lapply(lapply(snpid, function(x) strsplit(x,":|_")),unlist)
   Chr <- Pos <- NA
})

# obtain positions
cs <- within(cs, {
   for(i in seq(nrow(cs))) {
      Chr[i] <- as.integer(sp[[i]][1])
      Pos[i] <- as.integer(sp[[i]][2])
   }
})

# order by positions
ord <- with(cs, order(Chr,Pos))
cs <- within(cs[ord,],{i <- ssnpid <- sp <- NULL})
write.table(cs, file="jam.cs", row.names=FALSE, quote=FALSE)
if(file.exists("gcta-slct.csv"))
{
  require(openxlsx)
  slct <- read.csv("gcta-slct.csv",as.is=TRUE)
  m1 <- read.table("jam.out",as.is=TRUE,header=TRUE)
  cs <- read.table("jam.cs",header=TRUE,as.is=TRUE)
  gcta_m1 <- merge(m1[c("SNP","PostProb_model","PostProb","BF")],slct,by="SNP")
  gcta_cs <- merge(cs[c("snpid","PostProb","BF","Pos")],slct,by.x="snpid",by.y="SNP")
  ord <- with(gcta_cs,order(Chr,Pos))
  xlsx <- "gcta-jam.xlsx"
  wb <- createWorkbook(xlsx)
  addWorksheet(wb, "gcta")
  writeDataTable(wb, "gcta", slct)
  addWorksheet(wb, "m1")
  writeDataTable(wb, "m1", m1)
  addWorksheet(wb, "gcta-m1")
  writeDataTable(wb, "gcta-m1", gcta_m1)
  addWorksheet(wb, "cs")
  writeDataTable(wb, "cs", cs)
  addWorksheet(wb, "gcta-cs")
  writeDataTable(wb, "gcta-cs", gcta_cs[ord,])
  saveWorkbook(wb, file=xlsx, overwrite=TRUE)
}
