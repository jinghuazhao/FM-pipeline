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
