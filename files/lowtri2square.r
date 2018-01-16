# 16-1-2018 MRC-Epid JHZ

f <- Sys.getenv("f")

m <- read.table(paste0(f,".magic.map"),skip=11,header=TRUE,as.is=TRUE)
nr <- nrow(m)
z <- gzfile(paste0(f,".magic.ld.gz"))
ld <- read.table(z, fill=TRUE, col.names=paste0("v",1:nr))
ld[upper.tri(ld)] <- t(ld)[upper.tri(ld)]
diag(ld) <- 1
write.table(ld,file=paste0(f,".ld"),col.names=FALSE,row.names=FALSE,quote=FALSE)
