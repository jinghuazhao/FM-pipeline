# 16-1-2018 MRC-Epid JHZ

f <- Sys.getenv("f")
z <- gzfile(paste0(f,".magic.ld.gz"))
m <- read.table(z, fill=TRUE, col.names=paste0("v",))
m[upper.tri(m)] <- t(m)[upper.tri(m)]
write.table(m,file=paste0(f,".ld"),col.names=FALSE,row.names=FALSE,quote=FALSE)

