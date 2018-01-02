# 2-1-2018 MRC-Epid JHZ

options(echo=FALSE, width=200)

f <- Sys.getenv("f")
cat(f, "\n")
z <- read.table(paste0(f, ".z"), as.is=TRUE, col.names=c("snp","z"))
ld <- read.table(paste0(f, ".ld"))
snp <- read.table(paste0(f, ".snp"), as.is=TRUE, header=TRUE)
index <- with(subset(snp, snp_prob>0.01), index)
config <- read.table(paste0(f,".config"),as.is=TRUE,header=TRUE)
options(echo=TRUE)
subset(snp, snp_prob>0.01)
subset(config,config_prob>0.01)
chk <- cbind(z[index, ], snp[1:length(index), -2], ld[index, index])
chk
save(chk, file=paste0(f,".chk"))
