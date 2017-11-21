# 21-11-2017 MRC-Epid JHZ

options(echo=FALSE, width=200)

f <- Sys.getenv("f")
cat(f, "\n")
config <- read.table(paste0(f,".config"),as.is=TRUE,header=TRUE)
subset(config,config_prob>0.01)
z <- read.table(paste0(f, ".z"), as.is=TRUE, col.names=c("snp","z"))
snp <- read.table(paste0(f, ".snp"), as.is=TRUE, header=TRUE)
index <- with(subset(snp, snp_log10bf>0), index)
ld <- read.table(paste0(f, ".ld"))
options(echo=TRUE)
snp[1:length(index), ]
z[index, ]
ld[index, index]
