# 5-11-2017 MRC-Epid JHZ

options(width=200)

f <- Sys.getenv("f")
cat(f, "\n")
z <- read.table(paste0(f, ".z"), col.names=c("snp","z"))
snp <- read.table(paste0(f, ".snp"), header=TRUE)
index <- with(subset(snp, snp_log10bf>0), index)
ld <- read.table(paste0(f, ".ld"))
snp[1:length(index), ]
z[index, ]
ld[index, index]
