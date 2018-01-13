# 13-1-2018 MRC-Epid JHZ

options(digits=3, scipen=20, width=200)

f <- Sys.getenv("f")
cat(f, "\n")
z <- read.table(paste0(f, ".z"), as.is=TRUE, col.names=c("snp","z"))
ld <- read.table(paste0(f, ".ld"))
snp <- read.table(paste0(f, ".snp"), as.is=TRUE, header=TRUE)
index <- with(subset(snp, snp_prob>0.01), index)
config <- read.table(paste0(f,".config"),as.is=TRUE,header=TRUE)
subset(snp, snp_prob>0.01)
subset(config,config_prob>0.01)
chk <- cbind(z[index, ], snp[1:length(index), c(3,4)], ld[index, index])
chk
save(chk, file=paste0(f,".chk"))

## Adapted from code by Ji Chen 
## Find number of configs that account for 99% probability mass
end <- which( cumsum( config$config_prob ) >= 0.99 )[ 1 ]
head(config, end)

## Get variants in 99% credible set
credible_set <- unique( strsplit( paste( config$config[ seq( end ) ], collapse = ',' ), split = ',' )[[ 1 ]] )

index <- with(subset(snp, snp%in%credible_set), index)
cs <- cbind(z[index, ], snp[index, c(3,4)], ld[index, index])
cs

library(openxlsx)
xlsx <- paste0("chr", f, "-", rsid, ".xlsx")
unlink(xlsx, recursive = FALSE, force = TRUE)
wb <- createWorkbook(xlsx)
snp <- read.table(paste0("chr", f, ".snp"), as.is=TRUE, header=TRUE)
addWorksheet(wb, paste0(f, ".snp"))
writeDataTable(wb, paste0(f, ".snp"), snp)
config <- read.table(paste0("chr", f, ".config"), as.is=TRUE, header=TRUE)
addWorksheet(wb, paste0(f, ".cfg"))
writeDataTable(wb, paste0(f, ".cfg"), config)
load(paste0("chr", f, ".chk"))
addWorksheet(wb, paste0(f, ".chk"))
writeDataTable(wb, paste0(f, ".chk"), chk)
addWorksheet(wb, paste0(f, ".cs"))
writeDataTable(wb, paste0(f, ".cs"), cs)
saveWorkbook(wb, file=xlsx, overwrite=TRUE)
