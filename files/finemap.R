# 14-1-2018 MRC-Epid JHZ

options(digits=3, scipen=20, width=500)

f <- Sys.getenv("f")
cat(f, "\n")
z <- read.table(paste0(f, ".z"), as.is=TRUE, col.names=c("snp","z"))
ld <- read.table(paste0(f, ".ld"))
snp <- read.table(paste0(f, ".snp"), as.is=TRUE, header=TRUE)
config <- read.table(paste0(f,".config"),as.is=TRUE,header=TRUE)
subset(snp, snp_prob>0.01)
subset(config,config_prob>0.01)
id <- with(subset(snp, snp_prob>0.01), index)
chk <- cbind(z[id, ], with(subset(snp,index%in%id), c(snp_prob,snp_log10bf)), ld[id, id])
chk

## Code from Ji Chen for number of configs that account for % probability and variants in the credible set
end <- which( cumsum( config$config_prob ) >= 0.75 )[ 1 ]
end
credible_set <- unique( strsplit( paste( config$config[ seq( end ) ], collapse = ',' ), split = ',' )[[ 1 ]] )
head(config, end)
snplist <- with(subset(snp,(snp%in%credible_set)&snp_prob>0),snp)
length(snplist)
cs <- cbind(subset(z,snp%in%snplist),with(subset(snp, snp%in%snplist),cbind(snp_prob,snp_log10bf)))
cs

library(openxlsx)
xlsx <- paste0(f, ".xlsx")
unlink(xlsx, recursive = FALSE, force = TRUE)
wb <- createWorkbook(xlsx)
addWorksheet(wb, paste0(f, ".snp"))
writeDataTable(wb, paste0(f, ".snp"), snp)
addWorksheet(wb, paste0(f, ".cfg"))
writeDataTable(wb, paste0(f, ".cfg"), config)
addWorksheet(wb, paste0(f, ".chk"))
writeDataTable(wb, paste0(f, ".chk"), chk)
addWorksheet(wb, paste0(f, ".cs"))
writeDataTable(wb, paste0(f, ".cs"), cs)
saveWorkbook(wb, file=xlsx, overwrite=TRUE)
