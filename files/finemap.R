# 10-9-2019 JHZ

options(digits=3, scipen=20, width=500)

f <- Sys.getenv("f")
cat(f, "\n")
z <- read.table(paste0(f, ".fm.z"), as.is=TRUE, header=TRUE)
ld <- read.table(paste0(f, ".ld"))
snp <- read.table(paste0(f, ".snp"), as.is=TRUE, header=TRUE)
config <- read.table(paste0(f,".config"),as.is=TRUE,header=TRUE)
subset(snp, prob>0.01)
subset(config, prob>0.01)
z <- subset(snp, abs(z) > 5.73)
id <- as.integer(row.names(z))
ldt <- ld[id,id]
ldt[upper.tri(ldt, diag=TRUE)] <- NA
colnames(ldt) <- with(z, rsid)
chk <- cbind(z[c("rsid","z")], ldt)
chk

library(openxlsx)
xlsx <- paste0(f, "-finemap.xlsx")
unlink(xlsx, recursive = FALSE, force = TRUE)
wb <- createWorkbook(xlsx)
addWorksheet(wb, paste0(f, ".snp"))
writeDataTable(wb, paste0(f, ".snp"), snp)
addWorksheet(wb, paste0(f, ".cfg"))
writeDataTable(wb, paste0(f, ".cfg"), config)
addWorksheet(wb, paste0(f, ".chk"))
writeDataTable(wb, paste0(f, ".chk"), chk)
CredibleSet <- function()
{
## Code from Ji Chen for number of configs that account for % probability and variants in the credible set
  end <- which( cumsum( config$prob ) >= 0.75 )[ 1 ]
  end
  credible_set <- unique( strsplit( paste( config$config[ seq( end ) ], collapse = ',' ), split = ',' )[[ 1 ]] )
  head(config, end)
  snplist <- with(subset(snp,(snp%in%credible_set)&prob>0),snp)
  length(snplist)
  cs <- cbind(subset(z,snp%in%snplist),with(subset(snp, snp%in%snplist),cbind(prob,log10bf)))
  cs
  addWorksheet(wb, paste0(f, ".cs"))
  writeDataTable(wb, paste0(f, ".cs"), cs)
}
saveWorkbook(wb, file=xlsx, overwrite=TRUE)
