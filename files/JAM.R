# 9-2-2018 MRC-Epid JHZ

require(plink2R)
# require(snpStats)
require(R2BGLiMS)
require(methods)
require(openxlsx)

options(scipen=20, width=2000)

f <- Sys.getenv("f")
cat(f,"\n")
bed <- paste0(f,".bed")
bim <- paste0(f,".bim")
fam <- paste0(f,".fam")

# summary statistics
sumstats.name <- c("RS_ID","A1","A2","freqA1","b","se","P","N","chr","pos","SNP_ID")
sumstats <- read.table(paste0(f,".dat"), as.is=TRUE, col.names=sumstats.name)
beta <- with(sumstats, b)
rsid <- with(sumstats, RS_ID)
snpid <- with(sumstats, SNP_ID)

# reference panel with mean substitution for (small) proportion of missing data
p <- read_plink(f)
R <- with(p, as.data.frame(2-bed))
# p <- read.plink(bed,bim,fam)
# R <- as(with(p,genotypes),"numeric")
R[] <- lapply(R, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})
X.ref <- R

# JAM modeling
ssnpid <- paste0("snp", 1:length(beta))
names(beta) <- colnames(X.ref) <- ssnpid
priors <- list("a"=1, "b"=length(beta), "Variables"=ssnpid)
n <- 15234
j <- JAM(marginal.betas=beta, n=n, X.ref=X.ref, n.mil=5, tau=n, full.mcmc.sampling = FALSE, model.space.priors=priors)
save(j,file=paste0(f,".j"))
pst <- slot(j, "posterior.summary.table")
tm <- TopModels(j)
ssr <- data.frame(ssnpid=ssnpid, snpid=snpid, rsid=rsid)
cs <- CredibleSet(j, credible.percentile.threshold=0.75)
msbf <- ModelSizeBayesFactors(j)[[1]]
sink(paste0(f, ".jam"))
pst
ssr
cat("\nCredible set\n")
cs
cat("\nModel size Bayes Factors\n")
msbf
sink()
sink(paste0(f, ".top"))
tm
sink()
n.col <- ncol(tm)
n.snps <- n.col-1
post.prob <- tm[,n.col]
n.sel <- apply(tm[,1:n.snps],1,sum)
sink(paste0(f,".sum"))
cbind(n.sel,post.prob)
sink()
sink(paste0(f,".cs"))
cbind(subset(ssr,ssnpid%in%cs),subset(pst,rownames(pst)%in%cs))
sink()
if(identical(cs,character(0))) unlink(paste0(f,".cs"))
tm1 <- tm[1,-n.col]
selected <- names(tm1[tm1==1])
if(n.sel[1]>0&n.sel[1]!=n.snps)
{
   PostProb_model <- rep(post.prob[1],n.sel[1])
   t <- cbind(subset(ssr,ssnpid%in%selected), PostProb_model, subset(pst,rownames(pst)%in%selected))
   write.table(t,paste0(f,".sel"),row.names=FALSE,quote=FALSE)
}
png(paste0(f,".png"), units = 'in', width=18, height=12, res=300)
ManhattanPlot(j)
dev.off()

xlsx <- paste0(f,".xlsx")
wb <- createWorkbook(xlsx)
addWorksheet(wb, "ID")
writeDataTable(wb, "ID", ssr)
addWorksheet(wb, "TopModels")
writeDataTable(wb, "TopModels", as.data.frame(tm))
addWorksheet(wb, "Model.1")
PostProb_model <- rep(post.prob[1],n.sel[1])
writeDataTable(wb, "Model.1", cbind(subset(ssr,ssnpid%in%selected),PostProb_model,subset(pst,rownames(pst)%in%selected)))
addWorksheet(wb, "CredibleSet")
writeDataTable(wb, "CredibleSet", cbind(subset(ssr,ssnpid%in%cs),subset(pst,rownames(pst)%in%cs)))
addWorksheet(wb, "ModelSizeBayesFactors")
writeDataTable(wb, "ModelSizeBayesFactors", as.data.frame(msbf))
addWorksheet(wb, "posterior.summary.table")
writeDataTable(wb, "posterior.summary.table", cbind(ID=rownames(pst), as.data.frame(pst)))
addWorksheet(wb, "Manhattan.plot")
insertImage(wb, "Manhattan.plot", paste0(f, ".png"), width=18, height=12)
saveWorkbook(wb, file=xlsx, overwrite=TRUE)

# obsolete as it only deals with complete data
# cc <- complete.cases(t(R))
# beta <- beta[cc]
# X.ref <- R[,cc]
# ssnpid <- paste0("snp", 1:length(beta[cc]))
# ssr <- data.frame(ssnpid=ssnpid, snpid=snpid[cc], rsid=rsid[cc])
