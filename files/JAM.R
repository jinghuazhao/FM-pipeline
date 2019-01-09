# 9-1-2019 MRC-Epid JHZ

require(plink2R)
# require(snpStats)
require(R2BGLiMS)
require(methods)
require(openxlsx)

options(scipen=20, width=2000)

fp <- Sys.getenv("f")
f <- paste0(fp,"p")
cat(fp,"\n")
bed <- paste0(f,".bed")
bim <- paste0(f,".bim")
fam <- paste0(f,".fam")

# summary statistics
sumstats.name <- c("RS_ID","A1","A2","freqA1","b","se","P","N","chr","pos","SNP_ID")
sumstats <- read.table(paste0(fp,".dat"), as.is=TRUE, col.names=sumstats.name)

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

ss <- subset(sumstats,SNP_ID%in%p$bim$V2)
beta <- with(ss, b)
rsid <- with(ss, RS_ID)
snpid <- with(ss, SNP_ID)

# JAM modeling
ssnpid <- paste0("snp", 1:length(beta))
names(beta) <- colnames(X.ref) <- ssnpid
priors <- list("a"=1, "b"=length(beta), "Variables"=ssnpid)
n <- 15234
j <- JAM(marginal.betas=beta, n=n, X.ref=X.ref, n.mil=5, tau=n, full.mcmc.sampling = FALSE, model.space.priors=priors)
save(j,file=paste0(fp,".j"))
pst <- slot(j, "posterior.summary.table")
tm <- TopModels(j)
ssr <- data.frame(ssnpid=ssnpid, snpid=snpid, rsid=rsid)
cs <- CredibleSet(j, credible.percentile.threshold=0.75)
msbf <- ModelSizeBayesFactors(j)[[1]]
sink(paste0(fp, ".jam"))
pst
ssr
cat("\nCredible set\n")
cs
cat("\nModel size Bayes Factors\n")
msbf
sink()
sink(paste0(fp, ".top"))
tm
sink()
n.col <- ncol(tm)
n.snps <- n.col-1
post.prob <- tm[,n.col]
n.sel <- apply(tm[,1:n.snps],1,sum)
sink(paste0(fp,".sum"))
cbind(n.sel,post.prob)
sink()
sink(paste0(fp,".cs"))
cbind(subset(ssr,ssnpid%in%cs),subset(pst,rownames(pst)%in%cs))
sink()
if(identical(cs,character(0))) unlink(paste0(fp,".cs"))
tm1 <- tm[1,-n.col]
selected <- names(tm1[tm1==1])
if(n.sel[1]>0&n.sel[1]!=n.snps)
{
   PostProb_model <- rep(post.prob[1],n.sel[1])
   t <- cbind(subset(ssr,ssnpid%in%selected), PostProb_model, subset(pst,rownames(pst)%in%selected))
   write.table(t,paste0(fp,".sel"),row.names=FALSE,quote=FALSE)
}

xlsx <- paste0(fp,".xlsx")
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
saveWorkbook(wb, file=xlsx, overwrite=TRUE)

# obsolete as it only deals with complete data
# cc <- complete.cases(t(R))
# beta <- beta[cc]
# X.ref <- R[,cc]
# ssnpid <- paste0("snp", 1:length(beta[cc]))
# ssr <- data.frame(ssnpid=ssnpid, snpid=snpid[cc], rsid=rsid[cc])
