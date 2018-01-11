# 11-1-2018 MRC-Epid JHZ

require(plink2R)
require(R2BGLiMS)
require(methods)
options(scipen=20, width=2000)
f <- Sys.getenv("f")
cat(f,"\n")
# summary statistics
sumstats.name <- c("RS_ID","A1","A2","freqA1","b","se","P","N","chr","pos","SNP_ID")
sumstats <- read.table(paste0(f,".dat"), as.is=TRUE, col.names=sumstats.name)
beta <- with(sumstats, b)
rsid <- with(sumstats, RS_ID)
snpid <- with(sumstats, SNP_ID)

# reference panel
p <- read_plink(f)
R <- with(p, as.matrix(bed))
cc <- complete.cases(t(R))
beta <- beta[cc]
X.ref <- R[,cc]

# JAM modeling
ssnpid <- paste0("snp", 1:length(beta[cc]))
names(beta) <- colnames(X.ref) <- ssnpid
priors <- list("a"=1, "b"=length(beta), "Variables"=ssnpid)
n <- 15234
j <- JAM(marginal.betas=beta, n=n, X.ref=X.ref, n.mil=5, tau=n, full.mcmc.sampling = TRUE, model.space.priors=priors)
pst <- slot(j, "posterior.summary.table")
tm <- TopModels(j)
ssr <- data.frame(ssnpid=ssnpid, snpid=snpid[cc], rsid=rsid[cc])
sink(paste0(f, ".jam"))
pst
ssr
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
tm1 <- tm[1,-n.col]
selected <- names(tm1[tm1==1])
if(n.sel[1]>0&n.sel[1]!=n.snps)
{
   pp_model <- rep(post.prob[1],n.sel[1])
   t <- cbind(subset(ssr,ssnpid%in%selected), pp_model, subset(pst,rownames(pst)%in%selected))
   write.table(t,paste0(f,".sel"),row.names=FALSE,quote=FALSE)
}
