# 5-1-2018 MRC-Epid JHZ

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
ssr <- cbind(ssnpid, snpid[cc], rsid[cc])
sink(paste0(f, ".jam"))
pst
ssr
sink()
tm <- TopModels(j)
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
if(length(selected)>0&length(selected)!=n.col-1) cbind(selected,post.prob[1],ssr[ssr[,2]%in%selected,],pst[rownames(pst)%in%selected,])
