# 3-11-2017 MRC-Epid JHZ

require(plink2R)
require(R2BGLiMS)
require(methods)
options(scipen=20, width=2000)
fp <- Sys.getenv("fp")
cat(fp,"\n")
# summary statistics
sumstats.name <- c("RS_ID","A1","A2","freqA1","b","se","P","N","chr","pos","SNP_ID")
sumstats <- read.table(paste0(fp,".dat"), as.is=TRUE, col.names=sumstats.name)
beta <- with(sumstats, b)
rsid <- with(sumstats, RS_ID)
snpid <- with(sumstats, SNP_ID)

# reference panel
p <- read_plink(fp)
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
sink(paste0(fp, ".jam"))
slot(j, "posterior.summary.table")
cbind(ssnpid, snpid[cc], rsid[cc])
sink()
sink(paste0(fp, ".top"))
TopModels(j)
sink()
