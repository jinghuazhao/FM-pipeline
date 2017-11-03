#!/genetics/data/software/bin/Rscript --vanilla
# 3-11-2017 MRC-Epid JHZ

args <- commandArgs(trailingOnly = TRUE)
options(scipen=20,width=150)
getCredibleSNP <- function(snp, logProb, threshold=0.99){

 prob <- exp(logProb) 
 prob_normed <- prob/sum(prob)

 prob_cumsum <- cumsum(sort(prob_normed, decreasing=TRUE))
 nSNP <- as.numeric(which.max(prob_cumsum>threshold ))

 if(nSNP<=0){ nSNP=1 }

 credible_set <- snp[order(prob_normed, decreasing=TRUE)[1:nSNP]]
 select <- rep(FALSE, length(snp))
 select[order(prob_normed, decreasing=TRUE)[1:nSNP]] <- T 
 ret <- list(nSNP = nSNP, 
	prob_normed = prob_normed, 
	prob_cumsum = prob_cumsum[rank(-prob_normed, ties.method="random")], 
	credible_set=credible_set, 
	select=select 
 )
}
r <- args[1]
if(r==1)
{
  bed <- read.table(Sys.getenv("st.bed"), as.is=TRUE, col.names=c("chr","start","end","pos","rsid"),sep="\t")
  save(bed,file="bed.rda")
} else load("bed.rda")
f <- paste0("chr",bed[r,1],"_",bed[r,2],"_",bed[r,3])
dat <- read.table(paste0(f,".r"), as.is=TRUE,
       col.names=c("chr","pos","A","B","Freq1","Effect","StdErr","P","TOTALSAMPLESIZE","SNP",'SNPID'))
select <- with(dat,P)<0.1
max_SCZ_P <- min(dat$P[select])
max_SCZ_snp <- as.character(dat$SNP[select][which.min(dat$P[select])])
ret <- getCredibleSNP(as.character(dat$SNP[select]), qchisq(dat$P[select], 1, low=F)/2)

inCredible <- rep(NA, length(dat$P))
inCredible[match( dat$SNP[select], dat$SNP )] <- 0
inCredible[match( dat$SNP[select][ret$select], dat$SNP )] <- 1
prob_norm <- rep(NA, length(dat$P))
prob_norm[match(dat$SNP[select], dat$SNP )] <- ret$prob_normed
prob_cumsum <- rep(NA, length(dat$P))
prob_cumsum[match(dat$SNP[select], dat$SNP )] <- ret$prob_cumsum

result <- cbind(dat, inCredible=inCredible, probNorm=prob_norm, cumSum=prob_cumsum)

write.table(subset(result,!is.na(inCredible)&!is.na(prob_norm)&!is.na(prob_cumsum)), 
            file=paste0(f,".cre"), sep="\t", quote=F, col.names=T, row.names=F)
