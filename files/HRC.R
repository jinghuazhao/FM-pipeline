#!/genetics/data/software/bin/Rscript --vanilla
# 6-10-2017 MRC-Epid JHZ

options(scipen=20,width=150)

require(gap)
fm <- read.table("finemap.dat",as.is=TRUE,header=TRUE)
pdf("PPA.pdf")
par(xpd=TRUE,cex=0.6,srt=180)
ops <- mht.control(colors=rep(c("lightgray","lightblue"),11),srt=0,yline=2.5,xline=2,logscale=FALSE)
mhtplot2(data.frame(fm[,c("chr","pos","log10BF")],gene=NA,color=NA),ops,xlab="",ylab="",srt=0)
axis(2,at=c(3:13))
title("log10(BF) for strongly associated SNPs")
dev.off()

fm <- read.table("finemapp.dat",as.is=TRUE,header=TRUE)
pdf("PPAp.pdf")
par(xpd=TRUE,cex=0.6,srt=180)
ops <- mht.control(colors=rep(c("lightgray","lightblue"),11),srt=0,yline=2.5,xline=2,logscale=FALSE)
mhtplot2(data.frame(fm[,c("chr","pos","log10BF")],gene=NA,color=NA),ops,xlab="",ylab="",srt=0)
axis(2,at=c(3:13))
title("log10(BF) for strongly associated SNPs")
dev.off()

tophits <- read.table("/genetics/data/gwas/6-7-17/doc/region.dat",as.is=TRUE,col.names=c("rsid","chr","POS"),sep="\t")
bed <- read.table("/genetics/data/gwas/6-7-17/doc/st.bed",col.names=c("chr","start","end"))
r <- with(bed,table(chr))
R <- cbind(chr=1:22,n=0)
R[as.numeric(names(r)),] <- as.vector(r)
pdf("finemap_snp.pdf")
k <- 0
for(i in 1:22)
{
   if (R[i,2]==0) next
   par(mfrow=c(3,3))
   for(j in 1:R[i,2])
   {
   k <- k+1
   f <- paste(bed[k,1],bed[k,2],bed[k,3],sep="_")
   r2 <- read.table(paste0("chr",f,".r2"),as.is=TRUE,
         col.names=c("chr","pos","Freq1","Effect","StdErr","P-value","TOTALSAMPLESIZE","RS_ID","SNP_ID","z"))
   snp <- read.table(paste0("chr",f,".snp"),as.is=TRUE,header=TRUE)
   r2_snp <- merge(r2,snp,by.x="RS_ID",by.y="snp")
   r2_snp_tophits <- merge(r2_snp,tophits,by.x="RS_ID",by.y="rsid")
   head(r2);head(snp);head(tophits);head(r2_snp);head(r2_snp_tophits)
   ord <- with(r2_snp,order(index))
   snp <- r2_snp[ord,]
   with(r2_snp,plot(pos/100000,snp_log10bf,cex=0.3, xlab=paste0("chr",bed[k,1]," position (x100,000)"),ylab="log10(BF)",main=tophits[k,1]))
   with(r2_snp_tophits,points(pos/100000,snp_log10bf,col="red",pch=19))
   }
}
dev.off()

