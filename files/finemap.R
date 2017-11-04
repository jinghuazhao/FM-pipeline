#!/genetics/data/software/bin/Rscript --vanilla
# 4-11-2017 MRC-Epid JHZ

options(scipen=20, width=150)

require(gap)
fm <- read.table("finemap.dat", as.is=TRUE, header=TRUE)
pdf("PPA.pdf")
par(xpd=TRUE, cex=0.6, srt=180)
ops <- mht.control(colors=rep(c("lightgray","lightblue"),11), srt=0, yline=2.5, xline=2, logscale=FALSE)
mhtplot2(data.frame(fm[,c("chr","pos","log10BF")], gene=NA, color=NA), ops, xlab="", ylab="", srt=0)
axis(2,at=c(3:13))
title("log10(BF) for strongly associated SNPs")
dev.off()

fm <- read.table("finemapp.dat", as.is=TRUE, header=TRUE)
pdf("PPAp.pdf")
par(xpd=TRUE,cex=0.6,srt=180)
ops <- mht.control(colors=rep(c("lightgray","lightblue"),11), srt=0, yline=2.5, xline=2, logscale=FALSE)
mhtplot2(data.frame(fm[,c("chr","pos","log10BF")], gene=NA, color=NA), ops, xlab="", ylab="", srt=0)
axis(2,at=c(3:13))
title("log10(BF) for strongly associated SNPs")
dev.off()

bed <- read.table("st.bed", as.is=TRUE, header=TRUE)
r <- with(bed, table(chr))
R <- cbind(chr=1:22, n=0)
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
   f <- paste(bed[k,1:3], sep="_")
   dat <- read.table(paste0("chr",f,".dat"), as.is=TRUE,
          col.names=c("RS_ID","A1","A2","freqA1","b","se","P","N","chr","pos","SNP_ID"))
   snp <- read.table(paste0("chr",f,".snp"), as.is=TRUE, header=TRUE)
   dat_snp <- merge(dat, snp, by.x="RS_ID", by.y="snp")
   dat_snp_bed <- merge(dat_snp, bed, by.x="RS_ID", by.y="rsid")
   head(dat);head(snp);head(bed);head(dat_snp);head(dat_snp_bed)
   ord <- with(dat_snp,order(index))
   snp <- dat_snp[ord,]
   with(dat_snp,plot(pos/100000,snp_log10bf,cex=0.3, xlab=paste0("chr",bed[k,1]," position (x100,000)"),ylab="log10(BF)",main=tophits[k,1]))
   with(dat_snp_bed, points(pos/100000,snp_log10bf, col="red", pch=19))
   }
}
dev.off()
