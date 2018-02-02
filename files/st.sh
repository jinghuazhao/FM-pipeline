# 2-2-2018 MRC-Epid JHZ

# reformat summary statistics
gunzip -c meta_meno1kg_101115_filt_rsid_flags_KR.txt.gz | \
awk -vFS="\t" -vOFS="\t" '
(NR>1) {
  split($1,a,":")
  gsub(/chr/,"",a[1])
  if ($17==".") $17=$1
  print $17, $2, $3, $4, $8, $9, $10, $16, a[1], a[2]
}' | \
sort -k1,1 > repro.txt

export region=/genetics/data/gwas/6-7-17/doc/region.dat
cut -f1 $region > 97.snps

awk -vfl=250000 '{
  if(NR==1) print "chr start end pos rsid r"
  gsub(/,/,"",$3);l=$3-fl;u=$3+fl;if(l<0) l=1;print $2,l,u,$3,$1,NR
}' $region > st.bed

grep -v -f exclude.dat /gen_omics/data/EPIC-Norfolk/HRC/EPIC-Norfolk.sample > HRC.sample

cd /scratch/tempjhz22/LDcalc/HRC
seq 22|awk -vp=chr '{print p $1}' > merge-list
plink-1.9 --merge-list merge-list --maf 0.000072 --make-bed --out HRC
cd -
