# 3-2-2018 MRC-Epid JHZ

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

# Those at end of gcta-slct.sh are actually fine so thes following section using qctool_v2.0 is unnecessary.

export s=/gen_omics/data/EPIC-Norfolk/HRC
export o=/scratch/tempjhz22/LDcalc/HRC
export w=/genetics/data/gwas/1-11-17

seq 22 | parallel -j1 --env s --env o --env w -C' ' '
    sge "/genetics/bin/qctool_v2.0 -filetype bgen -g $s/chr{}.bgen -s $s/EPIC-Norfolk.sample -excl-samples $w/exclude.dat -snp-stats -osnp $o/chr{}.snpstats"'

seq 22 | parallel -j5 --env o -C' ' 'awk "NR>12 && (\$14<0.000072||\$18<0.4||\$19>=0.05){print \$2}" $o/chr{}.snpstats > $o/chr{}.excl'
seq 22 | parallel -j1 --env i --env b --env o -C' ' '
    qctool_v2.0 -filetype bgen -g $s/chr{}.bgen -s $s/EPIC-Norfolk.sample \
               -ofiletype binary_ped -og $o/chr{} -excl-samples exclude.dat \
               -excl-rsids $o/chr{}.excl -threshhold 0.9 -threads 5'

cd $o
cd /gen_omics/data/EPIC-Norfolk/HRC/binary_ped
ln -f $o/HRC.bed
ln -f $o/HRC.bim
ln -f $o/HRC.fam
ln -f $o/HRC.nosex
seq 22|awk -vp=chr '{print p $1}' > merge-list
plink-1.9 --merge-list merge-list --maf 0.000072 --make-bed --out HRC
cd -

qsub -S /bin/bash -V -N HRC -cwd -e HRC3.err -o HRC3.out -pe make 10 -q all.q /genetics/bin/gcta-slct.sh HRC
