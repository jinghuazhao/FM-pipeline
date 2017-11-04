#!/bin/bash
# 6-10-2017 MRC-Epid JHZ

export wd=/genetics/data/gwas/6-7-17
cd HRC

# finemap, Bayes Facttor > 20
echo "snpid region index snp_prob snp_log10bf rsid" > finemap.K20
seq 97|parallel -j5 --env wd -C' ' 'read chr start end <<<$(awk -vline={} "NR==line" $wd/doc/st.bed);\
f=chr${chr}_${start}_${end};cut -d" " -f10,11 ${f}.r > ${f}.tmp; \
awk "(NR>1&&\$3>0.8&&\$4>1.3){gsub(/.snp/,\"\",FILENAME);print FILENAME,\$0}" ${f}.snp|sort -k3,3|join -13 -22 - ${f}.tmp >> finemap.K20'
echo "chr pos log10BF prob snpid rsid region" > finemap.dat
awk '(NR>1){snpid=$1;gsub(/:|_/," ",$1);split($1,a," ");print a[1],a[2],$5,$4,snpid,$6,$2}' finemap.K20|sort -k1,1n -k2,2n >> finemap.dat
# The pruned data
echo "snpid region index snp_prob snp_log10bf rsid" > finemapp.K20   
seq 97|parallel -j5 --env wd -C' ' 'read chr start end <<<$(awk -vline={} "NR==line" $wd/doc/st.bed);\
f=chr${chr}_${start}_${end};cut -d" " -f10,11 ${f}.r > ${f}.tmp; \
awk "(NR>1&&\$3>0.8&&\$4>1.3){gsub(/.snp/,\"\",FILENAME);print FILENAME,\$0}" ${f}p.snp|sort -k3,3|join -13 -22 - ${f}.tmp >> finemapp.K20'
echo "chr pos log10BF prob snpid rsid region" > finemapp.dat
awk '(NR>1){snpid=$1;gsub(/:|_/," ",$1);split($1,a," ");print a[1],a[2],$5,$4,snpid,$6,$2}' finemapp.K20|sort -k1,1n -k2,2n >> finemapp.dat
$wd/HRC.R

# JAM
cat $wd/doc/st.bed|parallel -j6 --env wd -C' ' 'echo chr{1}_{2}_{3}; $wd/JAM.R {1} {2} {3} > chr{1}_{2}_{3}.log'

# FM-summary
seq 97|parallel -j15 --env wd -C' ' '$wd/getCredible.r {}'
echo "region chr pos A B Freq1 Effect StdErr P TOTALSAMPLESIZE SNP inCredible probNorm cumSum"|sed 's/ /\t/g' > FM.txt
cat $wd/doc/st.bed|parallel -j5 -C' ' '\
awk "!/SNP/{gsub(/\.cre/,\"\",FILENAME);print FILENAME, \$0}" OFS="\t" chr{1}_{2}_{3}.cre >> FM.txt'

# lcuszoom plots
echo "{OFS=\"\\t\";if(NR==1) print \"MarkerName\",\"P-value\",\"Weight\"; print \$8,\$6,\$7}" > lz.awk
ls *info | sed 's/\.info//g' | parallel -j10 -C' ' 'awk -f lz.awk {}.r2> {}.lz'
for i in $(seq 97)
do
  export refsnp=$(awk -vi=$i '(NR==i){print $1}' $wd/doc/region.dat)
  awk -vi=$i '(NR==i)' $wd/doc/st.bed|\
  parallel -j1 --env refsnp -C' ' 'locuszoom-1.3 --metal chr{1}_{2}_{3}.lz \
           --refsnp $refsnp --plotonly --no-date;pdftopng $refsnp.pdf -r 300 $refsnp'
done
cd -
