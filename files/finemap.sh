#!/bin/bash
# 4-11-2017 MRC-Epid JHZ

echo "snpid region index snp_prob snp_log10bf rsid" > finemap.K20
awk 'NR>1' st.bed | \
parallel -j${threads} -C' ' '\
    export f=chr{1}_{2}_{3}; \
    cut -d" " -f10,11 $f.r > $f.tmp; \
    awk "(NR>1&&\$3>0.8&&\$4>1.3){print ENVIRON["f"], \$0}" $f.snp | \
    sort -k3,3 | \
    join -13 -22 - $f.tmp >> finemap.K20'
echo "chr pos log10BF prob snpid rsid region" > finemap.dat
awk '(NR>1){snpid=$1;gsub(/:|_/," ",$1);split($1,a," ");print a[1],a[2],$5,$4,snpid,$6,$2}' finemap.K20 | \
sort -k1,1n -k2,2n >> finemap.dat
# The pruned data
echo "snpid region index snp_prob snp_log10bf rsid" > finemapp.K20
awk 'NR>1' | \
parallel -j${threads} -C' ' '\
    export f=chr{1}_{2}_{3}; \
    cut -d" " -f10,11 $f.r > $f.tmp; \
    awk "(NR>1&&\$3>0.8&&\$4>1.3){print ENVIRON["f"], \$0}" ${f}p.snp | \
    sort -k3,3 | \
    join -13 -22 - $f.tmp >> finemapp.K20'
echo "chr pos log10BF prob snpid rsid region" > finemapp.dat
awk '(NR>1){snpid=$1;gsub(/:|_/," ",$1);split($1,a," ");print a[1],a[2],$5,$4,snpid,$6,$2}' finemapp.K20 | \
sort -k1,1n -k2,2n >> finemapp.dat
