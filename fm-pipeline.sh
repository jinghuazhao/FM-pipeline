#!/bin/bash
# 1-11-2017 MRC-Epid JHZ

# GWAS summary statistics (the .sumstats file)
# filename containing list of lead SNPs
export snplist=2.snps
# GEN file
# in the form of chr{chr}_{start}_{end}.gen
# sample file
export sample_file=/gen_omics/data/EPIC-Norfolk/HRC/EPIC-Norfolk.sample
# sample for exclusion
export sample_to_exclude=/genetics/data/gwas/6-7-17/doc/exclude.dat
# -/+ flanking position
export flanking=25000
# number of threads
export threads=5
# software to be included in the analysis; change flags to 1 when available
# the outputs should be available individually from them
export locuszoom=0
export GCTA=0
export fgwas=0
export CAVIAR=0
export CAVIARBF=0
export finemap=1
export JAM=1
export FM_location=/genetics/bin/FM-pipeline

if [ $# -lt 1 ] || [ "$1" == "-h" ]; then
    echo "Usage: fm-pipeline.sh <input>"
    echo "where <input> is in sumstats format:"
    echo "SNP A1 A2 beta se N"
    echo "where SNP is RSid, A1 is effect allele"
    echo "and the outputs will be in <input>.out directory"
    exit
fi
if [ $(dirname $1) == "." ]; then
   dir=$(pwd)/$(basename $1).tmp
else
   dir=$1.tmp
fi
if [ ! -d $dir ]; then
   mkdir -p $dir
fi
cd $dir
cd /genetics/data/gwas/6-7-17/MAGIC
if $(test -f ${FM_location}/snp150.txt ); then
   echo "Chromosomal positions are ready to use"
   ln -sf ${FM_location}/snp150.txt
else
   echo "Obtaining chromosomal positions"
   wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp150Common.txt.gz
   gunzip -c snp150Common.txt.gz | cut -f2,4,5 | sort -k3,3 > snp150.txt
fi
echo supplement .sumstats with chromosomal positions
awk '{
  $2=toupper($2)
  $3=toupper($3)
}' $1 | join -11 -23 - snp150.txt | sed 's/chr//g' > $dir/$(basename $1).input
sort -k1,1 ${snplist} | join $dir/$(basename $1).input - > $dir/$(basename $1).lst
grep -w -f ${snplist} $dir/$(basename $1).input | awk -vs=$f{lanking} '{print $8,$9-s,$9+s}' > st.bed

cat $dir/$(basename $1).lst | parallel -j${threads} -C' ' \
'awk "(\$8==chr && \$9 >= pos-s && \$9 <= pos+s){\$2=toupper(\$2);\$3=toupper(\$3); \
 if(\$2<\$3) {a1=\$2; a2=\$3;} else {a1=\$3; a2=\$2}; \
 \$0=\$0 \" \" \$8 \":\" \$9 \"_\" a1 \"_\" a2;print}" chr={8} pos={9} s=${flanking} $dir/$(basename $1).input | \
 sort -k10,10 > chr{9}_{$({10}-${flanking})}_{$({10}+{flanking})}.dat'

echo "--> map/ped"
ls chr*.gen|sed 's/\.gen//g'|parallel -j${threads} --env wd -C' ' 'awk -f ${FM_location}/files/order.awk {}.gen > {}.ord;\
          gtool -G --g {}.ord --s ${sample_file} \
         --ped {}.ped --map {}.map --missing 0.05 --threshold 0.9 --log {}.log --snp --alleles \
         --chr $(echo {}|cut -d"_" -f1|sed "s/chr//g")'# --> auxiliary files
ls *.info|sed 's/\.info//g'|parallel -j${threads} -C' ' 'sort -k2,2 {}.map|join -110 -22 {}.dat -|sort -k10,10>{}.incl'
cat $wd/st.bed | parallel -j${threads} --env wd -C' ' 'f=chr{1}_{2}_{3};\
     awk "{print \$9,\$10,\$5,\$6,\$7,\$8,15234,\$11,\$1,\$6/\$7}" $f.incl > $f.r2;\
     cut -d" " -f9,10 $f.r2>$f.z;\
     awk "{print \$1}" $f.incl > $f.inc;\
     awk "{print \$1,\$4,\$3,\$14,\$15}" $f.incl > $f.a;\
     echo "RSID position chromosome A_allele B_allele" > $f.incl_variants;\
     awk "{print \$1,\$10,\$9,\$4,\$3}" $f.incl >> $f.incl_variants'
echo "--> bfile"
rm *bed *bim *fam
ls chr*.info|awk '(gsub(/\.info/,""))'|parallel -j${threads} --env wd -C' ' '\
         plink-1.9 --file {} --missing-genotype N --extract {}.inc --remove ${sample_to_exclude} \
         --make-bed --keep-allele-order --a2-allele {}.a 3 1 --out {}'
echo "--> bcor"
if [ $finemap -eq 1 ]; then
   ls *.info|sed 's/\.info//g'|parallel -j${threads} -C' ' '\
        ldstore --bcor {}.bcor --bplink {} --n-threads ${threads}; \  
        ldstore --bcor {}.bcor --merge ${threads}; \
        ldstore --bcor {}.bcor --matrix {}.ld --incl_variants {}.incl_variants; \
        sed -i -e "s/  */ /g; s/^ *//; /^$/d" {}.ld'
fi
echo "JAM, IPD"
ls chr*.info|awk '(gsub(/\.info/,""))'|parallel -j${threads} -C' ' '\
         plink-1.9 --bfile {} --indep-pairwise 500kb 5 0.80 --maf 0.05 --out {}; \
         grep -w -f {}.prune.in {}.a > {}.p; \
         grep -w -f {}.prune.in {}.dat > {}p.dat; \
         plink-1.9 --bfile {} --extract {}.prune.in --keep-allele-order --a2-allele {}.p 3 1 --make-bed --out {}p'
if [ $finemap -eq 1 ]; then
ls *.info|sed 's/\.info//g'|parallel -j${threads} -C' ' '\ 
         grep -w -f {}.prune.in {}.z > {}p.z; \
         ldstore --bcor {}p.bcor --bplink {}p --n-threads ${threads}; \
         ldstore --bcor {}p.bcor --merge ${threads}; \
         ldstore --bcor {}p.bcor --matrix {}p.ld; \
         sed -i -e "s/  */ /g; s/^ *//; /^$/d" {}p.ld'
echo "--> finemap"
   echo "z;ld;snp;config;log;n-ind" > finemap.cfg
   cat $wd/st.bed | parallel -j${threads} -C ' ' 'f=chr{1}_{2}_{3};sort -k7,7n $f.r2|tail -n1|cut -d" " -f7|\
   awk -vf=$f "{print sprintf(\"%s.z;%s.ld;%s.snp;%s.config;%s.log;%d\",f,f,f,f,f,int(\$1))}" >> finemap.cfg'
   finemap --sss --in-files finemap.cfg --n-causal-max 1 --corr-config 0.9
   sed 's/\./p\./g' finemap.cfg > finemapp.cfg
   finemap --sss --in-files finemapp.cfg --n-causal-max 1 --corr-config 0.9

   finemap --sss --in-files finemap.cfg --n-causal-max 3 --corr-config 0.9  
   finemap --sss --in-files finemapp.cfg --n-causal-max 3 --corr-config 0.9
fi

echo "--> JAM"
if [ $JAM -eq 1 ]; then
   cat $wd/st.bed|parallel -j${threads} --env wd -C' ' 'export fp=chr{1}_{2}_{3}p; R CMD BATCH --no-save ${FM_location}/files/JAM.R ${fp}.log'
fi
