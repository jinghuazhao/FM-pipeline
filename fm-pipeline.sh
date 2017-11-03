#!/bin/bash
# 3-11-2017 MRC-Epid JHZ

## settings -- change as apporopriate
# working directory
export wd=/genetics/data/gwas/6-7-17/MAGIC
# GWAS summary statistics (the .sumstats file)
export args=$wd/2hrglucose.txt
# filename containing list of lead SNPs
export snplist=$wd/2.snps
# GEN files, named chr{chr}_{start}_{end}.gen
export GEN_location=$wd/MAGIC
# sample file
export sample_file=/gen_omics/data/EPIC-Norfolk/HRC/EPIC-Norfolk.sample
# sample for exclusion
export sample_to_exclude=$wd/exclude.dat
# -/+ flanking position
export flanking=250000
# only generate st.bed containg chr, start, end triplets
export stbed=0
# N, study sample size as used by finemap
export N=15234
# number of threads
export threads=5
# software to be included in the analysis; change flags to 1 when available
# the outputs should be available individually
export CAVIAR=0
export CAVIARBF=0
export finemap=1
export JAM=1
export locuszoom=0
export fgwas=0
export GCTA=0
export FM_location=/genetics/bin/FM-pipeline

if [ $# -lt 1 ] || [ "$args" == "-h" ]; then
    echo "Usage: fm-pipeline.sh <input>"
    echo "where <input> is in sumstats format:"
    echo "SNP A1 A2 beta se N"
    echo "where SNP is RSid, A1 is effect allele"
    echo "and the outputs will be in <input>.out directory"
    exit
fi
if [ $(dirname $args) == "." ]; then
   dir=$(pwd)/$(basename $args).tmp
else
   dir=$args.tmp
fi
if [ ! -d $dir ]; then
   mkdir -p $dir
fi
cd $dir
if $(test -f ${FM_location}/snp150.txt ); then
   echo "Chromosomal positions are ready to use"
   ln -sf ${FM_location}/snp150.txt
else
   echo "Obtaining chromosomal positions"
   wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp150Common.txt.gz
   gunzip -c snp150Common.txt.gz | \
   cut -f2,4,5 | \
   sort -k3,3 > snp150.txt
fi
echo supplement .sumstats with chromosomal positions
export rt=$dir/$(basename $args)
awk '{
  $2=toupper($2)
  $3=toupper($3)
};1' $args | \
join -11 -23 - snp150.txt | \
sed 's/chr//g' > $rt.input
sort -k1,1 ${snplist} | \
join $dir/$(basename $args).input - > $rt.lst
grep -w -f ${snplist} $rt.input | \
awk -vs=${flanking} '{print $7,$8-s,$8+s, $8, $1}' > st.bed
if [ $stbed -eq 1 ]; then exit; fi

## region-specific data
echo generate region-specific data
cat $rt.lst | \
parallel -j${threads} -C' ' 'export f=chr{7}_$(({8}-${flanking}))_$(({8}+${flanking}));\
      awk "(\$7==chr && \$8 >= pos-s && \$8 <= pos+s){if(\$2<\$3) {a1=\$2; a2=\$3;} else {a1=\$3; a2=\$2};\
           \$0=\$0 \" \" \$7 \":\" \$8 \"_\" a1 \"_\" a2;print}" chr={7} pos={8} s=${flanking} $rt.input |\
           sort -k9,9 > $f.dat'
echo "--> map/ped"
ls $GEN_location/chr*.gen | \
xargs -e -n1 basename | \
sed 's/\.gen//g' | \
parallel -j${threads} --env wd -C' ' '\
      awk -f ${FM_location}/files/order.awk {}.gen > {}.ord;\
      gtool -G --g {}.ord --s ${sample_file} --ped {}.ped --map {}.map \
           --missing 0.05 --threshold 0.9 --log {}.log --snp --alleles --chr $(echo {}| \
             cut -d"_" -f1|sed "s/chr//g")'
echo "--> GWAS .sumstats auxiliary files"
ls $GEN_location/*.gen | \
sed 's/\.gen//g'| \
parallel -j${threads} -C' ' 'sort -k2,2 {}.map| \
     join -19 -22 {}.dat -| \
     sort -k10,10>{}.incl'
cat st.bed | \
parallel -j${threads} --env wd -C' ' 'f=chr{1}_{2}_{3}; \
     awk "{print \$8,\$9,\$3,\$4,\$5,\$6,\$7,\$2,\$1,\$5/\$6}" $f.incl > $f.r; \
     cut -d" " -f9,10 $f.r2>$f.z; \
     awk "{print \$1}" $f.incl > $f.inc; \
     awk "{print \$1,\$4,\$3,\$13,\$14}" $f.incl > $f.a; \
     echo "RSID position chromosome A_allele B_allele" > $f.incl_variants; \
     awk "{print \$1,\$9,\$8,\$4,\$3}" $f.incl >> $f.incl_variants'

## finemapping
echo "--> bfile"
rm *bed *bim *fam
if [ ${sample_to_exclude} == "" ]; then 
   export OPTs=""
else 
   export OPTs="--remove ${sample_to_exclude}"
fi
ls $GEN_location/chr*.gen | \
awk '(gsub(/\.gen/,""))' | \
parallel -j${threads} --env wd -C' ' '\
         plink-1.9 --file {} --missing-genotype N --extract {}.inc ${OPTs} \
         --make-bed --keep-allele-order --a2-allele {}.a 3 1 --out {}'
echo "JAM, IPD"
ls $GEN_location/chr*.gen | \
xargs -e -n1 basename | \
awk '(gsub(/\.gen/,""))' | \
parallel -j${threads} -C' ' '\
         plink-1.9 --bfile {} --indep-pairwise 500kb 5 0.80 --maf 0.05 --out {}; \
         grep -w -f {}.prune.in {}.a > {}.p; \
         grep -w -f {}.prune.in {}.dat > {}p.dat; \
         plink-1.9 --bfile {} --extract {}.prune.in --keep-allele-order --a2-allele {}.p 3 1 --make-bed --out {}p'
echo "--> finemap, bcor"
if [ $finemap -eq 1 ]; then
   ls $GEN_location/*.gen | \
   xargs -e -n1 basename | \
   sed 's/\.info//g' | \
   parallel -j3 -C' ' '\
        ldstore --bcor {}.bcor --bplink {} --n-threads ${threads}; \  
        ldstore --bcor {}.bcor --merge ${threads}; \
        ldstore --bcor {}.bcor --matrix {}.ld --incl_variants {}.incl_variants; \
        sed -i -e "s/  */ /g; s/^ *//; /^$/d" {}.ld'
   ls $GEN_location/*.gen | \
   xargs -e -n1 basename | \
   sed 's/\.info//g' | \
   parallel -j3 -C' ' '\ 
         grep -w -f {}.prune.in {}.z > {}p.z; \
         ldstore --bcor {}p.bcor --bplink {}p --n-threads ${threads}; \
         ldstore --bcor {}p.bcor --merge ${threads}; \
         ldstore --bcor {}p.bcor --matrix {}p.ld; \
         sed -i -e "s/  */ /g; s/^ *//; /^$/d" {}p.ld'
   echo "z;ld;snp;config;log;n-ind" > finemap.cfg
   cat st.bed | \
   parallel -j${threads} -C ' ' '\
         f=chr{1}_{2}_{3}; \
         sort -k7,7n $f.r | \
         tail -n1|cut -d" " -f7| \
         awk -vf=$f "{print sprintf(\"%s.z;%s.ld;%s.snp;%s.config;%s.log;%d\",f,f,f,f,f,int(\$1))}" >> finemap.cfg'
   finemap --sss --in-files finemap.cfg --n-causal-max 5 --corr-config 0.9
   cat st.bed | \
   parallel -j${threads} -C' ' '\
         f=chr{1}_{2}_{3};R --no-save < ${FM_location}/files/finemap-check.R > $f.chk'
   sed 's/\./p\./g' finemap.cfg > finemapp.cfg
   finemap --sss --in-files finemapp.cfg --n-causal-max 5 --corr-config 0.9
   cat st.bed | \
   parallel -j${threads} -C' ' 'f=chr{1}_{2}_{3}p;R --no-save < ${FM_location}/files/finemap-check.R > $f.chk'
fi
echo "--> JAM"
if [ $JAM -eq 1 ]; then
   cat st.bed | parallel -j${threads} -C' ' 'export fp=chr{1}_{2}_{3}p; R --no-save <${FM_location}/files/JAM.R > ${fp}.log'
fi
if [ $CAVIAR -eq 1 ]; then
   cat st.bed | parallel -j${threads} -C' ' '\
       export f=chr{1}_{2}_{3};CAVIAR -z ${f}.z -l ${f}.ld -r 0.9 -o ${f}'
fi
if [ $CAVIARBF -eq 1 ]; then
   cat st.bed | parallel -j${threads} -C' ' '\
       export f=chr{1}_{2}_{3};caviarbf -z ${f}.z -r ${f}.ld -n $N -t 0 -a 0.1 -c 3 --appr -o ${f}.caviarbf'
fi
