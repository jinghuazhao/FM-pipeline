#!/bin/bash
# 4-1-2019 JHZ

## SETTINGS

source fmp.ini

## ANALYSIS

if [ $# -lt 1 ] || [ "$args" == "-h" ]; then
    echo "Usage: fmp.sh <input>"
    echo "where <input> is in sumstats format:"
    echo "SNP A1 A2 freqA1 beta se P N chr pos"
    echo "where SNP is RSid, A1 is effect allele"
    exit
fi

export args=$1
if [ $(dirname $args) == "." ]; then dir=$(pwd)/$(basename $args).out; else dir=$args.out; fi
if [ ! -d $dir ]; then mkdir -p $dir; fi
export wd=$(pwd)
cd $dir
ln -sf $wd/$args
export rt=$dir/$(basename $args)
echo "--> $rt.input, st.bed"
awk '{$2=toupper($2);$3=toupper($3)};1' $args > $rt.input
if [ $clumping -eq 1 ]; then
   awk '
   {
      if (NR==1) print "snpid", "P"
      chr=$9;
      pos=$10;
      a1=$2;
      a2=$3;
      if (a1>a2) {
         snpid=chr ":" pos "_" a2 "_" a1;
      } else {
         snpid=chr ":" pos "_" a1 "_" a2;
      }
      print snpid, $7
   }' OFS='\t' $rt.input > $rt.tab
fi
ln -sf $wd/st.bed

echo "--> region-specific finemapping"
awk 'NR==10' st.bed | \
parallel -j${threads} -C' ' \
         --env FM_location \
         --env GEN_location  \
         --env CAVIAR \
         --env CAVIARBF \
         --env FM_summary \
         --env fgwas \
         --env finemap \
         --env GCTA \
         --env JAM \
         --env LD_MAGIC \
         --env LD_PLINK \
         --env LocusZoom \
          '$FM_location/fmp.subs {1} {2} {3} {4} {5}'

exit

# Genome-wide Complex Trait Analysis (GCTA)

if [ $GCTA -eq 1 ]; then
   echo "--> GCTA"
# --cojo-slct <==> jma.cojo, ldr.cojo
   echo "region SNP Chr bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r rsid" > gcta-slct.csv
   ls *.jma.cojo|sed 's/\.jma\.cojo//g' | parallel -j1 -C' ' '
       (
         echo "SNP Chr bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r rsid"
         sort -k2,2 {}.jma.cojo | \
         join -j2 - {}.tmp
       ) > {}.jma'
   awk 'NR>1' st.bed | parallel -j1 -C' ' '
       export f=chr{1}_{2}_{3}
       if [ -f $f.jma ]; then awk "!/SNP/{print f, \$0}" f=$f $f.jma >> gcta-slct.csv; fi'
   sed -i 's/ /,/g' gcta-slct.csv
# --cojo-cond <==> given.cojo, cma.cojo
   echo "region SNP Chr bp refA freq b se p n freq_geno bC bC_se pC rsid" > gcta-cond.csv
   ls *cma.cojo|sed 's/\.cma\.cojo//g' | parallel -j1 -C' ' '
       (
         echo "SNP Chr bp refA freq b se p n freq_geno bC bC_se pC rsid"
         sort -k2,2 {}.cma.cojo | \
         join -j2 - {}.tmp
       ) > {}.cma'
       awk 'NR>1' st.bed | parallel -j1 -C' ' '
       export f=chr{1}_{2}_{3}
       if [ -f $f.cma ]; then awk "!/SNP/{print f, \$0}" f=$f $f.cma >> gcta-cond.csv; fi'
   sed -i 's/ /,/g' gcta-cond.csv
# --cojo-top-SNPs <==> top.jma.cojo, top.ldr.cojo
   echo "region SNP Chr bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r rsid" > gcta-top.csv
   ls *top.jma.cojo | \
   sed 's/\.top\.jma\.cojo//g' | parallel -j1 -C' ' '
       (
         echo "SNP Chr bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r rsid"
         sort -k2,2 {}.top.jma.cojo | \
         join -j2 - {}.tmp
       ) > {}.top.jma'
   awk 'NR>1' st.bed | parallel -j1 -C' ' '
       export f=chr{1}_{2}_{3}
       if [ -f $f.top.jma ]; then awk "!/SNP/{print f, \$0}" f=$f $f.top.jma >> gcta-top.csv; fi'
   sed -i 's/ /,/g' gcta-top.csv
fi

#  Joint Analysis of Marginal SNP Effects (JAM)

if [ $JAM -eq 1 ]; then
   (
     ls *p.sum | sed 's/\p.\sum//g' | parallel -j1 -C' ' 'export f=chr{1}_{2}_{3};awk "NR==2&&\$2>0 {print f}" f=$f ${f}p.sum'
   ) > jam.top
   (
     cat jam.top | parallel -j1 -C' ' 'echo -e "\n" {};cat {}p.top {}p.jam'
   ) > jam.txt
   (
     echo "Region Chr SNP rsid PostProb_model PostProb Median CrI_Lower CrI_Upper Median_Present CrI_Lower_Present CrI_Upper_Present BF"
     ls *sel | parallel -j1 -C' ' '
         awk "NR>1{sub(/\p.sel/,\"\",FILENAME);split(\$2,a,\":\");\$1=FILENAME \" \" a[1];print}"' | \
         sort -k1,1n
   ) > jam.out
   R -q --no-save < ${FM_location}/files/JAM-cs.R > JAM-cs.log
fi

if [ $GCTA -eq 1 ] && [ $JAM -eq 1 ] && [ $finemap -eq 1 ] && [ $clumping -eq 1 ]; then
   R -q --no-save < ${FM_location}/files/gcta-jam-finemap.R > gcta-jam-finemap.log
   sort -k1,1 ld > ld.1
   gunzip -c /gen_omics/data/EPIC-Norfolk/HRC/binary_ped/id3.txt.gz | \
   awk '{snpid=$2;gsub(/:|\_/," ",snpid);split(snpid,a);Chr=a[1];Pos=a[2];print $0,Chr,Pos}' | \
   sort -k1,1 | \
   join -j1 - ld.1 | \
   sort -k4,4n -k5,5n | \
   awk '{print $2}' > ld.dat
   rm ld.1
   plink-1.9 --bfile $bfile --remove $remove_sample --exclude $exclude_snp \
             --extract ld.dat --r2 triangle spaces --out ld
   (
     awk 'NR>1{if(NR==2) printf $8;else printf " " $8}' id | \
     awk '1'
     sed -e 's/[[:space:]]*$//' ld.ld
   ) > ld.mat
   paste -d' ' id ld.mat > ld.txt
   R -q --no-save < ${FM_location}/files/ld.R > ld.log
fi

# functional genomic information with a genome-wide association study (fGWAS)

if [ $fgwas -eq 1 ]; then
   echo "--> fgwas"
   # obtain annotations
   seq 22 | parallel -j${threads} --env fgwas_location_1kg -C' ' '
       if [ ! -f $fgwas_location_1kg/chr{}.gene ]; then
          gunzip -c $fgwas_location_1kg/chr{}.annot.wdist.wcoding.gz | \
          awk "(NR>1){print \$1,\$2,\$3,\$(NF-6),\$(NF-5),\$(NF-2),\$(NF-1),\$NF}" | \
          sort -k2,2 > $fgwas_location_1kg/chr{}.gene
       fi'
   # specify regions
   # -/+ flanking position
   export flanking=250000
   awk -vfl=${flanking} 'NR>1{l=$2;u=$3;print $5,$1,$4,l,u,NR}' st.bed | \
   sort -k1,1 > fgwas.snplist
   cat fgwas.snplist | parallel -j${threads} --env FM_location --env fgwas_location_1kg -C' ' '
       export f=chr{2}_{4}_{5}
       awk -vsn={6} -f $FM_location/files/fgwas.awk $f.r | \
       sort -k3,3 | \
       join -13 -22 - $fgwas_location_1kg/chr{2}.gene | \
       awk -f $FM_location/files/gene.awk | \
       gzip -fc > $f.fgwas.gz'
   # tally for -fine option
   (
   sort -k6,6n fgwas.snplist | \
     parallel -j1 -C' ' '
     gunzip -c chr{2}_{4}_{5}.fgwas.gz | \
     awk "(NR>1 && !/INDEL/)"
   '
   ) > fgwas.tmp
   (
     echo "SNPID CHR POS Z F N ens_coding_exons ens_noncoding_exons tssdist syn nonsyn SEGNUMBER"
     sort -k12,12n -k2,2 -k3,3n fgwas.tmp
   ) > fgwas.fine
   gzip -f fgwas.fine
   # fgwas
   for an in ens_coding_exons ens_noncoding_exons tssdist syn nonsyn;
   do
       fgwas -i fgwas.fine.gz -fine -print -o fgwas-${an} -w ${an}
   done
   fgwas -i fgwas.fine.gz -fine -print -o fgwas -w ens_coding_exons+ens_noncoding_exons+syn+nonsyn
   gunzip -c fgwas.bfs.gz | \
   awk '(NR==1||$10>0.5)' > fgwas.PPA0.5
   # generate table for md document
   awk '(NR>1){print $1}' fgwas.PPA0.5 > fgwas.rsid
   grep -f fgwas.rsid st.bed | \
   awk '{print $5, "*"}' | \
   sort -k1,1 > fgwas.tmp
   (
     head -1 fgwas.PPA0.5 | \
     awk '{gsub(/ /,",",$0);$0=$0 "," "index"};1'
     awk 'NR>1' fgwas.PPA0.5 | \
     sort -k1,1 | \
     join -a1 - fgwas.tmp | \
     sed 's/chr//g' | \
     sort -k2,2n -k3,3n | \
     sed 's/ /,/g'
   ) > fgwas.csv
   # manipulations to trick conditional analysis
   cut -d' ' -f1 fgwas.snplist > fgwas.tmp
   (
     zgrep -f fgwas.tmp -v fgwas.fine.gz | \
     awk '{if(NR==1){print $0,"hit"}else{print $0,0}}'
     zgrep -f fgwas.tmp fgwas.fine.gz | \
     awk '{print $0,1}'
   ) > fgwas.cond
   (
     head -1 fgwas.cond
     awk 'NR>1' fgwas.cond | \
     sort -k12,12n -k3,3
   ) > fgwas.tmp
   awk '{$12="";print}' fgwas.tmp|gzip -cf > fgwas.tmp.gz
   fgwas -i fgwas.tmp.gz -k 500 -print -o fgwas.cond -w ens_coding_exons+ens_noncoding_exons+syn+nonsyn -cond hit
fi

# obsolete with gtool/plink-1.9 handling gen/ped
#   gtool -G --g $GEN_location/$f.ord --s ${sample_file} --ped $GEN_location/$f.ped --map $GEN_location/$f.map \
#         --missing 0.05 --threshold 0.9 --log $f.log --snp --alleles --chr {1}'
#   plink-1.9 --file $GEN_location/$f --missing-genotype N --extract $f.inc \
#         --make-bed --keep-allele-order --a2-allele $f.a 3 1 --out $f'
