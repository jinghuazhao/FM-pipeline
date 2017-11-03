#!/bin/bash
# 3-11-2017 MRC-Epid JHZ

## settings -- change as apporopriate
# working directory
export wd=/genetics/data/gwas/1-11-17
# GWAS summary statistics (the .sumstats file)
export args=$1
# filename containing list of lead SNPs
export snplist=$wd/1.snps
# GEN files, named chr{chr}_{start}_{end}.gen
export GEN_location=/genetics/data/gwas/6-7-17/HRC
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
export LocusZoom=1
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
ln -sf $wd/$args
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
echo Supplement .sumstats with chromosomal positions
export rt=$dir/$(basename $args)
awk '{
  $2=toupper($2)
  $3=toupper($3)
};1' $args | \
sort -k1,1 | \
join -11 -23 - snp150.txt | \
sed 's/chr//g' > $rt.input
sort -k1,1 ${snplist} | \
join $dir/$(basename $args).input - > $rt.lst
echo "chr start end pos rsid" > st.bed
grep -w -f ${snplist} $rt.input | \
awk -vs=${flanking} '{print $7,$8-s,$8+s, $8, $1}' >> st.bed
if [ $stbed -eq 1 ]; then
   echo "st.bed is generated"
   exit
fi

## region-specific data
echo Generate region-specific data
cat $rt.lst | \
parallel -j${threads} -C' ' 'export f=chr{7}_$(({8}-${flanking}))_$(({8}+${flanking}));\
    awk "(\$7==chr && \$8 >= pos-s && \$8 <= pos+s){if(\$2<\$3) {a1=\$2; a2=\$3;} else {a1=\$3; a2=\$2};\
         \$0=\$0 \" \" \$7 \":\" \$8 \"_\" a1 \"_\" a2;print}" chr={7} pos={8} s=${flanking} $rt.input |\
         sort -k9,9 > $f.dat'
echo "--> map/ped"
awk 'NR>1' st.bed | \
parallel -j${threads} --env wd -C' ' '\
    export f=chr{1}_{2}_{3}; \
    awk -f ${FM_location}/files/order.awk $GEN_location/$f.gen > $f.ord;\
    gtool -G --g $f.ord --s ${sample_file} --ped $f.ped --map $f.map \
         --missing 0.05 --threshold 0.9 --log $f.log --snp --alleles --chr $(echo $f| \
         cut -d"_" -f1|sed "s/chr//g")'
echo "--> GWAS .sumstats auxiliary files"
awk 'NR>1' st.bed | \
parallel -j${threads} -C' ' '\
    export f=chr{1}_{2}_{3}; \
    sort -k2,2 $f.map | \
    join -19 -22 $f.dat - | \
    sort -k10,10 > $f.incl'
awk 'NR>1' st.bed | \
parallel -j${threads} --env wd -C' ' '\
    export f=chr{1}_{2}_{3}; \
    awk "{print \$8,\$9,\$3,\$4,\$5,\$6,\$7,\$2,\$1,\$5/\$6}" $f.incl > $f.r; \
    cut -d" " -f9,10 $f.r > $f.z; \
    awk "{print \$1}" $f.incl > $f.inc; \
    awk "{print \$1,\$4,\$3,\$13,\$14}" $f.incl > $f.a; \
    echo "RSID position chromosome A_allele B_allele" > $f.incl_variants; \
    awk "{print \$1,\$9,\$8,\$4,\$3}" $f.incl >> $f.incl_variants'

## finemapping
echo "--> bfile"
if [ ${sample_to_exclude} == "" ]; then 
   export OPTs=""
else 
   export OPTs="--remove ${sample_to_exclude}"
fi
awk 'NR>1' st.bed | \
parallel -j${threads} -C' ' '\
    export f=chr{1}_{2}_{3}; \
    rm -f $f.bed $f.bim $f.fam
    plink-1.9 --file $f --missing-genotype N --extract $f.inc ${OPTs} \
    --make-bed --keep-allele-order --a2-allele $f.a 3 1 --out $f'
echo "JAM, IPD"
awk 'NR>1' st.bed | \
parallel -j${threads} -C' ' '\
    export f=chr{1}_{2}_{3}; \
    plink-1.9 --bfile $f --indep-pairwise 500kb 5 0.80 --maf 0.05 --out $f; \
    grep -w -f $f.prune.in $f.a > $f.p; \
    grep -w -f $f.prune.in $f.dat > ${f}p.dat; \
    plink-1.9 --bfile $f --extract $f.prune.in --keep-allele-order --a2-allele $f.p 3 1 --make-bed --out ${f}p'
echo "--> finemap, bcor"
if [ $finemap -eq 1 ]; then
   awk 'NR>1' st.bed | \
   parallel -j3 -C' ' '\
       export f=chr{1}_{2}_{3}; \
       ldstore --bcor $f.bcor --bplink $f --n-threads ${threads}; \  
       ldstore --bcor $f.bcor --merge ${threads}; \
       ldstore --bcor $f.bcor --matrix $f.ld --incl_variants $f.incl_variants; \
       sed -i -e "s/  */ /g; s/^ *//; /^$/d" $f.ld'
   awk 'NR>1' st.bed | \
   parallel -j3 -C' ' '\ 
       export f=chr{1}_{2}_{3}; \
       grep -w -f $f.prune.in $f.z > ${f}p.z; \
       ldstore --bcor ${f}p.bcor --bplink ${f}p --n-threads ${threads}; \
       ldstore --bcor ${f}p.bcor --merge ${threads}; \
       ldstore --bcor ${f}p.bcor --matrix ${f}p.ld; \
       sed -i -e "s/  */ /g; s/^ *//; /^$/d" ${f}p.ld'
   echo "z;ld;snp;config;log;n-ind" > finemap.cfg
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C ' ' '\
       export f=chr{1}_{2}_{3}; \
       sort -k7,7n $f.r | \
       tail -n1|cut -d" " -f7| \
       awk -vf=$f "{print sprintf(\"%s.z;%s.ld;%s.snp;%s.config;%s.log;%d\",f,f,f,f,f,int(\$1))}" >> finemap.cfg'
   finemap --sss --in-files finemap.cfg --n-causal-max 5 --corr-config 0.9
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '\
       export f=chr{1}_{2}_{3}; \
       R --no-save < ${FM_location}/files/finemap-check.R > $f.chk'
   sed 's/\./p\./g' finemap.cfg > finemapp.cfg
   finemap --sss --in-files finemapp.cfg --n-causal-max 5 --corr-config 0.9
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '\
       export f=chr{1}_{2}_{3}p; \
       R --no-save < ${FM_location}/files/finemap-check.R > $f.chk'
fi
echo "--> JAM"
if [ $JAM -eq 1 ]; then
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '\
       export f=chr{1}_{2}_{3}p; \
       R --no-save <${FM_location}/files/JAM.R > ${f}.log'
fi
if [ $CAVIAR -eq 1 ]; then
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '\
       export f=chr{1}_{2}_{3}; \
       CAVIAR -z ${f}.z -l ${f}.ld -r 0.9 -o ${f}'
fi
if [ $CAVIARBF -eq 1 ]; then
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '\
       export f=chr{1}_{2}_{3}; \
       caviarbf -z ${f}.z -r ${f}.ld -n $N -t 0 -a 0.1 -c 3 --appr -o ${f}.caviarbf'
fi
if [ $LocusZoom == 1 ]; then
   echo "{OFS=\"\\t\";if(NR==1) print \"MarkerName\",\"P-value\",\"Weight\"; print \$8,\$6,\$7}" > lz.awk
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '\
       export f=chr{1}_{2}_{3}; \
       export refsnp={5}; \
       awk -f lz.awk $f.r > $f.lz; \
       locuszoom-1.3 --metal $f.lz --refsnp $refsnp --plotonly --no-date;pdftopng $refsnp.pdf -r 300 $refsnp'
fi
if [ $GCTA -eq 1 ]; then
   # setup
   # ma for marginal effects used by GCTA
   rm *cojo *jma *cma
   echo "{if(NR==1) print \"SNP\",\"A1\",\"A2\",\"freq\",\"b\",\"se\",\"p\",\"N\";print \$1,\$4,\$5,\$6,\$7,\$8,\$9,int(\$10)}" > ma.awk
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '\
       export f=chr{1}_{2}_{3}; \
       sort -k4,4 {}_map | \
       join -111 -24 {}.r -|grep -f {}.inc | \
       awk -f ma.awk > {}.ma'
   # --cojo-slct
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '\
       export f=chr{1}_${2}_${3}; \
       gcta64 --bfile $f --cojo-file $f.ma --cojo-slct --out $f;'
   ls *.jma.cojo|sed 's/\.jma\.cojo//g' | \
   parallel -j${threads} -C' ' '\
       echo "SNP Chr bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r rsid" > {}.jma; \
       cut -d" " -f8,9 {}.r | \
       sort -k2,2 | \
       sed "s/ /\t/g">{}.tmp; \
       sort -k2,2 {}.jma.cojo | \
       join -j2 - {}.tmp >> {}.jma'
   echo "region SNP Chr bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r rsid" > gcta-top.csv
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '\
       export f=chr{1}_{2}_{3}; \
       awk "!/SNP/{gsub(/\.jma/,\"\",FILENAME);print FILENAME, \$0}" ${f}.jma >> gcta-slct.csv'
   sed -i 's/ /,/g' gcta-slct.csv

   # --cojo-top-SNPs
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '\
       export f=chr${chr}_${start}_${end}; \
       gcta64 --bfile $f --cojo-file $f.ma --cojo-top-SNPs 3 --out ${f}.top'
   ls *top.jma.cojo | \
   sed 's/\.top\.jma\.cojo//g' | \
   parallel -j${threads} -C' ' '\
       echo "SNP Chr bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r rsid" > {}top.jma; \
       cut -d" " -f8,9 {}.r | \
       sort -k2,2 | \
       sed "s/ /\t/g">{}.tmp; \
       sort -k2,2 {}.top.jma.cojo | \
       join -j2 - {}.tmp >> {}.top.jma'
   echo "region SNP Chr bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r rsid" > gcta-top.csv
   awk 'NR>1' | \
   parallel -j${threads} -C' ' '\
       export f=chr{1}_{2}_{3}; \
       awk "!/SNP/{gsub(/\.top\.jma/,\"\",FILENAME);print FILENAME, \$0}" ${f}.top.jma >> gcta-top.csv'
   sed -i 's/ /,/g' gcta-top.csv

   # --cojo-cond
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '\
       export f=chr{1}_{2}_{3}; \
       grep {5} ${f}.r | \
       cut -d" " -f11 > ${f}.snpid; \
       gcta64 --bfile $f --cojo-file $f.ma --cojo-cond ${f}.snpid --out ${f}'
   ls *cma.cojo|sed 's/\.cma\.cojo//g' | \
   parallel -j${threads} -C' ' '\
       echo "SNP Chr bp refA freq b se p n freq_geno bC bC_se pC rsid" > {}.cma; \
       cut -d" " -f8,9 {}.r | \
       sort -k2,2 | \
       sed "s/ /\t/g" > {}.tmp; \
       sort -k2,2 {}.cma.cojo | \
       join -j2 - {}.tmp >> {}.cma'
   echo "region SNP Chr bp refA freq b se p n freq_geno bC bC_se pC rsid" > gcta-cond.csv
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '\
       export f=chr{1}_{2}_{3}; \
       awk "!/SNP/{gsub(/\.cma/,\"\",FILENAME);print FILENAME, \$0}" ${f}.cma >> gcta-cond.csv'
   sed -i 's/ /,/g' gcta-cond.csv

   # dosage format
   rt=/gen_omics/data/EPIC-Norfolk/Dosage
   chr=22
   # gcta64 --dosage-mach-gz ${rt}/chr${chr}.dosage.gz ${rt}/chr${chr}.mlinfo.gz --make-grm --thread-num 10 --out chr${chr}
fi

