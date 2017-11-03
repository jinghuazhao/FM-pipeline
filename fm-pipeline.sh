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
# set to only generate st.bed containg chr, start, end, pos, rsid, r sextuplets
export stbed=0
# number of threads
export threads=5
# software to be included in the analysis; change flags to 1 when available
# the outputs should be available individually
export CAVIAR=0
export CAVIARBF=0
export finemap=1
export JAM=1
export fm_summary=1
export LocusZoom=1
export GCTA=1
export fgwas=0
export fgwas_location_1kg=/genetics/data/software/fgwas/1000-genomes
export FM_location=/genetics/bin/FM-pipeline
if [ $# -lt 1 ] || [ "$args" == "-h" ]; then
    echo "Usage: fm-pipeline.sh <input>"
    echo "where <input> is in sumstats format:"
    echo "SNP A1 A2 freqA1 beta se P N"
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
grep -w -f ${snplist} $rt.input | \
awk -vs=${flanking} '{print $9,$10-s,$10+s, $10, $1}' > st.dat
echo "chr start end pos rsid r" > st.bed
awk '{$0=$0 " " NR};1' st.dat >> st.bed
rm st.dat
if [ $stbed -eq 1 ]; then
   echo "st.bed is generated"
   exit
fi

## region-specific data
echo Generate region-specific data
cat $rt.lst | \
parallel -j${threads} -C' ' '\
    export f=chr{9}_$(({10}-${flanking}))_$(({10}+${flanking}));\
    awk "(\$9==chr && \$10 >= pos-s && \$10 <= pos+s){if(\$2<\$3) {a1=\$2; a2=\$3;} else {a1=\$3; a2=\$2};\
         \$0=\$0 \" \" \$9 \":\" \$10 \"_\" a1 \"_\" a2;print}" chr={9} pos={10} s=${flanking} $rt.input |\
         sort -k11,11 > $f.dat'
echo "--> map/ped"
awk 'NR>1' st.bed | \
parallel -j${threads} --env wd -C' ' '\
    export f=chr{1}_{2}_{3}; \
    awk -f ${FM_location}/files/order.awk $GEN_location/$f.gen > $GEN_location/$f.ord;\
    gtool -G --g $GEN_location/$f.ord --s ${sample_file} --ped $GEN_location/$f.ped --map $GEN_location/$f.map \
         --missing 0.05 --threshold 0.9 --log $f.log --snp --alleles --chr {1}'
echo "--> GWAS .sumstats auxiliary files"
awk 'NR>1' st.bed | \
parallel -j${threads} -C' ' '\
    export f=chr{1}_{2}_{3}; \
    sort -k2,2 $GEN_location/$f.map | \
    join -111 -22 $f.dat - | \
    sort -k11,11 > $f.incl'
awk 'NR>1' st.bed | \
parallel -j${threads} --env wd -C' ' '\
    export f=chr{1}_{2}_{3}; \
    awk "{print \$10,\$11,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$2,\$1,\$6/\$7}" $f.incl > $f.r; \
    cut -d" " -f11,12 $f.r > $f.z; \
    awk "{print \$1}" $f.incl > $f.inc; \
    awk "{print \$1,\$4,\$3,\$15,\$16}" $f.incl > $f.a; \
    echo "RSID position chromosome A_allele B_allele" > $f.incl_variants; \
    awk "{print \$1,\$11,\$10,\$4,\$3}" $f.incl >> $f.incl_variants'

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
    plink-1.9 --file $GEN_location/$f --missing-genotype N --extract $f.inc ${OPTs} \
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
       sort -k9,9g $f.r | \
       tail -n1|cut -d" " -f9| \
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
       R --no-save < ${FM_location}/files/JAM.R > ${f}.log'
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
       caviarbf -z ${f}.z -r ${f}.ld -n $(sort -k9,9g $f.r | \
       tail -n1 | cut -d" " -f9) -t 0 -a 0.1 -c 3 --appr -o ${f}.caviarbf'
fi
if [ $LocusZoom -eq 1 ]; then
   echo "{OFS=\"\\t\";if(NR==1) print \"MarkerName\",\"P-value\",\"Weight\"; print \$10,\$8,\$9}" > lz.awk
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '\
       export f=chr{1}_{2}_{3}; \
       export refsnp={5}; \
       awk -f lz.awk $f.r > $f.lz; \
       locuszoom-1.3 --metal $f.lz --refsnp $refsnp --flank 250kb --plotonly \
                     --source 1000G_March2012 --build hg19 --pop EUR --no-date; \
       pdftopng $refsnp.pdf -r 300 $refsnp'
fi
if [ $fm_summary -eq 1 ]; then
   echo "region chr pos A B Freq1 Effect StdErr P N SNP inCredible probNorm cumSum" | \
   sed 's/ /\t/g' > FM-summary.txt
   awk 'NR>1' st.bed | \
   parallel -j${threads} --env FM_location -C' ' '\
       export f=chr{1}_{2}_{3}; \
       $FM_location/files/getCredible.r {6}; \
       awk "!/SNP/{gsub(/\.cre/,\"\",FILENAME);print FILENAME, \$0}" OFS="\t" $f.cre >> FM-summary.txt'
fi
if [ $GCTA -eq 1 ]; then
   echo "--> GCTA"
   awk 'NR>1' st.bed | \
   parallel -j${threads} --env GEN_location -C' ' '\
       export f=chr{1}_{2}_{3}; \
       awk -vchr={1} -f $FM_location/files/info.awk $GEN_location/$f.info | \
       sort -k2,2 > $f.tmp; \
       sort -k2,2 $GEN_location/$f.map | \
       join -j2 $f.tmp - | \
       awk -vOFS="\t" "{print \$8,\$7,0,\$1,\$11,\$12,\$3}" > ${f}_map'
   # ma for marginal effects used by GCTA
   rm -f *cojo *jma *cma
   echo "{if(NR==1) print \"SNP\",\"A1\",\"A2\",\"freq\",\"b\",\"se\",\"p\",\"N\";\
         print \$1,\$4,\$5,\$6,\$7,\$8,\$9,int(\$10)}" > ma.awk
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '\
       export f=chr{1}_{2}_{3}; \
       sort -k4,4 ${f}_map | \
       join -111 -24 $f.r - | \
       grep -f $f.inc | \
       awk -f ma.awk > $f.ma'
   # --cojo-slct
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '\
       export f=chr{1}_{2}_{3}; \
       gcta64 --bfile $f --cojo-file $f.ma --cojo-slct --out $f'
   ls *.jma.cojo|sed 's/\.jma\.cojo//g' | \
   parallel -j${threads} -C' ' '\
       echo "SNP Chr bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r rsid" > {}.jma; \
       cut -d" " -f10,11 {}.r | \
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
       export f=chr{1}_{2}_{3}; \
       gcta64 --bfile $f --cojo-file $f.ma --cojo-top-SNPs 3 --out ${f}.top'
   ls *top.jma.cojo | \
   sed 's/\.top\.jma\.cojo//g' | \
   parallel -j${threads} -C' ' '\
       echo "SNP Chr bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r rsid" > {}top.jma; \
       cut -d" " -f10,11 {}.r | \
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
       cut -d" " -f10,11 {}.r | \
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
   # rt=/gen_omics/data/EPIC-Norfolk/Dosage
   # chr=22
   # gcta64 --dosage-mach-gz ${rt}/chr${chr}.dosage.gz ${rt}/chr${chr}.mlinfo.gz --make-grm --thread-num 10 --out chr${chr}
fi
if [ $fgwas -eq 1 ]; then
  echo "--> fgwas"
  if [ ! -d "fgwas" ]; then
     mkdir -p fgwas
  fi
  cd fgwas
  # obtain annotations
  seq 22 | \
  parallel -j${threads} -C' ' '\
      if [ ! -f $fgwas_location_1kg/chr{}.gen ]; then
      gunzip -c $fgwas_location_1kg/chr{}.annot.wdist.wcoding.gz | \
      awk "(NR>1){\$2=sprintf(\"%010d\",\$2); print \$1,\$2,\$3,\$(NF-6),\$(NF-5),\$(NF-2),\$(NF-1),\$NF}"|\
      sort -k2,2 > $fgwas_location_1kg/chr{}.gene \
      fi'

  # specify regions
  awk -vfl=${flanking} '{l=$2;u=$3;if(l<0) l=1;print $5,$1,$4,l,u,NR}' st.bed | \
  sort -k1,1 > fgwas.snplist

  # the standard fgwas data
  awk 'NR>1' fgwas.snplist | \
  parallel -j${threads} --env GEN_location -C' ' '\
  read rsid chr pos start end sn <<<$(awk -vline={} "NR==line" fgwas.snplist);\
      export f=chr{2}_{3}_{4};\
      awk -vsn={6} -f $GEN_location/files/fgwas.awk $f.r | \
      join -13 -22 - $fgwas_location_1kg/chr{1}.gene | \
      awk -f $GEN_location/files/gene.awk | \
      gzip -fc > $f.fgwas.gz'

  # tally for -fine option
  echo "SNPID CHR POS Z F N ens_coding_exons ens_noncoding_exons tssdist syn nonsyn SEGNUMBER"> fgwas.fine
  sort -k6,6n fgwas.snplist|awk '{
    cmd=sprintf("gunzip -c chr%d_%d_%d.fgwas.gz | awk \x27NR>1\x27|awk \x27!/INDEL/\x27>> fgwas.fine",$2,$4,$5,$6)
    system(cmd)
  }'
  gzip -f fgwas.fine

  # f-gwas
  for an in ens_coding_exons ens_noncoding_exons tssdist syn nonsyn; do
     fgwas -i fgwas.fine.gz -fine -print -o fgwas-${an} -w ${an}
  done
  fgwas -i fgwas.fine.gz -fine -print -o fgwas -w ens_coding_exons+ens_noncoding_exons+syn+nonsyn
  gunzip -c fgwas.bfs.gz | awk '(NR==1||$10>0.5)' > fgwas.PPA0.5

  # generate table for md document
  awk '(NR>1){print $1}' fgwas.PPA0.5 > fgwas.rsid
  grep -f fgwas.rsid $wd/doc/region.dat|awk '{print $1, "*"}' |sort -k1,1 > fgwas.tmp
  head -1 fgwas.PPA0.5 | awk '{gsub(/ /,",",$0);$0=$0 "," "index"};1' > fgwas.csv
  awk 'NR>1' fgwas.PPA0.5 | \
  sort -k1,1 | \
  join -a1 - fgwas.tmp | \
  sed 's/chr//g' | \
  sort -k2,2n -k3,3n | \
  sed 's/ /,/g' >> fgwas.csv

  # manipulations to trick conditional analysis
  cut -d' ' -f1 fgwas.snplist > fgwas.tmp
  zgrep -f fgwas.tmp -v fgwas.fine.gz|awk '{if(NR==1){print $0,"hit"}else{print $0,0}}' > fgwas.cond
  zgrep -f fgwas.tmp fgwas.fine.gz | \
  awk '{print $0,1}' >> fgwas.cond
  head -1 fgwas.cond > fgwas.tmp
  awk 'NR>1' fgwas.cond|sort -k12,12n -k3,3 >> fgwas.tmp
  awk '{$12="";print}' fgwas.tmp|gzip -cf > fgwas.tmp.gz
  fgwas -i fgwas.tmp.gz -k 500 -print -o fgwas.cond -w ens_coding_exons+ens_noncoding_exons+syn+nonsyn -cond hit
fi
