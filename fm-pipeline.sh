#!/bin/bash
# 9-11-2017 MRC-Epid JHZ

## settings -- change as apporopriate
# working directory
export wd=/genetics/data/gwas/1-11-17
# GWAS summary statistics (the .sumstats file)
export args=$1
# filename containing list of lead SNPs
export snplist=$wd/97.snps
# GEN files, named chr{chr}_{start}_{end}.gen
export GEN_location=/genetics/data/gwas/6-7-17/HRC
# sample file
export sample_file=/gen_omics/data/EPIC-Norfolk/HRC/EPIC-Norfolk.sample
# sample for exclusion
export sample_to_exclude=$wd/exclude.dat
# -/+ flanking position
export flanking=250000
# to generate st.bed containg chr, start, end, pos, rsid, r sextuplets
export use_UCSC=0
# regenerate ped/map instead of using existing ones
export gen_to_ped=0
# force to use rsid with finemap
export force_rsid=0
# number of threads
export threads=5
# also to obtain results after LD pruning
export allow_prune=0
# software to be included in the analysis; change flags to 1 when available
# the outputs should be available individually
export CAVIAR=0
export CAVIARBF=0
export FM_summary=0
export GCTA=0
export JAM=0
export LocusZoom=1
export fgwas=0
export finemap=1
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
   dir=$(pwd)/$(basename $args).out
else
   dir=$args.out
fi
if [ ! -d $dir ]; then
   mkdir -p $dir
fi
cd $dir
ln -sf $wd/$args
export rt=$dir/$(basename $args)
if [ $use_UCSC -eq 1 ]; then
   echo Supplement .sumstats with chromosomal positions
   if $(test -f ${FM_location}/snp150.txt ); then
      echo "Chromosomal positions are ready to use"
      ln -sf ${FM_location}/snp150.txt
   else
      echo "Obtaining chromosomal positions"
      wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp150.txt.gz
      gunzip -c snp150.txt.gz | \
      awk '{split($2,a,"_");sub(/chr/,"",a[1]);print a[1],$4,$5}' | \
      sort -k3,3 > snp150.txt
   fi
   awk '{
     $2=toupper($2)
     $3=toupper($3)
   };1' $args | \
   sort -k1,1 | \
   join -11 -23 - snp150.txt > $rt.input
   grep -w -f ${snplist} $rt.input | \
   awk -vs=${flanking} '{
      l=$10-s
      if(l<0) l=1
      print $9, l, $10+s, $10, $1
   }' > st.dat
   echo "chr start end pos rsid r" > st.bed
   sort -k1,1n -k4,4n st.dat | \
   awk '{$0=$0 " " NR};1' >> st.bed
   rm st.dat
else
   awk '{
     $2=toupper($2)
     $3=toupper($3)
   };1' $args > $rt.input
   ln -sf $wd/st.bed
fi

echo Generate region-specific data
if [ $use_UCSC -eq 1 ]; then
   sort -k1,1 ${snplist} | \
   join $dir/$(basename $args).input - > $rt.lst
   cat $rt.lst | \
   parallel -j${threads} -C' ' '
       export l=$(({10}-${flanking})); \
       if [ $l -le 0 ]; then export l=1; fi; \
       export u=$(({10}+${flanking})); \
       export f=chr{9}_${l}_${u}; \
       awk "(\$9==chr && \$10 >= l && \$10 <= u){if(\$2<\$3) {a1=\$2; a2=\$3;} else {a1=\$3; a2=\$2};\
            \$0=\$0 \" \" \$9 \":\" \$10 \"_\" a1 \"_\" a2;print}" chr={9} l=$l u=$u $rt.input |\
            sort -k11,11 > $f.dat'
else
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}; \
       awk "(\$9==chr && \$10 >= l && \$10 <= u){if(\$2<\$3) {a1=\$2; a2=\$3;} else {a1=\$3; a2=\$2};\
            \$0=\$0 \" \" \$9 \":\" \$10 \"_\" a1 \"_\" a2;print}" chr={1} l={2} u={3} $rt.input |\
            sort -k11,11 > $f.dat'
fi
if [ $gen_to_ped -eq 1 ]; then
   echo "--> map/ped"
   awk 'NR>1' st.bed | \
   parallel -j${threads} --env wd -C' ' '
       export f=chr{1}_{2}_{3}; \
       awk -f ${FM_location}/files/order.awk $GEN_location/$f.gen > $GEN_location/$f.ord;\
       gtool -G --g $GEN_location/$f.ord --s ${sample_file} --ped $GEN_location/$f.ped --map $GEN_location/$f.map \
            --missing 0.05 --threshold 0.9 --log $f.log --snp --alleles --chr {1}'
fi
echo "--> GWAS .sumstats auxiliary files"
awk 'NR>1' st.bed | \
parallel -j${threads} -C' ' '
    export f=chr{1}_{2}_{3}; \
    sort -k2,2 $GEN_location/$f.map | \
    join -111 -22 $f.dat - | \
    sort -k11,11 > $f.incl; \
    cut -d" " -f10,11 $f.r > $f.rsid; \
    awk "{print \$10,\$11,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$2,\$1,\$6/\$7}" $f.incl > $f.r'
if [ $force_rsid -eq 1 ]; then
   awk 'NR>1' st.bed | \
   parallel -j${threads} --env wd -C' ' '
       export f=chr{1}_{2}_{3}; \
       cut -d" " -f10,12 $f.r > ${f}rsid.z; \
       awk "{print \$2,\$4,\$3,\$15,\$16}" $f.incl > ${f}rsid.a; \
       echo "RSID position chromosome A_allele B_allele" > ${f}rsid.incl_variants; \
       awk "{print \$2,\$11,\$10,\$4,\$3}" $f.incl >> ${f}rsid.incl_variants'
else
   awk 'NR>1' st.bed | \
   parallel -j${threads} --env wd -C' ' '
       export f=chr{1}_{2}_{3}; \
       cut -d" " -f11,12 $f.r > $f.z; \
       awk "{print \$1}" $f.incl > $f.inc; \
       awk "{print \$1,\$4,\$3,\$15,\$16}" $f.incl > $f.a; \
       echo "RSID position chromosome A_allele B_allele" > $f.incl_variants; \
       awk "{print \$1,\$11,\$10,\$4,\$3}" $f.incl >> $f.incl_variants'
fi

## finemapping
echo "--> bfile"
if [ ${sample_to_exclude} == "" ]; then 
   export OPTs=""
else 
   export OPTs="--remove ${sample_to_exclude}"
fi
if [ $force_rsid -eq 1 ]; then
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}; \
       plink-1.9 --file $GEN_location/$f --missing-genotype N --extract $f.inc ${OPTs} \
       --make-bed --keep-allele-order --a2-allele $f.a 3 1 --update-name $f.rsid --out ${f}rsid'
else
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}; \
       plink-1.9 --file $GEN_location/$f --missing-genotype N --extract $f.inc ${OPTs} \
       --make-bed --keep-allele-order --a2-allele $f.a 3 1 --out $f'
fi
if [ $allow_prune -eq 1 ]; then
   echo "JAM, IPD"
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}; \
       plink-1.9 --bfile $f --indep-pairwise 500kb 5 0.80 --maf 0.05 --out $f; \
       grep -w -f $f.prune.in $f.a > $f.p; \
       grep -w -f $f.prune.in $f.dat > ${f}p.dat; \
       plink-1.9 --bfile $f --extract $f.prune.in --keep-allele-order --a2-allele $f.p 3 1 --make-bed --out ${f}p'
fi

if [ $CAVIAR -eq 1 ]; then
   echo "--> CAVIAR"
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}; \
       CAVIAR -z $f.z -l $f.ld -r 0.9 -o $f'
fi

if [ $CAVIARBF -eq 1 ]; then
   echo "--> CAVIARBF"
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}; \
       caviarbf -z $f.z -r $f.ld -n $(sort -k9,9g $f.r | \
       tail -n1 | cut -d" " -f9) -t 0 -a 0.1 -c 3 --appr -o $f.caviarbf'
fi

if [ $FM_summary -eq 1 ]; then
   ecoh "--> FM-summary"
   echo "region chr pos A B Freq1 Effect StdErr P N SNP inCredible probNorm cumSum" | \
   sed 's/ /\t/g' > FM-summary.txt
   awk 'NR>1' st.bed | \
   parallel -j${threads} --env FM_location -C' ' '
       export f=chr{1}_{2}_{3}; \
       $FM_location/files/getCredible.r {6}; \
       awk "!/SNP/{print ENVIRON[\"f\"], \$0}" OFS="\t" $f.cre >> FM-summary.txt'
fi

if [ $GCTA -eq 1 ]; then
   echo "--> GCTA"
   awk 'NR>1' st.bed | \
   parallel -j${threads} --env GEN_location -C' ' '
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
   parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}; \
       sort -k4,4 ${f}_map | \
       join -111 -24 $f.r - | \
       grep -f $f.inc | \
       awk -f ma.awk > $f.ma'
   # --cojo-slct
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}; \
       gcta64 --bfile $f --cojo-file $f.ma --cojo-slct --out $f'
   ls *.jma.cojo|sed 's/\.jma\.cojo//g' | \
   parallel -j${threads} -C' ' '
       echo "SNP Chr bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r rsid" > {}.jma; \
       cut -d" " -f10,11 {}.r | \
       sort -k2,2 | \
       sed "s/ /\t/g">{}.tmp; \
       sort -k2,2 {}.jma.cojo | \
       join -j2 - {}.tmp >> {}.jma'
   echo "region SNP Chr bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r rsid" > gcta-slct.csv
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}; \
       awk "!/SNP/{print ENVIRON[\"f\"], \$0}" $f.jma >> gcta-slct.csv'
   sed -i 's/ /,/g' gcta-slct.csv
   # --cojo-top-SNPs
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}; \
       gcta64 --bfile $f --cojo-file $f.ma --cojo-top-SNPs 3 --out $f.top'
   ls *top.jma.cojo | \
   sed 's/\.top\.jma\.cojo//g' | \
   parallel -j${threads} -C' ' '
       echo "SNP Chr bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r rsid" > {}.top.jma; \
       cut -d" " -f10,11 {}.r | \
       sort -k2,2 | \
       sed "s/ /\t/g" > {}.tmp; \
       sort -k2,2 {}.top.jma.cojo | \
       join -j2 - {}.tmp >> {}.top.jma'
   echo "region SNP Chr bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r rsid" > gcta-top.csv
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}; \
       awk "!/SNP/{print ENVIRON[\"f\"], \$0}" $f.top.jma >> gcta-top.csv'
   sed -i 's/ /,/g' gcta-top.csv
   # --cojo-cond
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}; \
       grep {5} $f.r | \
       cut -d" " -f11 > $f.snpid; \
       gcta64 --bfile $f --cojo-file $f.ma --cojo-cond $f.snpid --out $f'
   ls *cma.cojo|sed 's/\.cma\.cojo//g' | \
   parallel -j${threads} -C' ' '
       echo "SNP Chr bp refA freq b se p n freq_geno bC bC_se pC rsid" > {}.cma; \
       cut -d" " -f10,11 {}.r | \
       sort -k2,2 | \
       sed "s/ /\t/g" > {}.tmp; \
       sort -k2,2 {}.cma.cojo | \
       join -j2 - {}.tmp >> {}.cma'
   echo "region SNP Chr bp refA freq b se p n freq_geno bC bC_se pC rsid" > gcta-cond.csv
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}; \
       awk "!/SNP/{print ENVIRON[\"f\"], \$0}" $f.cma >> gcta-cond.csv'
   sed -i 's/ /,/g' gcta-cond.csv
fi

if [ $JAM -eq 1 ]; then
   echo "--> JAM"
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}p; \
       R --no-save < ${FM_location}/files/JAM.R > $f.log'
fi

if [ $LocusZoom -eq 1 ]; then
   echo "--> LocusZoom"
   awk 'NR>1' st.bed | \
   parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}; \
       awk "{OFS=\"\\t\";if(NR==1) print \"MarkerName\",\"P-value\",\"Weight\"; print \$10,\$8,\$9}" $f.r > $f.lz'
   awk 'NR>1' st.bed | \
   parallel -j1 -C' ' '
       rm -f ld_cache.db; \
       locuszoom-1.4 --source 1000G_March2012 --build hg19 --pop EUR --metal chr{1}_{2}_{3}.lz --plotonly --chr {1} --start {2} --end {3} --no-date; \
       pdftopng chr{1}_{2}-{3}.pdf -r 300 {5}'
fi

if [ $fgwas -eq 1 ]; then
   echo "--> fgwas"
   # obtain annotations
   seq 22 | \
   parallel -j${threads} -C' ' '
       if [ ! -f $fgwas_location_1kg/chr{}.gen ]; then
       gunzip -c $fgwas_location_1kg/chr{}.annot.wdist.wcoding.gz | \
       awk "(NR>1){print \$1,\$2,\$3,\$(NF-6),\$(NF-5),\$(NF-2),\$(NF-1),\$NF}" | \
       sort -k2,2 > $fgwas_location_1kg/chr{}.gene \
       fi'
   # specify regions
   awk -vfl=${flanking} '{l=$2;u=$3;if(l<0) l=1;print $5,$1,$4,l,u,NR}' st.bed | \
   sort -k1,1 > fgwas.snplist
   # the standard fgwas data
   awk 'NR>1' fgwas.snplist | \
   parallel -j${threads} --env GEN_location -C' ' '
   read rsid chr pos start end sn <<<$(awk -vline={} "NR==line" fgwas.snplist); \
        export f=chr{2}_{3}_{4}; \
        awk -vsn={6} -f $GEN_location/files/fgwas.awk $f.r | \
        join -13 -22 - $fgwas_location_1kg/chr{1}.gene | \
        awk -f $GEN_location/files/gene.awk | \
        gzip -fc > $f.fgwas.gz'
   # tally for -fine option
   echo "SNPID CHR POS Z F N ens_coding_exons ens_noncoding_exons tssdist syn nonsyn SEGNUMBER" > fgwas.fine
   sort -k6,6n fgwas.snplist | \
   awk '{
     cmd=sprintf("gunzip -c chr%d_%d_%d.fgwas.gz | awk \x27NR>1\x27 | awk \x27!/INDEL/\x27 >> fgwas.fine",$2,$4,$5,$6)
     system(cmd)
   }'
   gzip -f fgwas.fine
   # fgwas
   for an in ens_coding_exons ens_noncoding_exons tssdist syn nonsyn;
   do
       fgwas -i fgwas.fine.gz -fine -print -o fgwas-${an} -w ${an}
   done
   fgwas -i fgwas.fine.gz -fine -print -o fgwas -w ens_coding_exons+ens_noncoding_exons+syn+nonsyn
   gunzip -c fgwas.bfs.gz | 
   awk '(NR==1||$10>0.5)' > fgwas.PPA0.5
   # generate table for md document
   awk '(NR>1){print $1}' fgwas.PPA0.5 > fgwas.rsid
   grep -f fgwas.rsid st.bed | \
   awk '{print $5, "*"}' | \
   sort -k1,1 > fgwas.tmp
   head -1 fgwas.PPA0.5 | \
   awk '{gsub(/ /,",",$0);$0=$0 "," "index"};1' > fgwas.csv
   awk 'NR>1' fgwas.PPA0.5 | \
   sort -k1,1 | \
   join -a1 - fgwas.tmp | \
   sed 's/chr//g' | \
   sort -k2,2n -k3,3n | \
   sed 's/ /,/g' >> fgwas.csv
   # manipulations to trick conditional analysis
   cut -d' ' -f1 fgwas.snplist > fgwas.tmp
   zgrep -f fgwas.tmp -v fgwas.fine.gz | \
   awk '{if(NR==1){print $0,"hit"}else{print $0,0}}' > fgwas.cond
   zgrep -f fgwas.tmp fgwas.fine.gz | \
   awk '{print $0,1}' >> fgwas.cond
   head -1 fgwas.cond > fgwas.tmp
   awk 'NR>1' fgwas.cond | \
   sort -k12,12n -k3,3 >> fgwas.tmp
   awk '{$12="";print}' fgwas.tmp|gzip -cf > fgwas.tmp.gz
   fgwas -i fgwas.tmp.gz -k 500 -print -o fgwas.cond -w ens_coding_exons+ens_noncoding_exons+syn+nonsyn -cond hit
fi

if [ $finemap -eq 1 ]; then
   echo "--> finemap"
   if [ $force_rsid -eq 1 ]; then
      awk 'NR>1' st.bed | \
      parallel -j${threads} -C' ' '
          export f=chr{1}_{2}_{3}; \
          ldstore --bcor ${f}rsid.bcor --bplink ${f}rsid --n-threads ${threads}; \  
          ldstore --bcor ${f}rsid.bcor --merge ${threads}; \
          ldstore --bcor ${f}rsid.bcor --matrix ${f}rsid.ld --incl_variants ${f}rsid.incl_variants; \
          sed -i -e "s/  */ /g; s/^ *//; /^$/d" ${f}rsid.ld'
      echo "z;ld;snp;config;log;n-ind" > finemap_rsid.cfg
      awk 'NR>1' st.bed | \
      parallel -j${threads} -C ' ' '
          export f=chr{1}_{2}_{3}; \
          sort -k9,9g $f.r | \
          tail -n1|cut -d" " -f9 | \
          awk -vf=$f "{print sprintf(\"%srsid.z;%srsid.ld;%srsid.snp;%srsid.config;%srsid.log;%d\",f,f,f,f,f,int(\$1))}" >> finemap_rsid.cfg'
      finemap --sss --in-files finemap_rsid.cfg --n-causal-max 5 --corr-config 0.9
   elif [ $allow_prune -eq 1 ]; then
      awk 'NR>1' st.bed | \
      parallel -j${threads} -C' ' '
          export f=chr{1}_{2}_{3}; \
          grep -w -f $f.prune.in $f.z > ${f}p.z; \
          ldstore --bcor ${f}p.bcor --bplink ${f}p --n-threads ${threads}; \
          ldstore --bcor ${f}p.bcor --merge ${threads}; \
          ldstore --bcor ${f}p.bcor --matrix ${f}p.ld; \
          sed -i -e "s/  */ /g; s/^ *//; /^$/d" ${f}p.ld'
      sed 's/\./p\./g' finemap.cfg > finemapp.cfg
      finemap --sss --in-files finemapp.cfg --n-causal-max 5 --corr-config 0.9
      echo "snpid region index snp_prob snp_log10bf rsid" > snpp.K20
      echo "rank config config_prob config_log10bf region" > configp.dat
      awk 'NR>1' st.bed | \
      parallel -j${threads} -C' ' '
          export f=chr{1}_{2}_{3}p; \
          awk '(NR>1 && NR<5){print $0, sub(/p.config/,"",FILENAME)}' >> configp.dat
          R --no-save < ${FM_location}/files/finemap-check.R > $f.chk
          export f=chr{1}_{2}_{3}; \
          cut -d" " -f10,11 $f.r > $f.tmp; \
          awk "(NR>1&&\$3>0.8&&\$4>1.3){print ENVIRON[\"f\"], \$0}" ${f}p.snp | \
          sort -k3,3 | \
          join -13 -22 - $f.tmp >> snpp.K20'
      echo "chr pos log10BF prob snpid rsid region" > snpp.dat
      awk '(NR>1){snpid=$1;gsub(/:|_/," ",$1);split($1,a," ");print a[1],a[2],$5,$4,snpid,$6,$2}' snpp.K20 | \
      sort -k1,1n -k2,2n >> snpp.dat
      export f="snpp"
      export p="PPAp.pdf"
      R --no-save < ${FM_location}/files/finemap-plot.R > finemapp-plot.log
   else
      awk 'NR>1' st.bed | \
      parallel -j${threads} -C' ' '
          export f=chr{1}_{2}_{3}; \
          ldstore --bcor $f.bcor --bplink $f --n-threads ${threads}; \  
          ldstore --bcor $f.bcor --merge ${threads}; \
          ldstore --bcor $f.bcor --matrix $f.ld --incl_variants $f.incl_variants; \
          sed -i -e "s/  */ /g; s/^ *//; /^$/d" $f.ld'
      echo "z;ld;snp;config;log;n-ind" > finemap.cfg
      awk 'NR>1' st.bed | \
      parallel -j${threads} -C ' ' '
          export f=chr{1}_{2}_{3}; \
          sort -k9,9g $f.r | \
          tail -n1|cut -d" " -f9 | \
          awk -vf=$f "{print sprintf(\"%s.z;%s.ld;%s.snp;%s.config;%s.log;%d\",f,f,f,f,f,int(\$1))}" >> finemap.cfg'
      finemap --sss --in-files finemap.cfg --n-causal-max 5 --corr-config 0.9
      echo "snpid region index snp_prob snp_log10bf rsid" > snp.K20
      echo "rank config config_prob config_log10bf region" > config.dat
      awk 'NR>1' st.bed | \
      parallel -j${threads} -C' ' '
          export f=chr{1}_{2}_{3}; \
          awk '(NR>1 && NR<5){print $0, sub(/.config/,"",FILENAME)}' >> config.dat
          R --no-save < ${FM_location}/files/finemap-check.R > $f.chk
          cut -d" " -f10,11 $f.r > $f.tmp; \
          awk "(NR>1&&\$3>0.8&&\$4>1.3){print ENVIRON[\"f\"], \$0}" $f.snp | \
          sort -k3,3 | \
          join -13 -22 - $f.tmp >> snp.K20'
      echo "chr pos log10BF prob snpid rsid region" > snp.dat
      awk '(NR>1){snpid=$1;gsub(/:|_/," ",$1);split($1,a," ");print a[1],a[2],$5,$4,snpid,$6,$2}' snp.K20 | \
      sort -k1,1n -k2,2n >> snp.dat
      export f="snp"
      export p="PPA.pdf"
      R --no-save < ${FM_location}/files/finemap-plot.R > finemap-plot.log
      if [ $LocusZoom -eq 1 ]; then
         awk 'NR>1' st.bed | \
         parallel -j${threads} -C' ' '
             export f=chr{1}_{2}_{3}; \
             awk "{if(NR==1) \$0=\$0 \" order\"; else \$0=\$0 \" \" NR-1;print}" $f.snp > $f.sav; \
             awk "NR==1" $f.sav | \
             awk "{print \$0 \" rsid\"}" > $f.snp; \
             awk "(NR>1)" $f.sav | \
             sort -k2,2 | \
             join -j2 - $f.rsid | \
             sort -k5,5n | \
             awk "{t=\$1;\$1=\$2;\$2=t};1" >> $f.snp'
         R --no-save < ${FM_location}/files/finemap-xlsx.R > finemap-xlsx.log
      fi
   fi
fi
