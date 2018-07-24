#!/bin/bash
# 24-7-2018 MRC-Epid JHZ

## SETTINGS

export PATH=/genetics/bin:/usr/local/bin:$PATH:/genetics/data/software/bin
export R_LIBS=/genetics/bin/R:/usr/local/lib64/R/library:/genetics/data/software/R
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib64/R/lib:/genetics/data/software/lib

## OPTIONS

# nonempty value to skip parallel sessions for data handling and go directly to analysis

export dry_run=

# software for analysis; set flags to 1 to enable

export clumping=0
export CAVIAR=0
export CAVIARBF=0
export FM_summary=0
export GCTA=0
export JAM=1
export LocusZoom=0
export fgwas=0
export finemap=1
export fgwas_location_1kg=/genetics/data/software/fgwas/1000-genomes
export FM_location=/genetics/bin/FM-pipeline

export wd=$(pwd)
# GEN files named chr{chr}_{start}_{end}.gen.gz
export GEN_location=$FM_location/1KG/LD-blocks
# sample file
export sample_file=$FM_location/1KG/EUR.sample
# wholegenome genotype file
export HRC=/gen_omics/data/EPIC-Norfolk/HRC/binary_ped
export bfile=$HRC/HRC
export remove_sample=$HRC/exclude.id
export exclude_snp=$HRC/exclude.snps
# number of threads
export threads=1
export LD_MAGIC=0
export LD_PLINK=0

## ANALYSIS

if [ $# -lt 1 ] || [ "$args" == "-h" ]; then
    echo "Usage: fmp.sh <input>"
    echo "where <input> is in sumstats format:"
    echo "SNP A1 A2 freqA1 beta se P N chr pos"
    echo "where SNP is RSid, A1 is effect allele"
    exit
fi

export args=$1
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
echo "--> $rt.input, st.bed"
awk '{
  $2=toupper($2)
  $3=toupper($3)
};1' $args > $rt.input
ln -sf $wd/st.bed

export OPTs=""
if [ ! -z "$dry_run" ]; then OPTs="--dry-run"; fi

echo "--> binary_ped"
awk 'NR>1' st.bed | parallel $OPTs -j${threads} --env sample_file --env FM_location --env GEN_location -C' ' '
    export f=chr{1}_{2}_{3}; \
    gunzip -c $GEN_location/$f.gen.gz | \
    awk -f $FM_location/files/order.awk chr={1} > $GEN_location/$f.ord;\
    qctool -filetype gen -g $GEN_location/$f.ord -s ${sample_file} -ofiletype binary_ped -og $GEN_location/$f \
          -threads $threads -threshold 0.9 -log $f.log -assume-chromosome {1}'
echo "region-specific data"
awk 'NR>1' st.bed | parallel $OPTs -j${threads} -C' ' '
    export f=chr{1}_{2}_{3}; \
    awk "(\$9==chr && \$10 >= l && \$10 <= u){if(\$2<\$3) {a1=\$2; a2=\$3;} else {a1=\$3; a2=\$2};\
         \$0=\$0 \" \" \$9 \":\" \$10 \"_\" a1 \"_\" a2;print}" chr={1} l={2} u={3} $rt.input | \
         sort -k11,11 > $f.txt'
echo "--> GWAS auxiliary files"
awk 'NR>1' st.bed | parallel $OPTs -j${threads} --env GEN_location -C' ' '
    export f=chr{1}_{2}_{3}; \
    sort -k2,2 $GEN_location/$f.bim | \
    join -111 -22 $f.txt - | \
    sort -k11,11 > $f.incl; \
    awk "{print \$10,\$11,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$2,\$1,\$6/\$7}" $f.incl > $f.r; \
    cut -d" " -f10,11 $f.r > $f.rsid'
awk 'NR>1' st.bed | parallel $OPTs -j${threads} --env wd -C' ' '
    export f=chr{1}_{2}_{3}; \
    cut -d" " -f11,12 $f.r > $f.z; \
    awk "{print \$1}" $f.incl > $f.inc; \
    awk "{print \$1,\$4,\$3,\$15,\$16}" $f.incl > $f.a; \
    echo "RSID position chromosome A_allele B_allele" > $f.incl_variants; \
    awk "{print \$1,\$11,\$10,\$4,\$3}" $f.incl >> $f.incl_variants'
awk 'NR>1' st.bed | parallel $OPTs -j${threads} -C' ' '
    export f=chr{1}_{2}_{3}; \
         grep -f $f.inc $f.txt | \
         sort -k11,11 > $f.dat'

echo "--> bfile"
awk 'NR>1' st.bed | parallel $OPTs -j${threads} --env GEN_location -C' ' '
    export f=chr{1}_{2}_{3}; \
    plink-1.9 --bfile $GEN_location/$f --extract $f.inc \
    --make-bed --keep-allele-order --a2-allele $f.a 3 1 --out $f'

if [ $LD_MAGIC -eq 1 ]; then
    awk 'NR>1' st.bed | parallel -j${threads} --env threads --env sample_file --env FM_location --env GEN_location -C' ' '
    export f=chr{1}_{2}_{3}; \
    gunzip -c $GEN_location/$f.gen.gz | \
    awk -f $FM_location/files/order.awk chr={1} > $GEN_location/$f.ord;\
    qctool_v2.0 -filetype gen -g $GEN_location/$f.ord -s ${sample_file} -ofiletype gen -og $GEN_location/$f.magic.gen \
          -threads $threads -threshhold 0.9 -log $f.log -omit-chromosome;\
    awk -f $FM_location/files/info.awk c=2 $GEN_location/$f.info > $GEN_location/$f.magic.info; \
    gzip -f $GEN_location/$f.magic.gen; \
    Rscript --vanilla $FM_location/files/computeCorrelationsImpute2forFINEMAP.r \
            $GEN_location/$f.magic.info $GEN_location/$f.magic.gen.gz {1} {2} {3} 0.01 0.4 $f.magic $threads; \
    Rscript --vanilla $FM_location/files/lowtri2square.r'
fi

if [ $LD_PLINK -eq 1 ]; then
   awk 'NR>1' st.bed | parallel -j${threads} --env threads -C' ' '
       export f=chr{1}_{2}_{3}; \
       plink-1.9 --bfile $f --maf 0.0001 --freq --threads 3 --out $f; \
       awk "(\$5<0.0001){print \$2}" $f.frq > $f.excl; \
       cp $f.z $f.sav; \
     # grep -w -v -f $f.excl $f.sav > $f.z; \
       plink-1.9 --bfile $f --maf 0.0001 --r square --threads 3 --out $f; \
       sed "s/\t/ /g" $f.ld > $f.plink'
     # grep -w -v -f $f.excl $f.r below
fi

if [ $clumping -eq 1 ]; then
   echo "--> clumping"
   awk '
   {
      if (NR==1) print "snpid", "P"
      chr=$9
      pos=$10
      a1=$2
      a2=$3
      if (a1>a2) {
         snpid=chr ":" pos "_" a2 "_" a1
      } else {
         snpid=chr ":" pos "_" a1 "_" a2
      }
      print snpid, $7
   }' OFS='\t' $rt.input > $rt.tab
   plink-1.9 --bfile $bfile --remove $remove_sample --exclude $exclude_snp --clump $rt.tab \
             --clump-field P \
             --clump-kb 500 \
             --clump-p1 5e-08 \
             --clump-r2 0.1 \
             --clump-snp-field snpid \
             --out $rt.clump
   awk 'NR>1' st.bed | parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}; \
       awk "{if (NR==1) print \"snpid\", \"P\"; print \$11,\$7}" OFS="\t" $f.dat > $f.tab; \
       plink-1.9 --bfile $f --clump $f.tab \
            --clump-field P \
            --clump-kb 500 \
            --clump-p1 5e-08 \
            --clump-r2 0.1 \
            --clump-snp-field snpid \
            --out $f.clump'
fi

if [ $CAVIAR -eq 1 ] || [ $CAVIARBF -eq 1 ] || [ $finemap -eq 1 ]; then
   awk 'NR>1' st.bed | parallel -j${threads} --env threads -C' ' '
       export f=chr{1}_{2}_{3}; \
       ldstore --bcor $f.bcor --bplink $f --n-threads ${threads}; \
       ldstore --bcor $f.bcor --merge ${threads}; \
       ldstore --bcor $f.bcor --matrix $f.ld --incl_variants $f.incl_variants; \
       sed -i -e "s/  */ /g; s/^ *//; /^$/d" $f.ld'
fi

# CAusal Variants Identication in Associated Regions (CAVIAR)

if [ $CAVIAR -eq 1 ]; then
   echo "--> CAVIAR"
   awk 'NR>1' st.bed | parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}; \
       CAVIAR -z $f.z -l $f.ld -r 0.9 -o $f'
fi

# CAusal Variants Identication in Associated Regions BF (CAVIARBF)

if [ $CAVIARBF -eq 1 ]; then
   echo "--> CAVIARBF"
   awk 'NR>1' st.bed | parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}; \
       caviarbf -z $f.z -r $f.ld -n $(sort -k9,9g $f.r | \
       tail -n1 | cut -d" " -f9) -t 0 -a 0.1 -c 3 --appr -o $f.caviarbf'
fi

# Fine-mapping method only using summary statistics (FM-summary)

if [ $FM_summary -eq 1 ]; then
   echo "--> FM-summary"
   echo "region chr pos A B Freq1 Effect StdErr P N SNP inCredible probNorm cumSum" | \
   sed 's/ /\t/g' > FM-summary.txt
   awk 'NR>1' st.bed | parallel -j${threads} --env FM_location -C' ' '
       export f=chr{1}_{2}_{3}; \
       $FM_location/files/getCredible.r; \
       awk "!(/SNP/&&/inCredible/){print f, \$0}" OFS="\t" f=$f $f.cre >> FM-summary.txt'
fi

# Genome-wide Complex Trait Analysis (GCTA)

if [ $GCTA -eq 1 ]; then
   echo "--> GCTA"
   awk 'NR>1' st.bed | parallel -j${threads} --env FM_location --env GEN_location -C' ' '
       export f=chr{1}_{2}_{3}; \
       awk -f $FM_location/files/info.awk c=1 chr={1} $GEN_location/$f.info | \
       sort -k2,2 > $f.tmp; \
       sort -k2,2 $GEN_location/$f.bim | \
       join -j2 $f.tmp - | \
       awk -vOFS="\t" "{print \$7,\$6,0,\$2,\$10,\$11,\$9}" > ${f}_map; \
       sort -k4,4 ${f}_map | \
       join -111 -24 $f.r - | \
       grep -f $f.inc | \
       awk -f $FM_location/files/ma.awk > $f.ma; \
       gcta64 --bfile $f --cojo-file $f.ma --cojo-joint --cojo-collinear 0.9 --out $f; \
       gcta64 --bfile $f --cojo-file $f.ma --cojo-slct --maf 0.000072 --out $f; \
       grep {5} $f.r | \
       cut -d" " -f11 > $f.snpid; \
       gcta64 --bfile $f --cojo-file $f.ma --cojo-cond $f.snpid --out $f; \
       gcta64 --bfile $f --cojo-file $f.ma --cojo-top-SNPs 1 --out $f.top; \
       cut -d" " -f10,11 $f.r | \
       sort -k2,2 | \
       sed "s/ /\t/g">$f.tmp'
# --cojo-slct <==> jma.cojo, ldr.cojo
   echo "region SNP Chr bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r rsid" > gcta-slct.csv
   ls *.jma.cojo|sed 's/\.jma\.cojo//g' | parallel -j1 -C' ' '
       echo "SNP Chr bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r rsid" > {}.jma; \
       sort -k2,2 {}.jma.cojo | \
       join -j2 - {}.tmp >> {}.jma'
   awk 'NR>1' st.bed | parallel -j1 -C' ' '
       export f=chr{1}_{2}_{3}; \
       if [ -f $f.jma ]; then awk "!/SNP/{print f, \$0}" f=$f $f.jma >> gcta-slct.csv; fi'
   sed -i 's/ /,/g' gcta-slct.csv
# --cojo-cond <==> given.cojo, cma.cojo
   echo "region SNP Chr bp refA freq b se p n freq_geno bC bC_se pC rsid" > gcta-cond.csv
   ls *cma.cojo|sed 's/\.cma\.cojo//g' | parallel -j1 -C' ' '
       echo "SNP Chr bp refA freq b se p n freq_geno bC bC_se pC rsid" > {}.cma; \
       sort -k2,2 {}.cma.cojo | \
       join -j2 - {}.tmp >> {}.cma'
       awk 'NR>1' st.bed | parallel -j1 -C' ' '
       export f=chr{1}_{2}_{3}; \
       if [ -f $f.cma ]; then awk "!/SNP/{print f, \$0}" f=$f $f.cma >> gcta-cond.csv; fi'
   sed -i 's/ /,/g' gcta-cond.csv
# --cojo-top-SNPs <==> top.jma.cojo, top.ldr.cojo
   echo "region SNP Chr bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r rsid" > gcta-top.csv
   ls *top.jma.cojo | \
   sed 's/\.top\.jma\.cojo//g' | parallel -j1 -C' ' '
       echo "SNP Chr bp refA freq b se p n freq_geno bJ bJ_se pJ LD_r rsid" > {}.top.jma; \
       sort -k2,2 {}.top.jma.cojo | \
       join -j2 - {}.tmp >> {}.top.jma'
   awk 'NR>1' st.bed | parallel -j1 -C' ' '
       export f=chr{1}_{2}_{3}; \
       if [ -f $f.top.jma ]; then awk "!/SNP/{print f, \$0}" f=$f $f.top.jma >> gcta-top.csv; fi'
   sed -i 's/ /,/g' gcta-top.csv
fi

#  Joint Analysis of Marginal SNP Effects (JAM)

if [ $JAM -eq 1 ]; then
   echo "--> JAM"
   awk 'NR>1' st.bed | parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}; \
       grep {5} $f.r | \
       cut -d" " -f11 > $f.snpid; \
       plink-1.9 --bfile $f --exclude $f.snpid --indep-pairwise 500kb 1 0.80 --maf 0.0001 --out $f; \
       cat $f.snpid >> $f.prune.in; \
       grep -w -f $f.prune.in $f.a > $f.p; \
       grep -w -f $f.prune.in $f.dat > ${f}p.dat; \
       plink-1.9 --bfile $f --extract $f.prune.in --keep-allele-order --a2-allele $f.p 3 1 --make-bed --out ${f}p'
   awk 'NR>1' st.bed | parallel -j${threads} --env FM_location -C' ' '
       export f=chr{1}_{2}_{3}p; \
       R -q --no-save < ${FM_location}/files/JAM.R > $f.log'
   rm -f jam.top jam.txt
   touch jam.top jam.txt
   awk 'NR>1' st.bed | parallel -j1 -C' ' 'export f=chr{1}_{2}_{3};awk "NR==2&&\$2>0 {print f}" f=$f ${f}p.sum' >> jam.top
   cat jam.top | parallel -j1 -C' ' 'echo -e "\n" {} >> jam.txt;cat {}p.top {}p.jam >> jam.txt'
   echo "Region Chr SNP rsid PostProb_model PostProb Median CrI_Lower CrI_Upper Median_Present CrI_Lower_Present CrI_Upper_Present BF" > jam.out
   ls *sel | parallel -j1 -C' ' '
       awk "NR>1{sub(/\p.sel/,\"\",FILENAME);split(\$2,a,\":\");\$1=FILENAME \" \" a[1];print}"' | \
       sort -k1,1n >> jam.out
   R -q --no-save < ${FM_location}/files/JAM-cs.R > JAM-cs.log
fi

# LocusZoom

if [ $LocusZoom -eq 1 ]; then
   echo "--> LocusZoom"
   awk 'NR>1' st.bed | parallel -j${threads} -C' ' '
       export f=chr{1}_{2}_{3}; \
       awk "{OFS=\"\\t\";if(NR==1) print \"MarkerName\",\"P-value\",\"Weight\"; print \$10,\$8,\$9}" $f.r > $f.lz'
   awk 'NR>1' st.bed | parallel -j1 -C' ' '
       rm -f ld_cache.db; \
       locuszoom-1.4 --source 1000G_March2012 --build hg19 --pop EUR --metal chr{1}_{2}_{3}.lz --plotonly --chr {1} --start {2} --end {3} --no-date; \
       pdftopng chr{1}_{2}-{3}.pdf -r 300 {5}'
   R -q --no-save < ${FM_location}/files/lz.R > lz.log
fi

# functional genomic information with a genome-wide association study (fGWAS)

if [ $fgwas -eq 1 ]; then
   echo "--> fgwas"
   # obtain annotations
   seq 22 | parallel -j${threads} --env fgwas_location_1kg -C' ' '
       if [ ! -f $fgwas_location_1kg/chr{}.gene ]; then
          gunzip -c $fgwas_location_1kg/chr{}.annot.wdist.wcoding.gz | \
          awk "(NR>1){print \$1,\$2,\$3,\$(NF-6),\$(NF-5),\$(NF-2),\$(NF-1),\$NF}" | \
          sort -k2,2 > $fgwas_location_1kg/chr{}.gene \
       fi'
   # specify regions
   # -/+ flanking position
   export flanking=250000
   awk -vfl=${flanking} 'NR>1{l=$2;u=$3;print $5,$1,$4,l,u,NR}' st.bed | \
   sort -k1,1 > fgwas.snplist
   cat fgwas.snplist | parallel -j${threads} --env FM_location --env fgwas_location_1kg -C' ' '
       export f=chr{2}_{4}_{5}; \
       awk -vsn={6} -f $FM_location/files/fgwas.awk $f.r | \
       sort -k3,3 | \
       join -13 -22 - $fgwas_location_1kg/chr{2}.gene | \
       awk -f $FM_location/files/gene.awk | \
       gzip -fc > $f.fgwas.gz'
   # tally for -fine option
   rm -f fgwas.tmp
   touch fgwas.tmp
   sort -k6,6n fgwas.snplist | parallel -j1 -C' ' '
     gunzip -c chr{2}_{4}_{5}.fgwas.gz | \
     awk "(NR>1 && !/INDEL/)" >> fgwas.tmp
   '
   echo "SNPID CHR POS Z F N ens_coding_exons ens_noncoding_exons tssdist syn nonsyn SEGNUMBER" > fgwas.fine
   sort -k12,12n -k2,2 -k3,3n fgwas.tmp >> fgwas.fine
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

# finemap

if [ $finemap -eq 1 ]; then
   echo "--> finemap"
   echo "z;ld;snp;config;log;n_samples" > finemap.cfg
   awk 'NR>1' st.bed | parallel -j${threads} -C ' ' '
       export f=chr{1}_{2}_{3}; \
       sort -k9,9g $f.r | \
       tail -n1 | \
       cut -d" " -f9 | \
       awk -vf=$f "{print sprintf(\"%s.z;%s.ld;%s.snp;%s.config;%s.log;%d\",f,f,f,f,f,int(\$1))}" >> finemap.cfg'
   finemap --sss --in-files finemap.cfg --n-causal-snps 5 --corr-config 0.9
   awk 'NR>1' st.bed | parallel -j1 --env FM_location -C' ' '
       export f=chr{1}_{2}_{3}; \
       awk "{if(NR==1) \$0=\$0 \" order\"; else \$0=\$0 \" \" NR-1;print}" $f.snp > $f.sav; \
       awk "NR==1" $f.sav | \
       awk "{print \$0 \" rsid\"}" > $f.snp; \
       awk "(NR>1)" $f.sav | \
       sort -k2,2 | \
       join -j2 - $f.rsid | \
       sort -k5,5n | \
       awk "{t=\$1;\$1=\$2;\$2=t};1" >> $f.snp; \
       R -q --no-save < ${FM_location}/files/finemap.R > $f.out'
fi

# cherry-picking of information on z, ld, or possibly other things for signals with available data

if [ $GCTA -eq 1 ] && [ $JAM -eq 1 ] && [ $finemap -eq 1 ]; then
   R -q --no-save < ${FM_location}/files/gcta-jam-finemap.R > gcta-jam-finemap.log
   sort -k1,1 ld > ld.1
   # As before the same setup for gcta-slct.sh is used here but may complain otherwise
   gunzip -c /gen_omics/data/EPIC-Norfolk/HRC/binary_ped/id3.txt.gz | \
   awk '{snpid=$2;gsub(/:|\_/," ",snpid);split(snpid,a);Chr=a[1];Pos=a[2];print $0,Chr,Pos}' | \
   sort -k1,1 | \
   join -j1 - ld.1 | \
   sort -k4,4n -k5,5n | \
   awk '{print $2}' > ld.dat
   rm ld.1
   plink-1.9 --bfile $bfile --remove $remove_sample --exclude $exclude_snp \
             --extract ld.dat --r2 triangle spaces --out ld
   awk 'NR>1{if(NR==2) printf $8;else printf " " $8}' id | \
   awk '1' > ld.mat
   sed -e 's/[[:space:]]*$//' ld.ld >> ld.mat
   paste -d' ' id ld.mat > ld.txt
   R -q --no-save < ${FM_location}/files/ld.R > ld.log
fi

# However, we could opt to recalculate LD from the lists formed from GCTA and JAM results.

# obsolete with gtool/plink-1.9 handling gen/ped
#   gtool -G --g $GEN_location/$f.ord --s ${sample_file} --ped $GEN_location/$f.ped --map $GEN_location/$f.map \
#         --missing 0.05 --threshold 0.9 --log $f.log --snp --alleles --chr {1}'
#   plink-1.9 --file $GEN_location/$f --missing-genotype N --extract $f.inc \
#         --make-bed --keep-allele-order --a2-allele $f.a 3 1 --out $f'
