#!/bin/bash
# 16-1-2018 MRC-Epid JHZ

setup() {
stata <<END
gzuse /gen_omics/data/EPIC-Norfolk/HRC/SNPinfo
gen A1A2=cond(A1<A2,"_"+A1+"_"+A2,"_"+A2+"_"+A1)
gen snpid=string(chr)+":"+string(pos,"%12.0f")+A1A2
sort snpid
gen maf=cond(FreqA2<=0.5, FreqA2, 1-FreqA2)
gen MAC=2*21044*maf
outsheet snpid if (MAC<3 | info<0.4) using exclude.snps, noname noquote replace
keep rsid RSnum snpid
outsheet using SNPinfo.txt, delim(" ") noname noquote replace
END
export GEN=/gen_omics/data/EPIC-Norfolk/HRC
export sample=/gen_omics/data/EPIC-Norfolk/HRC/EPIC-Norfolk.sample
cd /scratch/tempjhz22/LDcalc/HRC
seq 22 | parallel --env GEN --env sample -C' ' 'sge "/genetics/bin/plink2 --bgen $GEN/chr{}.bgen --sample $sample --chr {} --make-bed --out chr{}"'
rm merge-list
touch merge-list
for i in $(seq 22); do echo chr${i} >> merge-list; done
/genetics/bin/plink-1.9 --merge-list merge-list --make-bed --out HRC
}

echo "SNP A1 A2 freq b se p N" > gcta.dat
sort -k9,9n -k10,10n $1 | awk '
{
  a1=toupper($2)
  a2=toupper($3)
  chr=$9
  pos=$10
  if (a1<a2) snpid=chr ":" pos "_" a1 "_" a2
  else snpid=chr ":" pos "_" a2 "_" a1
  $1=snpid
  $2=a1
  $3=a2
  print $1,$2,$3,$4,$5,$6,$7,$8
}' | sort -k1,1 | join -13 -21 SNPinfo.txt - | awk '{$1=$2="";print}' | awk '{$1=$1};1' >> gcta.dat

gcta64 --bfile /scratch/tempjhz22/LDcalc/HRC/HRC --remove exclude.id --exclude exclude.snps --cojo-file gcta.dat --cojo-slct --thread-num 10 --out gcta

keep_note() {
  grep -w -f 97.snps repro.txt | sort -k9,9n -k10,10n | awk '{
    a1=toupper($2)
    a2=toupper($3)
    chr=$9
    pos=$10
    if (a1<a2) snpid=chr ":" pos "_" a1 "_" a2
    else snpid=chr ":" pos "_" a2 "_" a1
    print snpid
  }' > 97.snpid
  gcta64 --bfile $GEN/HRC --remove exclude.id --exclude exclude.snps --cojo-file gcta.dat --cojo-cond 97.snpid --thread-num 10 --out gcta
  export GEN=/scratch/tempjhz22/LDcalc/HRC
  gcta64 --bfile $GEN/HRC --remove exclude.id --exclude exclude.snps --cojo-file gcta.dat --cojo-top-SNPs 5 --thread-num 10 --out gcta.top
  export GENI=/gen_omics/data/EPIC-Norfolk/HRC
  export GENO=/scratch/tempjhz22/LDcalc/HRC
  export sample_file=$GENI/EPIC-Norfolk.sample
  export FM_pipeline=/genetics/bin/FM-pipeline
  parallel --env GENI --env GENO --env sample_file --env FM_pipeline -C' ' '
    sge "export f=chr{1}_{2}; \
    gunzip -c $GENI/\$f.gen.gz | \
    awk -f $FM_pipeline/files/order.awk chr={1} | \
    gzip -f > $GENO/\$f.gen.gz; \
    /genetics/bin/gtool -G --g $GENO/\$f.gen.gz --s ${sample_file} --ped $GENO/\$f.ped --map $GENO/\$f.map \
          --missing 0.05 --threshold 0.9 --log \$f.log --snp --alleles --chr {1}"' ::: $(seq 22) ::: $(seq 30)
  parallel -C' ' '
  sge "cd /scratch/tempjhz22/LDcalc/HRC;\
  zcat chr{}_1.ped.gz chr{}_2.ped.gz chr{}_3.ped.gz chr{}_4.ped.gz chr{}_5.ped.gz chr{}_6.ped.gz chr{}_7.ped.gz chr{}_8.ped.gz chr{}_9.ped.gz chr{}_10.ped.gz chr{}_11.ped.gz chr{}_12.ped.gz chr{}_13.ped.gz chr{}_14.ped.gz chr{}_15.ped.gz chr{}_16.ped.gz chr{}_17.ped.gz chr{}_18.ped.gz chr{}_19.ped.gz chr{}_20.ped.gz chr{}_21.ped.gz chr{}_22.ped.gz chr{}_23.ped.gz chr{}_24.ped.gz chr{}_25.ped.gz chr{}_26.ped.gz chr{}_27.ped.gz chr{}_28.ped.gz chr{}_29.ped.gz chr{}_30.ped.gz > HRC{}.ped;\
  zcat chr{}_1.map.gz chr{}_2.map.gz chr{}_3.map.gz chr{}_4.map.gz chr{}_5.map.gz chr{}_6.map.gz chr{}_7.map.gz chr{}_8.map.gz chr{}_9.map.gz chr{}_10.map.gz chr{}_11.map.gz chr{}_12.map.gz chr{}_13.map.gz chr{}_14.map.gz chr{}_15.map.gz chr{}_16.map.gz chr{}_17.map.gz chr{}_18.map.gz chr{}_19.map.gz chr{}_20.map.gz chr{}_21.map.gz chr{}_22.map.gz chr{}_23.map.gz chr{}_24.map.gz chr{}_25.map.gz chr{}_26.map.gz chr{}_27.map.gz chr{}_28.map.gz chr{}_29.map.gz chr{}_30.map.gz > HRC{}.map"' ::: $(seq 22)
}
