#!/bin/bash
# 2-2-2018 MRC-Epid JHZ

export PATH=/genetics/bin:/usr/local/bin:$PATH
export rt=/gen_omics/data/EPIC-Norfolk/HRC/binary_ped
export bfile=$rt/HRC
export idfile=$rt/id3.txt.gz
export remove_sample=$rt/exclude.id
export exclude_snp=$rt/exclude.snps
export threads=10

echo "1. Data preparation"
gunzip -c $idfile | sort -k1,1 > id3.txt
echo "SNP A1 A2 freq b se p N" > $1.dat
sort -k9,9n -k10,10n $1 | \
awk '!(/SNP/&&/A1/&&/A2/&&/freq/&&/b/&&/se/&&/p/&&/N/){
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
}' | sort -k1,1 | join -j1 id3.txt - | awk '{$1=$3="";print}' | awk '{$1=$1};1' >> $1.dat

echo "2. GCTA --cojo analysis"
export OPT1=""
if [ -f $remove_sample ] && [ ! -z "$remove_sample" ]; then export OPT1="--remove $remove_sample"; fi
export OPT2=""
if [ -f $exclude_snp ] && [ ! -z "$exclude_snp" ]; then export OPT2="--exclude $exclude_snp"; fi

gcta64 --bfile $bfile $OPT1 $OPT2 --cojo-file $1.dat --cojo-slct --maf 0.000072 --thread-num $threads --out $1

echo "3. Adding snpid and rsid"

awk 'NR==1{print $0,"\tsnpid","\trsid"}' $1.jma.cojo > $1.jma.out
sort -k2,2 id3.txt > id3.tmp
awk 'NR>1{print $0, NR}' $1.jma.cojo | sort -k2,2 | \
join -j2 - id3.tmp | sort -k15,15n | awk '{$15="";print}' | awk '{t=$1;$1=$2;$2=t;gsub(/ /,"\t",$0)};1' >> $1.jma.out

rm id3.tmp id3.txt

echo "Done!"

how_to_setup() {
stata <<END
gzuse /gen_omics/data/EPIC-Norfolk/HRC/SNPinfo
gen snpid=string(chr)+":"+string(pos,"%12.0f")+cond(A1<A2,"_"+A1+"_"+A2,"_"+A2+"_"+A1)
sort snpid
gen maf=cond(FreqA2<=0.5, FreqA2, 1-FreqA2)
gen MAC=2*21044*maf
outsheet snpid if (MAC<3 | info<0.4) using exclude.snps, noname noquote replace
keep snpid rsid RSnum
order snpid rsid RSnum
outsheet using /gen_omics/data/EPIC-Norfolk/HRC/binary_ped/id3.txt, delim(" ") noname noquote replace
!gzip -f /gen_omics/data/EPIC-Norfolk/HRC/binary_ped/id3.txt
END
export GEN=/gen_omics/data/EPIC-Norfolk/HRC

## The following is now updated as in files/st.sh
export sample=/gen_omics/data/EPIC-Norfolk/HRC/EPIC-Norfolk.sample
cd /gen_omics/data/EPIC-Norfolk/HRC/binary_ped
seq 22 | parallel --env GEN --env sample -C' ' 'sge "/genetics/bin/plink2 --bgen $GEN/chr{}.bgen --sample $sample --chr {} --make-bed --out chr{}"'
rm -f merge-list
touch merge-list
for i in $(seq 22); do echo chr${i} >> merge-list; done
/genetics/bin/plink-1.9 --merge-list merge-list --make-bed --out HRC
}
