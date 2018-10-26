#!/bin/bash
# 26-10-2018 MRC-Epid JHZ

source gcta-slct.ini

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

if [ $region -eq 1 ]; then
   export flanking=250
   awk 'NR>1' st.bed | parallel -j1 --env flanking --env threads -C' ' '
       export f=chr{1}_{2}_{3}; \
       gcta64 --bfile $bfile $OPT1 $OPT2 --cojo-file $1.dat --cojo-slct --maf 0.000072 --extract-region-bp {1} {4} $flanking --thread-num $threads --out $f'
fi

echo "3. Adding snpid and rsid"

awk 'NR==1{print $0,"\tsnpid","\trsid"}' $1.jma.cojo > $1.jma.out
sort -k2,2 id3.txt > id3.tmp
awk 'NR>1{print $0, NR}' $1.jma.cojo | sort -k2,2 | \
join -j2 - id3.tmp | sort -k15,15n | awk '{$15="";print}' | awk '{t=$1;$1=$2;$2=t;gsub(/ /,"\t",$0)};1' >> $1.jma.out

rm id3.tmp id3.txt

echo "Done!"
