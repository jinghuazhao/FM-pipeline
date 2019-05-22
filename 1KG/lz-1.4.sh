#!/bin/bash
# 3-5-2019 JHZ

sbatch --wait $HOME/FM-pipeline/doc/lz-1.4.sb

for i in `seq 22`; do echo EUR1KG-$i; done > lz-1.4.list
plink --merge-list lz-1.4.list --make-bed --out lz-1.4

qctool -filetype binary_ped -g lz-1.4.bed -ofiletype gen -og lz-1.4.gen.gz 

sbatch --wait $HOME/FM-pipeline/1KG/extract.sb

// generate .sample file

(
  echo "ID_1 ID_2 missing sex phenotype"
  echo "0 0 0 D B"
  awk '{print $1,$2,$3,$4,"NA"}' lz-1.4.fam
) > lz-1.4.sample

// generate .info files

plink-1.9 --bfile LocusZoom --freq --out lz-1.4
awk -vOFS="\t" '(NR>1){print $2,$5}' lz-1.4.frq > lz-1.4.dat

export TMPDIR=$HOME/FM-pipeline

stata -b do lz-1.4.do
