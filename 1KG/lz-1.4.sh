#!/bin/bash
# 11-1-2019 JHZ

sbatch --wait $HOME/FM-pipeline/doc/lz-1.4.sb

for i in `seq 22`; do echo EUR1KG-$i; done > EUR.list
plink --merge-list EUR.list --make-bed --out EUR

qctool -filetype binary_ped -g EUR.bed -ofiletype gen -og EUR.gen.gz 

sbatch --wait $HOME/FM-pipeline/1KG/extract.sb

// generate .info files

plink-1.9 --bfile EUR --freq --out EUR
awk -vOFS="\t" '(NR>1){print $2,$5}' EUR.frq > EUR.dat

export TMPDIR=$HOME/FM-pipeline

stata -b do lz-1.4.do
