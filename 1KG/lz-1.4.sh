#!/bin/bash
# 10-1-2019 JHZ

sbatch --wait lz-1.4.sb

for i in `seq 22`; do echo EUR1KG-$i; done > EUR.list
plink --merge-list EUR.list --make-bed --out EUR

qctool -filetype binary_ped -g EUR.bed -ofiletype gen -og EUR.gen.gz 

sbatch --wait extract.sb
