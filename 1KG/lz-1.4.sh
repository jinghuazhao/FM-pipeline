#!/bin/bash
# 10-1-2019 JHZ

sbatch --wait $HOME/FM-pipeline/doc/lz-1.4.sb

for i in `seq 22`; do echo EUR1KG-$i; done > EUR.list
plink --merge-list EUR.list --make-bed --out EUR

qctool -filetype binary_ped -g EUR.bed -ofiletype gen -og EUR.gen.gz 

sbatch --wait $HOME/FM-pipeline/doc/extract.sb

// generate SNPinfo.dta

plink-1.9 --bfile EUR --freq --out EUR
awk -vOFS="\t" '(NR>1){print $2,$5}' EUR.frq > EUR.dat
stata <<END
  insheet rsid FreqA2 using EUR.dat, case
  sort rsid
  gzsave EUR, replace
  insheet chr rsid m pos A1 A2 using EUR.bim, case clear
  gen RSnum=rsid
  gen info=1
  gen type=2
  sort rsid
  gzmerge using EUR
  gen snpid=string(chr)+":"+string(pos,"%12.0f")+"_"+cond(A1<A2,A1,A2)+"_"+cond(A1<A2,A2,A1)
  sort chr pos
  drop _merge
  gen maf=cond(FreqA2<=0.5, FreqA2, 1-FreqA2)
  sort chr pos snpid
  gen MAC=2*489*maf
  rename FreqA2 exp_freq_a1
  order snpid pos exp_freq_a1 info type
  gzsave SNPinfo, replace
END

// extract GEN file for each LD region in 1000G data

set more off

gzuse SNPinfo.dta.gz, clear
tempfile f0
forval k=1/22 {
   preserve
   keep if chr==`k'
   save `f0', replace
   import delimited using EUR.bed, asdouble clear
   keep chr start end
   drop if end<=start
   gen chr=substr(chrom,4,2)
   keep if chr==`k'
   drop chr
   sort start
   count
   local nclus=r(N)
   merge 1:1 _n using `f0', nogen
   forval j=1/`nclus' {
      local lowr=start[`j']
      local uppr=end[`j']
      local f="chr`k'_`lowr'_`uppr'"
      outsheet snpid pos exp_freq_a1 info type rsid if pos>=`lowr' & pos<=`uppr' using LocusZoom-1.4/`f'.info, names noquote replace nolab delim(" ")
   }
   restore
}
