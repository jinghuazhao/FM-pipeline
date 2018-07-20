// extract GEN file for each LD region in EPIC-Omics HRC imputed data

set more off

local F /genetics/bin/FUSION/LDREF
local T /scratch/tempjhz22/LDcalc/1KG

gzuse `F'/SNPinfo.dta.gz, clear
gen maf=cond(FreqA2<=0.5, FreqA2, 1-FreqA2)
sort chr pos snpid
gen MAC=2*489*maf
rename FreqA2 exp_freq_a1
order snpid pos exp_freq_a1 info type

!rm -f `T'/Extract.sh
tempfile f0
forval k=1/22 {
   preserve
   keep if chr==`k'
   save `f0', replace
   import delimited using st.bed, asdouble delim(" ") clear
   keep chr start end
   drop if end<=start
   destring chr, replace
   keep if chr==`k'
   drop chr
   sort start
   count
   local nclus=r(N)
   merge 1:1 _n using `f0', nogen
   forval j=1/`nclus' {
      local lowr=start[`j']
      local uppr=end[`j']
      outsheet snpid if pos>=`lowr' & pos<=`uppr' & (MAC<3 | info<0.4) using `T'/exc`k'_`lowr'_`uppr'.txt, nonames noquote replace nolab
      outsheet snpid pos exp_freq_a1 info type rsid if pos>=`lowr' & pos<=`uppr' & MAC>=3 & info>=0.4 using `T'/chr`k'_`lowr'_`uppr'.info, names noquote replace nolab delim(" ")
      !echo -e "sge \"/genetics/bin/qctool -g `F'/chr`k'.gen.gz -og chr`k'_`lowr'_`uppr'.gen -incl-range `lowr'-`uppr' -omit-chromosome -excl-rsids exc`k'_`lowr'_`uppr'.txt -sort;gzip -f chr`k'_`lowr'_`uppr'.gen \"" >> `T'/Extract.sh
   }
   restore
}

cd `T'
!chmod u+x Extract.sh
!./Extract.sh
