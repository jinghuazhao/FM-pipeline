// extract GEN file for each LD region in EPIC-Omics HRC imputed data

set more off

local F /gen_omics/data/EPIC-Norfolk/HRC
local T /scratch/tempjhz22/LDcalc/HRC

gzuse if chr!=25 using `F'/SNPinfo.dta.gz, clear
gen maf=cond(FreqA2<=0.5, FreqA2, 1-FreqA2)
sort chr pos rsid
gen MAC=2*21044*maf
rename FreqA2 exp_freq_a1
order rsid pos exp_freq_a1 info type RSnum

tempfile f0
!rm -f `T'/Extract.sh
forval k=1/22 {
   preserve
   keep if chr==`k'
   save `f0', replace
   import delimited using st.bed, asdouble delim(" ") clear
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
      local f="`k'_`lowr'_`uppr'"
      outsheet rsid if pos>=`lowr' & pos<=`uppr' & (MAC<3 | info<0.4) using `T'/exc`f'.txt, nonames noquote replace nolab
      outsheet rsid pos exp_freq_a1 info type RSnum if pos>=`lowr' & pos<=`uppr' & MAC>=3 & info>=0.4 using `T'/chr`f'.info, names noquote replace nolab delim(" ")
      !echo -e "sge \"/genetics/bin/qctool -g `F'/chr`k'.gen.gz -og chr`f'.gen.gz -incl-range `lowr'-`uppr' -omit-chromosome -excl-rsids exc`f'.txt -sort\"" >> `T'/Extract.sh
   }
   restore
}

cd `T'
!chmod u+x Extract.sh
!./Extract.sh
