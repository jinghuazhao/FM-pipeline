// extract GEN file for each LD region in EPIC-Omics HRC imputed data

set more off

local F /gen_omics/data/EPIC-Norfolk/HRC
local T /scratch/tempjhz22/LDcalc/HRC

gzuse if chr!=25 using `F'/SNPinfo.dta.gz, clear
sort chr pos rsid
gen maf=cond(FreqA2<=0.5, FreqA2, 1-FreqA2)
gen MAC=2*21044*maf
rename FreqA2 exp_freq_a1
order rsid pos exp_freq_a1 info type RSnum

tempfile f0
!rm -f `T'/Extract.sh
!touch `T'/Extract.sh
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
      local f="chr`k'_`lowr'_`uppr'"
      local qctool="/genetics/bin/qctool-1.4"
      local samples="/gen_omics/data/EPIC-Norfolk/HRC/EPIC-Norfolk.sample"
      local excl_samples="/genetics/data/gwas/1-11-17/exclude.dat"
      outsheet rsid pos exp_freq_a1 info type RSnum if pos>=`lowr' & pos<=`uppr' & MAC>=3 & info>=0.4 using `T'/`f'.txt, names noquote replace nolab delim(" ")
      !echo -e "sge \"`qctool' -g `F'/chr`k'.gen.gz -og `f'.gen -s `samples' -excl-samples `excl_samples' -incl-range `lowr'-`uppr' -omit-chromosome -snp-missing-call-rate 0.100001 -snp-missing-rate 0.05 -maf 0.01 1 -info 0.4 1 -sort; cut -d' ' -f1 `f'.gen > `f'.snpid; grep -w -f `f'.snpid `f'.txt > `f'.info; gzip -f `f'.gen\"" >> `T'/Extract.sh
   }
   restore
}

cd `T'
!chmod u+x Extract.sh
!./Extract.sh

// origina scripts
//    outsheet rsid if pos>=`lowr' & pos<=`uppr' & (MAC<3 | info<0.4) using `T'/exc`f'.txt, nonames noquote replace nolab
//    outsheet rsid pos exp_freq_a1 info type RSnum if pos>=`lowr' & pos<=`uppr' & MAC>=3 & info>=0.4 using `T'/chr`f'.info, names noquote replace nolab delim(" ")
//    !echo -e "sge \"/genetics/bin/qctool -g `F'/chr`k'.gen.gz -og chr`f'.gen -incl-range `lowr'-`uppr' -omit-chromosome -excl-rsids exc`f'.txt -sort; gzip -f chr`f'.gen\"" >> `T'/Extract.sh
