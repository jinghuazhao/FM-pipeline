// SNPinfo.dta.gz

insheet rsid FreqA2 using lz-1.4.dat, case
sort rsid
gzsave lz-1.4, replace
insheet chr rsid m pos A1 A2 using lz-1.4.bim, case clear
gen RSnum=rsid
gen info=1
gen type=2
sort rsid
gzmerge using lz-1.4
gen snpid=string(chr)+":"+string(pos,"%12.0f")+"_"+cond(A1<A2,A1,A2)+"_"+cond(A1<A2,A2,A1)
sort chr pos
drop _merge
gen maf=cond(FreqA2<=0.5, FreqA2, 1-FreqA2)
sort chr pos snpid
gen MAC=2*503*maf
rename FreqA2 exp_freq_a1
order snpid pos exp_freq_a1 info type
gzsave "lz-1.4.dta.gz", replace

// .info

global home: env HOME
set more off

gzuse "lz-1.4.dta.gz", clear
tempfile f0
forval k=1/22 {
   preserve
   keep if chr==`k'
   save `f0', replace
   import delimited using $home/FM-pipeline/1KG/EUR.bed, asdouble clear
   keep chrom start end
   drop if end<=start
   gen chr=substr(chrom,4,2)
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
      outsheet snpid pos exp_freq_a1 info type rsid if pos>=`lowr' & pos<=`uppr' using `f'.info, names noquote replace nolab delim(" ")
   }
   restore
}
