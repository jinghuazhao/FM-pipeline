// extract GEN file for each LD region in 1000G data

set more off

gzuse FUSION.dta.gz, clear
!rm -f extract.sh
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
      local f="chr`k'_`lowr'_`uppr'"
      outsheet snpid pos exp_freq_a1 info type rsid if pos>=`lowr' & pos<=`uppr' using `f'.info, names noquote replace nolab delim(" ")
      !echo -e "qctool -g chr`k'.gen.gz -og `f'.gen -incl-range `lowr'-`uppr' -omit-chromosome -sort;gzip -f `f'.gen" >> extract.sh
   }
   restore
}
cd `T'
!chmod u+x extract.sh
!./extract.sh
