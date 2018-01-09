/* extract GEN file for each LD region in EPIC-Omics HRC imputed data */
set more off
tempfile f0

local DIRGEN /gen_omics/data/EPIC-Norfolk/HRC
local DIRBGEN /scratch/tempjhz22/LDcalc/LD_MAGIC

gzuse `DIRGEN'/SNPinfo-chr1-22.dta.gz, clear
drop if chr==25
//gen maf=cond(FreqA2<=0.5, FreqA2, 1-FreqA2)
sort chr pos rsid
gen MAC=2*21044*maf

rename probid snp_id
rename rsid rs_id
rename pos position
rename FreqA2 exp_freq_a1
order snp_id rs_id position exp_freq_a1 info type
tostring chr, gen(CHR)
replace snp_id=CHR if snp_id!="---" & snp_id!=CHR

!rm -f `DIRBGEN'/Extract.sh

forval k=1/22 {
	preserve
	keep if chr==`k'
	save `f0', replace
	import delimited using doc/st.sav, varnames(nonames) asdouble delim(" ") clear
	keep chr start end
	drop if v3<=v2
//	replace v1="23" if v1=="X"
	destring v1, replace
	keep if v1==`k'
	drop v1
	sort v2 
	rename v2 St
	rename v3 En
	count
	local nclus=r(N)
	merge 1:1 _n using `f0', nogen
	forval j=1/`nclus' {
		local lowr=St[`j']
		local uppr=En[`j']
		outsheet rs_id if position>=`lowr' & position<=`uppr' & (MAC<3 | info<0.4) using `DIRBGEN'/exc`k'_`lowr'_`uppr'.txt, nonames noquote replace nolab
		outsheet snp_id rs_id position exp_freq_a1 info type RSnum if position>=`lowr' & position<=`uppr' & MAC>=3 & info>=0.4 using `DIRBGEN'/chr`k'_`lowr'_`uppr'.info, names noquote replace nolab delim(" ")
		!echo -e "sge \"/genetics/bin/qctool -g `DIRGEN'/chr`k'.gen.gz -og `DIRBGEN'/chr`k'_`lowr'_`uppr'.gen -incl-range `lowr'-`uppr' -omit-chromosome -excl-rsids exc`k'_`lowr'_`uppr'.txt -sort;gzip -f chr`k'_`lowr'_`uppr'.gen\"" >> `DIRBGEN'/Extract.sh
	}

	restore
}

cd `DIRBGEN'
!chmod u+x Extract.sh
!./Extract.sh
