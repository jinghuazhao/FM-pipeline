// extract GEN file for whole chromosome in EPIC-Omics HRC imputed data

set more off

local F /gen_omics/data/EPIC-Norfolk/HRC
local T /scratch/tempjhz22/LDcalc/HRC

gzuse if chr!=25 using `F'/SNPinfo.dta.gz, clear
sort chr pos rsid
gen maf=cond(FreqA2<=0.5, FreqA2, 1-FreqA2)
gen MAC=2*21044*maf
rename FreqA2 exp_freq_a1
order rsid pos exp_freq_a1 info type RSnum

!rm -f `T'/Extract.sh

forval k=1/22 {
   preserve
   keep if chr==`k'
   outsheet rsid if (MAC<3 | info<0.4) using `T'/exc`k'.txt, nonames noquote replace nolab
   outsheet rsid pos exp_freq_a1 info type RSnum if MAC>=3 & info>=0.4 using `T'/chr`k'.info, names noquote replace nolab delim(" ")
   !echo -e "sge \"/genetics/bin/qctool -g `F'/chr`k'.gen.gz -og chr`k'.gen.gz -omit-chromosome -excl-rsids exc`k'.txt -sort" >> `T'/Extract.sh
   restore
}

cd `T'
!chmod u+x Extract.sh
!./Extract.sh
