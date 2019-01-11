# 11-1-2019 JHZ

wget -qO- https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2 | tar xfj - --strip-components=1
seq 22|parallel -j4 -C' ' 'plink-1.9 --bfile 1000G.EUR.{} --recode oxford gen-gz --out chr{}'
seq 22|awk -vp=1000G.EUR. '{print p $1}' > merge-list
plink-1.9 --merge-list merge-list --make-bed --out EUR
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
  gzsave FUSION, replace
END
