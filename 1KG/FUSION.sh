# 1-3-2019 JHZ

wget -qO- https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2 | tar xfj - --strip-components=1
seq 22|parallel -j4 -C' ' 'plink-1.9 --bfile 1000G.EUR.{} --recode oxford gen-gz --out chr{}'
seq 22|awk -vp=1000G.EUR. '{print p $1}' > merge-list
plink-1.9 --merge-list merge-list --make-bed --out FUSION
plink-1.9 --bfile FUSION --freq --out FUSION
(
  echo "ID_1 ID_2 missing sex phenotype"
  echo "0 0 0 D B"
  awk '{print $1,$2,$3,$4,"NA"}' FUSION.fam
) > FUSION.sample
awk -vOFS="\t" '(NR>1){print $2,$5}' FUSION.frq > FUSION.dat

stata -b do FUSION.do
