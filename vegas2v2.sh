# 17-1-2018 MRC-Epid JHZ

export HRC=/gen_omics/data/EPIC-Norfolk/HRC/binary_ped/HRC
export glist=/genetics/bin/glist-hg19
export wd=$PWD

code_to_setup() {
# isolate associate genes through galaxy
sort -k1,1n -k2,2n glist-hg19 | \
awk '{if(NR==1) print "#chrom","start","end","gene";print "chr" $1,$2,$3,$4}' OFS="\t" > glist-hg19.bed
awk 'NR>1{
   if (NR==2) print "#chrom","Start","End","oxfordid","snpid","rsid"
   print "chr" $1,$3-1,$3,$2,$15,$16
}' FS="\t" OFS="\t" $1.jma.out > $1.jma.bed
# now we reformat galaxy data
echo chr pos oxfordid snpid rsid chrom start end gene > $1.jma.gene
awk 'NR>1{
  gsub(/chr/,"",$1)
  $2=""
  print
}' OFS="\t" $1.jma.galaxy | awk '{$1=$1};1' | sort -k1,1n -k2,2n >> $1.jma.gene
sed -i 's/ /\t/g' $1.jma.gene
}

awk '$NF!="."{print $NF}' $1.jma.gene > $1.genelist
awk -vFS="\t" -vOFS="\t" 'NR>1{print $2,$13}' $1.jma.out > $1.snpandp
awk '{if($1=="X") $1=23;if($1=="Y") $1=24;if($1!="XY") print}' $glist > glist
vegas2v2 -G -snpandp $1.snpandp -custom $HRC -glist glist -genelist $1.genelist -out $1

#######################################################################################
# Notes on setting up VEGAS2v2:
#
# 0. assumming that $1.out is available from FM-pipeline
# 1. # PLINK2 alpha does not support --noweb and we switch to 1.9 beta instead
#    sed -i 's/plink2/plink-1.9/g' vegasv2
# 2. We also reformat glist since VEGASv2 does not recognise X, Y, XY
# 3. it is very computer-intensive with -glist, i.e., 
#    vegas2v2 -G -snpandp $1.snpandp -custom $HRC -glist glist -out $1
#    so we use -genelist which is perhaps more pertinent
#    Gene-based analysis has -max=1E6 resampling by default
#    Pathway-based analysis setup which also has -maxsample=1E6 resampling by default
#    awk '{print $2,$8}' $1.out | grep -v Gene | sed 's/"//g' $1.geneandp
#    vegas2v2 -P -geneandp $1.geneandp -geneandpath $1.vegas2pathSYM -glist glist
#    The pathway definitions somehow require some work on other databases
#######################################################################################

rm glist
