# 24-1-2018 MRC-Epid JHZ

export PATH=$PATH:/genetics/data/software/bin
export bfile=/gen_omics/data/EPIC-Norfolk/HRC/binary_ped/HRC
export glist=/genetics/bin/glist-hg19
export wd=$PWD

# isolate associate genes through Galaxy
sort -k1,1n -k2,2n glist-hg19 | \
awk '{if(NR==1) print "#chrom","start","end","gene";print "chr" $1,$2,$3,$4}' OFS="\t" > glist-hg19.bed
awk 'NR>1{
   if (NR==2) print "#chrom","Start","End","oxfordid","snpid","rsid"
   print "chr" $1,$3-1,$3,$2,$15,$16
}' FS="\t" OFS="\t" $1.jma.out > $1.jma.bed
# Once uploaded to https://usegalaxy.org, we can reformat data from "Operate on Genomic Intervals",
# "Join the intervals of two datasets side-by-side", "All records of first dataset (fill null with ".")"

# As it is interactive access, we use bedtools instead:
intersectBed -a $1.jma.bed -b glist-hg19.bed -loj > $1.jma.bedtools
echo chr pos oxfordid snpid rsid chrom start end gene > $1.jma.gene
awk '{
  gsub(/chr/,"",$1)
  $2=""
  print
}' OFS="\t" $1.jma.bedtools | awk '{$1=$1};1' | sort -k1,1n -k2,2n >> $1.jma.gene
sed -i 's/ /\t/g;s/-1/./g' $1.jma.gene

awk 'NR>1 && $NF!="."{print $NF}' $1.jma.gene | sort -k1,1 | uniq > $1.genelist
awk -vFS="\t" -vOFS="\t" 'NR>1{print $2,$13}' $1.jma.out > $1.snpandp
awk '{if($1=="X") $1=23;if($1=="Y") $1=24;if($1!="XY") print}' $glist > glist
vegas2v2 -G -snpandp $1.snpandp -custom $bfile -glist glist -genelist $1.genelist -out $1

#######################################################################################
# Notes on setting up VEGAS2v2:
#
# 1. # PLINK2 alpha does not support --noweb and we switch to 1.9 beta instead
#    sed -i 's/plink2/plink-1.9/g' vegas2v2
#    and mask the following statement: system("rm -R $time_stamp");
# 2. We reformat glist-hg19 since VEGAS2v2 does not recognise X, Y, XY
# 3. it is very computer-intensive with -glist, i.e., 
#    vegas2v2 -G -snpandp $1.snpandp -custom $bfile -glist glist -out $1
#    so we use -genelist which is perhaps more pertinent
#    Gene-based analysis has -max=1E6 resampling by default
#    Gene-specific p values can be obtained from -G option above as follows,

awk '{print $2,$8}' $1.gene-basedoutput | grep -v Gene | sed 's/"//g' > $1.geneandp

#    Onice pathway database is ready, the following command can be called:
#
#    vegas2v2 -P -geneandp $1.geneandp -geneandpath $1.vegas2pathSYM -glist glist
#
#    Pathway-based analysis setup which also has -maxsample=1E6 resampling by default
#######################################################################################

rm glist
