# FM-pipeline

This is a pipeline for finemapping using GWAS summary statistics, implemented in Bash as a series of steps to furnish an incremental analysis. As 
sketched in the diagram below ![one](files/fm-pipeline.png) the procedure attempts identify causal variants from region(s) showing significant SNP-trait 
association.

The process involves the following steps,
1. Extraction of effect (beta)/z statistics from GWAS summary statistics (.sumstats), 
2. Extraction of correlation from the reference panel among overlapped SNPs from 1 and the reference panel containing individual level data. 
3. Information from 1 and 2 above is then used as input for finemapping.

The measure of evidence is typically (log10) Bayes factor (BF) and associate SNP probability in the causal set.

Name | function | Input | Output | Reference
-----|----------|-------|--------|----------
JAM | Finemapping | beta, individual reference data | Bayes Factor of being causal | Newcombe, et al. (2016)
Finemap | Finemapping | z, correlation matrix | Causal SNPs+configuration | Benner, et al. (2016)
CAVIAR | Finemapping | z, correlation matrix | Causal sets and probabilities | Hormozdiari, et al. (2014)
CAVIARBF | finemapping | z, correlation matrix | BF+probabilities for all configurations | Chen, et al. (2015)
FM-summary | finemapping | .sumstats Association results | GitHub download
GCTA | joint/conditional analysis | .sumstats, reference data | Association results | Yang, et al. (2012)
fgwas | annotation 

Note that JAM requires Java 1.8 so call to Java -jar inside the function needs to 
reflect this, not straightforward with install_github() from devtools but one needs to 
clone the package, modify the R source code and then use 'R CMD INSTALL R2BGLiMS'.

## INSTALLATION

On many occasions, the pipeline takes advantage of the [GNU parallel](http://www.gnu.org/software/parallel/).

Besides (sub)set of software listed in the table above, the pipeline requires [GTOOL](http://www.well.ox.ac.uk/%7Ecfreeman/software/gwas/gtool.html),
[PLINK](https://www.cog-genomics.org/plink2) 1.9, and the companion program LDstore from finemap's websiet need to be installed. 
[LocusZoom](http://locuszoom.sph.umich.edu/) is also helpful with graphics.

The pipeline itself can be installed in the usual way,
```
git clone https://github.com/jinghuazhao/FM-pipeline
```
The setup is in line with summary statistics from consortia where only RSid are given for the fact that their chromosomal position may be changed
over different builds. To remedy this, we use information from UCSC, i.e.,
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp150Common.txt.gz
gunzip -c snp150Common.txt.gz | cut -f2,4,5 | sort -k3,3 > snp150.txt
```
The software included in this pipeline range from descriptive analysis via fgwas, locuszoom, GCTA to those dedicat3ed to finemapping including CAVIAR, 
CAVIARBF, finemap, R2BGLiMS/JAM. An adapted version of FM-summary is also given.

## USAGE

The syntax of pipeline is simply
```
bash fm-pipeline.sh <input>
```
You will need to change the configurations at the beginning of the script before execution.

## Inputs

### GWAS summary statistics and lead SNPs

The **first input file** will be GWAS summary statistics with the following columns,

SNP | A1 | A2 | beta | se | N
-----|----|----|------|----|--
RSid | Effect allele | Other allele | effect estimate | standard error of effect | sample size

The **second input file** is a list of SNPs for which finemapping will be conducted.

No header is required for neither file.

### Reference panel

The pipeline uses a reference panel in a .GEN format, taking into account directions of effect in both the GWAS summary statistics and the reference panel. Its 
development will facilitate summary statistics from a variety of consortiua as with reference panels such as the HRC and 1000Genomes.

A .GEN file is required for each region, named such that chr${chr}_{start}_{end}.gen, together with a sample file. A [utility program in Stata](files/p0.do) is 
provided to generated such files from their whole chromosome counterpart. The SNPinfo.dta.gz has the following information,

chr |        rsid  |       RSnum |    pos |    FreqA2 |    info  | type |  A1  | A2
----|--------------|-------------|--------|-----------|----------|------|------|----
 1  | 1:54591_A_G  | rs561234294 |  54591 |  .0000783 |  .33544  |    0 |   A  |  G  
 1  | 1:55351_T_A  | rs531766459 |  55351 |  .0003424 |   .5033  |    0 |   T  |  A  
... | ... | ... | ... | ... | ... | ... | ... | ... |

Given these, one can do away with Stata and work on a text version for instance SNPinfo.txt.

Optionally, specification of file containing sample to be excluded from the reference panel is also indicated.

## Outputs

The output will involve counterpart(s) from individual software, i.e., .cavibf, 
caviarbf, .snp/.config, .jam/.top

Software | file extension name | Description
---------|---------------------|------------
CAVIAR   | .caviar
CAVIARBF | .caviarbf
finemap  | .snp/.config | The top SNPs with largest log10(BF) and top configurations as with their log10(BF)
JAM      | .jam/.top | the posterior summary table and top models containing selected SNPs

Sometimes it is helpful to examine directions of effects together with the correlation of them, e.g., for use with finemap, the code [here](files/finemap-check.R) 
serves for illustration.

## EXAMPLE

We use GWAS on 2-hr glucose level as reported by the MAGIC consortium, Saxena, et al. (2010). The data is obtained as follows,
```
wget ftp://ftp.sanger.ac.uk/pub/magic/MAGIC_2hrGlucose_AdjustedForBMI.txt
gzip -f MAGIC_2hrGlucose_AdjustedForBMI.txt
gunzip -c MAGIC_2hrGlucose_AdjustedForBMI.txt.gz | awk -vOFS="\t" -vN=15234 '(NR>1){print $1, $2, $3, $5, $6, N}' | sort -k1,1 > 2hrglucose.txt
```
and the command to call is
```
bash fm-pipeline.sh 2hrglucose.txt
```
A list of two SNPs is contained in [2.snps](files/2.snps). The Stata program 
[p0.do](files/p0.do) generates [Extract.sh](files/Extract.sh) excluding SNPs in 
exc3_122844451_123344451.txt and exc3_122881254_123381254.txt.

## SOFTWARE AND REFERENCES

**[FM-summary](https://github.com/hailianghuang/FM-summary)**

Huang H, et al (2017). Fine-mapping inflammatory bowel disease loci to single-variant resolution. Nature 547, 173–178, doi:10.1038/nature22969

**[fgwas](https://github.com/joepickrell/fgwas)**

Pickrell JK (2014) Joint analysis of functional genomic data and genome-wide association studies of 18 human traits. bioRxiv 10.1101/000752

**[GCTA](cnsgenomics.com/software/gcta/)**

Yang J, et al. (2012). Conditional and joint multiple-SNP analysis of GWAS summary statistics identifies additional variants influencing complex traits. Nat Genet 
44:369-375

**[CAVIAR](https://github.com/fhormoz/caviar)**

Hormozdiari F, et al. (2014). Identifying Causal Variants at Loci with Multiple Signals of Association. Genetics, 44, 725–731

**[CAVIARBF](https://bitbucket.org/Wenan/caviarbf)**

Chen W, et al. (2015). Fine Mapping Causal Variants with an Approximate Bayesian Method Using Marginal Test Statistics. Genetics 200:719-736.

Kichaev G, et al (2014). Integrating functional data to prioritize causal variants in statistical fine-mapping studies." PLoS Genetics 10:e1004722;

Kichaev, G., Pasaniuc, B. (2015). Leveraging Functional-Annotation Data in Trans-ethnic Fine-Mapping Studies. Am. J. Hum. Genet. 97, 260–271.

**[finemap](http://www.christianbenner.com/#)**

Benner C, et al. (2016) FINEMAP: Efficient variable selection using summary data from genome-wide association studies. Bioinformatics 32, 1493-1501        

**[JAM](https://github.com/pjnewcombe/R2BGLiMS)**

Newcombe PJ, et al. (2016). JAM: A Scalable Bayesian Framework for Joint Analysis of Marginal SNP Effects. Genet Epidemiol 40:188–201

**MAGIC paper**

Saxena R, et al. (2010). Genetic variation in GIPR influences the glucose and insulin responses to an oral glucose challenge. Nat Genet 42:142-148
