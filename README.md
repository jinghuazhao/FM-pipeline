# FM-pipeline
A pipeline for GWAS finemapping

This is a collection of functions for an incremental approach to GWAS finemapping involving summary statistics.

## Input

The input will be GWAS summary statistics with the following columns,

SNP | A1 | A2 | beta | se | N
-----|----|----|------|----|--
RSid | Effect allele | Other allele | effect estimate | standard error of effect | sample size

## Output

The output will involve output from a variety of software whose list is given below.

## Setup

This is in line with summary statistics from consortia where only RSid are given for the fact that their chromosomal position may be changed over different builds. TO remedy this, we use information from UCSC.
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp150Common.txt.gz
gunzip -c snp150Common.txt.gz | cut -f2,4,5 | sort -k3,3 > snp150.txt
```
The software included in this pipeline range from descriptive analysis via fgwas, locuszoom, GCTA to those dedicat3ed to finemapping including CAVIAR, CAVIARBF, finemap, R2BGLiMS/JAM. An adapted version of FM-summary is also given.

## References

**[FM-summary](https://github.com/hailianghuang/FM-summary)**

Huang H, et al (2017). Fine-mapping inflammatory bowel disease loci to single-variant resolution. Nature 547, 173–178, doi:10.1038/nature22969

**[fgwas](https://github.com/joepickrell/fgwas)**

Pickrell JK (2014) Joint analysis of functional genomic data and genome-wide association studies of 18 human traits. bioRxiv 10.1101/000752

**[CAVIAR](https://github.com/fhormoz/caviar)""

Hormozdiari F, et al. (2014). Identifying Causal Variants at Loci with Multiple Signals of Association. Genetics, 44, 725–731

**[CAVIARBF](https://bitbucket.org/Wenan/caviarbf)**

Chen W, et al. (2015). Fine Mapping Causal Variants with an Approximate Bayesian Method Using Marginal Test Statistics. Genetics 200:719-736.

Kichaev G, et al (2014). Integrating functional data to prioritize causal variants in statistical fine-mapping studies." PLoS Genetics 10:e1004722;

Kichaev, G., Pasaniuc, B. (2015). Leveraging Functional-Annotation Data in Trans-ethnic Fine-Mapping Studies. Am. J. Hum. Genet. 97, 260–271.

**[finemap](http://www.christianbenner.com/#)**

Benner C, et al. (2016) FINEMAP: Efficient variable selection using summary data from genome-wide association studies. Bioinformatics 32, 1493-1501        

**[JAM](https://github.com/pjnewcombe/R2BGLiMS)**

Newcombe PJ, et al. (2016). JAM: A Scalable Bayesian Framework for Joint Analysis of Marginal SNP Effects. Genet Epidemiol 40:188–201
