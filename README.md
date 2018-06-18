# GUS-LD Version 0.1.1 
Genotyping Uncertainty with Sequencing data - Linkage Disequilibrium (GUS-LD).

R code that implements the method introduced by Bilton et al. (2018) that estimates linkage disequilibrium using low coverage sequencing data that have been generated using high-throughput sequencing multiplexing methods, such as genotyping-by-sequencing, restriction-site association DNA or exome-capture. The novelty of this methodology is its ability to account for genotyping errors resulting from low read depths (e.g., allelic dropout). 

[![gplv3+](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl.html)

### Data Format:

The key data to use this package are:
- genon: Matrix of genotype calls, which are coded 2 = homozygous genotype for the reference allele, 1 = heterozygous genotype, 0 = homozygous genotype for the alternate allele. The rows must correspond to the number of individuals and columns correspond to the number of SNPs.
- depth_Ref: Matrix of allele counts for the reference allele associated with each genotype call. This matrix corresponds to the matrix of genotype calls (e.g., the i<sup>th</sup> row and j<sup>th</sup> column in depth is the read depth corresponding to the i<sup>th</sup> row and j<sup>th</sup> column in the genon matrix).
- depth_Alt: Matrix of allele counts for the alternate allele associated with each genotype call. This matrix corresponds to the matrix of genotype calls (e.g., the i<sup>th</sup> row and j<sup>th</sup> column in depth is the read depth corresponding to the i<sup>th</sup> row and j<sup>th</sup> column in the genon matrix).

To convert data from VCF format into the required format, use the function 'readVCF' which is available in the script readVCF.R. For example,
```
source("readVCF.R")
data <- readVCF("vcffilename.vcf")
```

### Files:

- [GUS-LD.R](GUS_LD.R): R script which contains a function for computing all the pairwise LD estimates in the data. Note that parallelization is possible in this function but requires that the [foreach](https://cran.r-project.org/web/packages/foreach/index.html) and [doSNOW](https://cran.r-project.org/web/packages/doSNOW/index.html) R packages are installed.
- [GUS-LD-example.R](GUS_LD-example.R): R script with some example code using the deer dataset found in Bilton et al. (2017).
- [deer.RData](deer.RData): Matrix of genotype calls and allele counts for reference and alternate alleles for the deer dataset used in Bilton et al. (2017).
- [readVCF.R](readVCF.R): R script with functions for converting a VCF file into the required format for GUS-LD.

### Example:

The best way to understand how GUS_LD works is to run through the file [GUS_LD-example.R](GUS-LD-example.R).

### References:

Bilton, T.P., McEwan, J.C., Clarke, S.M., Brauning, R., Van Stijn, T.C., Rowe, S.J., Dodds, K.G. (2018). Linkage disequilibrium estimation in low coverage high-throughput sequencing data. *Genetics*, *209*(2), 289-400. doi:[10.1534/genetics.118.300831](https://doi.org/10.1534/genetics.118.300831)

### License
GUS-LD is Copyright (C) 2017-2018 AgResearch Ltd., released under the GNU General Public License version 3.
 
