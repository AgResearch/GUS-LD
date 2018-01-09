# GUS-LD 
Genotyping Uncertainty with Sequencing data - Linkage Disequilibrium (GUS-LD).

R code that implements the method introduced by Bilton et al. (2017) that estimates linkage disequilibrium using low coverage sequencing data that have been generated using high-throughput sequencing multiplexing methods, such as genotyping-by-sequencing, restriction-site association DNA or exome-capture. The novelty of this methodology is its ability to account for genotyping errors resulting from low read depths (e.g., allelic dropout). 

[![gplv3+](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl.html)

### Data Format:

The key data to use this package are:
- genon: Matrix of genotype calls, which are coded 2 = homozygous genotype for the reference allele, 1 = heterozygous genotype, 0 = homozygous genotype for the alternate allele. The rows must correspond to the number of individuals and columns correspond to the number of SNPs.
- depth: Matrix of read depths associated with each genotype call. This matrix corresponds to the matrix of genotype calls (e.g., the i<sup>th</sup> row and j<sup>th</sup> column in depth is the read depth corresponding to the i<sup>th</sup> row and j<sup>th</sup> column in the genon matrix).

To convert data from VCF format into the required format, the [KGD](https://www.github.com/AgResearch/KGD) software (Dodds et al., 2015) can be used which also includes some quality control checks.

### Files:

- [GUS_LD.R](GUS_LD.R): R script which contains a function for computing all the pairwise LD estimates in the data. Note that parallelization is possible in this function but requires that the [foreach](https://cran.r-project.org/web/packages/foreach/index.html) and [doSNOW](https://cran.r-project.org/web/packages/doSNOW/index.html) R packages are installed.
- [GUS_LD-example.R](GUS_LD-example.R): R script with some example code using the deer dataset found in Bilton et al. (2017).
- [deer_depth.txt](deer_depth.txt): Matrix of read depths for the deer dataset used in Bilton et al. (2017).
- [deer_genon.txt](deer_genon.txt): Matrix of genotype calls for the deer dataset used in Bilton et al. (2017).
- Conder folder: This folder contains files for running GUS-LD using [HTconder](https://research.cs.wisc.edu/htcondor/) on a HPC system (if available to the user). This is particularly useful for computing all the pairwise LD estimates when the data set size is very large. Note: there is no guarantee that this code will work for your installation of condor as these things are often configured differently. The key files needed are:  
  - [GUS_LD_Condor.sh](Condor/GUS_LD_Condor.sh): The environment variables in this bash script need to be specified. These are:    
    - GENON_NAME: File name of the genon object saved as a .txt file.    
    - DEPTH_NAME: File name of the depth object saved as a .txt file.   
    - DATA_FOLDER: Relative path to the folder where GENON_NAME and DEPTH_NAME are located.  
    - RUN_NAME: Name used for the run. Note that a folder with this name will be created which will contain all the output from the generated from this script.    
    - NBLOCKS: Number of blocks to split the genon and depth matrices into. This allows each block to be run as a separate job on Conder thus giving parallelization.    
    - CONDOR_MEMORY: Amount of memory requested for each condor job.     
    - R_LIBRARY: Path to any additional R libraries that may be needed by the used. Note that foreach and doSNOW are also used here to add further parallelization within each condor job.  
  - [GUS_LD_ComputeMat.sh](Condor/GUS_LD_ComputeMat.sh): Bash script for putting the matrices created from GUS_LD_Condor.sh back together into one matrix. Note: this script should only be run after all the condor jobs have finished. 

### Example:

The best way to understand how GUS_LD works is to run through the file [GUS_LD-example.R](GUS_LD-example.R).

### References:

Bilton, T.P., McEwan, J.C., Clarke, S.M., Brauning, R., Van Stijn, T.C., Rowe, S.J., Dodds, K.G. (2017). Linkage disequilibrium estimation in low coverage high-throughput sequencing data. doi:[10.1101/235937](https://doi.org/10.1101/235937)

Dodds, K.G., McEwan, J.C., Brauning, R., Anderson, R.M., Van Stijn, T.C., Kristj&#225;nsson, T., Clarke, S.M. (2015) Construction of relatedness matrices using genotyping-by-sequencing data. BMC Genomics, **16**:1047. doi:[10.1186/s12864-015-2252-3](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-2252-3)

### License
GUS-LD is Copyright (C) 2017 AgResearch Ltd., released under the GNU General Public License version 3.

