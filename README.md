# GUS-LD

[![gplv3+](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl.html)

Genotyping Uncertainty with Sequencing data - Linkage Disequilibrium

An R package for estimating linkage disequilibrium using low and/or high coverage sequencing data without requiring filtering with respect to read depth.

### Installation:

The easiest way to install GUSMap in R is using the devtools package.

```
install.packages("devtools")
library(devtools)
install_github("tpbilton/GUSbase")
install_github("tpbilton/GUS-LD")
```

Note: Some of the functions are coded in C and therefore an appropriate C compiler is needed for the package to work. For windows OS, Rtools (https://cran.r-project.org/bin/windows/Rtools/) provides a compiler. 

### Tutorials:

The best way to understand how to use GUSLD is to read through the [Introduction](http://htmlpreview.github.io/?https://github.com/AgResearch/GUS-LD/blob/master/inst/doc/Introduction.html) Tutorial contained in the package.

### Development:

There are plans to contiune developing this package. If you have any suggests for improvement or additional features you like to see, I'd suggest posting a question under the Issues. I'm also happy for people to suggest changes via a pull request provided it fits within the favour of the package.

### Citation:

To cite this R package:

Bilton, T.P., McEwan, J.C., Clarke, S.M., Brauning, R., Van Stijn, T.C., Rowe, S.J., & Dodds, K.G. (2018). Linkage disequilibrium estimation in low coverage high-throughput sequencing data. *Genetics*, *209*(2), 289-400. doi:[10.1534/genetics.118.300831](https://doi.org/10.1534/genetics.118.300831)

### Funding:

The initial development of this package was partially funded by the Ministry of Business, Innovation and Employment via its funding of the “Genomics for Production & Security in a Biological Economy” programme (Contract ID C10X1306).
