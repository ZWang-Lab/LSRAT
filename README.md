# LSRAT: Longitudinal variant-Set Retrospective Association Test

LSRAT is an R program that performs retrospective variant-set testing for longitudinal traits in GWA studies.


## Installation Instructions:

### Required software and packages

1. R(http://www.r-project.org/)
  
2. R Packages: snpStats, LBRAT, GENESIS, plyr, geepack, lme4, nlme, mvtnorm
  
3. PLINK 1.0 or 1.9 (https://www.cog-genomics.org/plink2)

Install the required R package before you install LSRAT package. Install the **LSRAT** using the following steps.


### Install LSRAT on LUNIX or Mac OSX

```
git clone https://github.com/ZWang-Lab/LSRAT.git

R CMD INSTALL LSRAT

```
Alternatively, one can install it in R using the following code.
### Install LSRAT in R
```
library(devtools)
devtools::install_github(repo = 'ZWang-Lab/LSRAT')

```

## Usage instructions:

Details for functions and examples in the LSRAT are available in the manual (LSRAT-manual.pdf)
