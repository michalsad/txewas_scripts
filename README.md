[![DOI](https://zenodo.org/badge/662841914.svg)](https://zenodo.org/badge/latestdoi/662841914)

Code accompanying:

* M. Sadowski, M. Thompson, J. Mefford, T. Haldar, A. Oni-Orisan, R. Border, A. Pazokitoroudi, N. Cai, J. F. Ayroles, S. Sankararaman, A. W. Dahl, and N. Zaitlen, [Characterizing the genetic architecture of drug response using gene-context interaction methods](https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00351-3), Cell Genomics, 2024.

## `h2res.R`: Estimate the heritability of drug response
Estimates the SNP heritability of the phenotype change after treatment ($h^2_\text{response}$) from single time point measurements of this phenotype in drug users and non-users. Utilizes [GxEMM](https://github.com/andywdahl/gxemm). Run: `Rscript h2res.R -h` to see usage information.

## `txewas.R`: Gene-environment interaction test
Given a phenotype, an environmental factor, and imputed expression of a gene for a set of individuals, tests if the interaction between the expression of this gene and the environmental variable is associated with the phenotype. Run: `Rscript txewas.R -h` to see usage information.

## `hFDR.R`: Hierarchical FDR control
Given TxEWAS associations in multiple tissues, performs the hierarchical FDR correction. Uses [TreeQTL](http://bioinformatics.org/treeqtl). Run: `Rscript hFDR.R -h` to see usage information.

## `create_bigSNP.R`: Create a bigSNP object from a BGEN file
Reads in a BGEN file and creates a bigSNP object that enables fast access to imputed allele dosages from R. Uses [bigsnpr](https://privefl.github.io/bigsnpr/index.html). Run: `Rscript create_bigSNP.R -h` to see usage information.

## `impute.R`: Impute gene expression based on genotypes
Reads in a gene expression prediction model and imputes gene expression into genotypes from a bigSNP object. Run: `Rscript impute.R -h` to see usage information.

## `score.R`: Generate a score file for individual level prediction
Reads in a gene expression prediction model and generates a score file for individual level prediction with PLINK. Run: `Rscript score.R -h` to see usage information.

