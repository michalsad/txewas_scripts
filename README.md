[![DOI](https://zenodo.org/badge/662841914.svg)](https://zenodo.org/badge/latestdoi/662841914)

This repository contains scripts used to perform the analysis presented in:

* M. Sadowski, M. Thompson, J. Mefford, T. Haldar, A. Oni-Orisan, R. Border, A. Pazokitoroudi, J. F. Ayroles, S. Sankararaman, A. Dahl, and N. Zaitlen, “Characterizing the genetic architecture of drug response using gene-context interaction methods”, medRxiv, 2024. https://www.medrxiv.org/content/10.1101/2024.03.29.24305093v1.

## `txewas.R`: Gene-environment interaction test
Given a phenotype, an environmental factor, and imputed expression of a gene for a set of individuals, tests if the interaction between the expression of this gene and the environmental variable is associated with the phenotype. Run: `Rscript txewas.R -h` to see usage information. To compute imputed expression for your data, please refer to the TWAS website: http://gusevlab.org/projects/fusion.

## `hFDR.R`: Hierarchical error control with TreeQTL
Given TxEWAS associations in multiple tissues, performs the hierarchical FDR correction. Run: `Rscript hFDR.R -h` to see usage information. It assumes that all TxEWAS associations for a given tissue are saved in one file.

## `PRS`: Custom R package for calculating PRS
The package can be installed using the R package manager: `install.packages("PRS_0.0.0.9000.tar.gz")`. Scripts in `runPRS` were used to work with the `PRS` package&mdash;usage information for every script can be accessed with the `-h` option.

For the heritability analysis tools, please see: https://github.com/andywdahl/gxemm.

