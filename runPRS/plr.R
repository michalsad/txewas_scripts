#!/usr/bin/env Rscript

parser <- argparser::arg_parser("Run penalized regression", hide.opts = TRUE)

parser <- argparser::add_argument(
  parser,
  "pheno",
  help = paste0("Tab-delimited file containing responses. The column with ",
                "responses must be preceeded by two columns containing ",
                "family and within-family IDs, respectively."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "covar",
  help = paste0("Tab-delimited file containing covariates. The first two ",
                "columns of the file must contain family and within-family ",
                "IDs, respectively."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "geno",
  help = "RDS file containing genotypes.",
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "sample",
  help = paste0("Tab-delimited file of sample IDs accompanying the genotype ",
                "file. The file must contain two columns with family and ",
                "within-family IDs, respectively."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "train_sample",
  help = paste0("Tab-delimited file of sample IDs to be used for fitting ",
                "the model. The file must contain two columns with family and ",
                "within-family IDs, respectively."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--out",
  help = "Name of the output file.",
  nargs = 1,
  default = "plr.rds"
)

parser <- argparser::add_argument(
  parser,
  "--cores",
  help = "Number of computing cores.",
  nargs = 1,
  default = 1
)

parser <- argparser::add_argument(
  parser,
  "--binary",
  help = "Fit penalized logistic regression.",
  flag = TRUE
)


argv <- argparser::parse_args(parser)

bigparallelr::set_blas_ncores(1)
options(bigstatsr.check.parallel.blas = FALSE)

mod <- PRS::plr(argv$pheno, argv$covar, argv$geno, argv$sample,
                argv$train_sample, argv$binary, argv$cores)
saveRDS(mod, argv$out)
