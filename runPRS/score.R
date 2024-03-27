#!/usr/bin/env Rscript

parser <- argparser::arg_parser("Calculate PRS", hide.opts = TRUE)

parser <- argparser::add_argument(
  parser,
  "genotype",
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
  "model",
  help = "RDS file containing a fitted regression model.",
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "test_sample",
  help = paste0("Tab-delimited file of sample IDs to be used for prediction ",
                "The file must contain two columns with family and ",
                "within-family IDs, respectively."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--out",
  help = "Name of the output file.",
  nargs = 1,
  default = "scores.tsv"
)

parser <- argparser::add_argument(
  parser,
  "--cores",
  help = "Number of computing cores.",
  nargs = 1,
  default = 1
)


argv <- argparser::parse_args(parser)

bigparallelr::set_blas_ncores(1)
options(bigstatsr.check.parallel.blas = FALSE)

score <- PRS::calcPRS(argv$genotype, argv$sample, argv$model, argv$test_sample,
                      argv$cores)
PRS::save.plink.file(score, out.file = argv$out)
