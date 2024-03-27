#!/usr/bin/env Rscript

library(PRS)

parser <- argparser::arg_parser("Get PGS statistics", hide.opts = TRUE)

parser <- argparser::add_argument(
  parser,
  "score",
  help = paste0("Tab-delimited file containing a polygenic score. The file ",
                "must have at least three columns: FID, IID, sumPRS, in this ",
                "order. A avgPRS column has to be added if the --avgPRS ",
                "option is used."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "pheno",
  help = paste0("Tab-delimited file containing a response. The first two ",
                "columns must be FID and IID, respectively."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "cov",
  help = paste0("Tab-delimited file containing covariates. The first two ",
                "columns must be FID and IID, respectively. The file must ",
                "have a baseline column."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--sample",
  help = paste0("Tab-delimited file of sample IDs to be used in the ",
                "analysis. The file must contain two columns with family and ",
                "within-family IDs, respectively. If not specified, an ",
                "intersection of samples from score, pheno, and cov ",
                "will be used."),
  type = "character",
  nargs = 1
)

parser <- argparser::add_argument(
  parser,
  "--out",
  help = "Name of the output file.",
  type = "character",
  nargs = 1,
  default = "out.tsv"
)

parser <- argparser::add_argument(
  parser,
  "--avgPRS",
  help = "Use average PGS.",
  flag = TRUE
)


argv <- argparser::parse_args(parser)

score <- read.plink.file(argv$score)
pheno <- read.plink.file(argv$pheno)
covs <- read.plink.file(argv$cov)

if (is.na(argv$sample)) {
  samples <- intersect(rownames(score),
                       intersect(rownames(pheno), rownames(covs)))
  print(sprintf(
    paste0("%d samples present in all input files will be used. If this is ",
           "not desired, please rerun specifing a sample file."),
    length(samples)))
} else {
  samples <- read.sample.file(argv$sample)
}

if (argv$avgPRS) {
  score <- score[samples, "avgPRS"]
} else {
  score <- score[samples, "sumPRS"]
}

pheno <- pheno[samples,]
covs <- bigstatsr::covar_from_df(covs[samples,])
covs.no.baseline <- covs[, colnames(covs) != "baseline"]

r2 <- incR2(pheno, score, covs)
se <- incR2.bootSE(pheno, score, covs)
r2.no.baseline <- incR2(pheno, score, covs.no.baseline)
se.no.baseline <- incR2.bootSE(pheno, score, covs.no.baseline)

mod <- lm(pheno ~ covs + score)
coefs <- coef(summary(mod))
coefs <- c(coefs["score", 1:2], coefs["covsbaseline", 1:2])

mod.no.baseline <- lm(pheno ~ covs.no.baseline + score)
coefs.no.baseline <- coef(summary(mod.no.baseline))
coefs.no.baseline <- coefs.no.baseline["score", 1:2]

res <- matrix(
  c(r2.no.baseline, se.no.baseline, r2, se, coefs.no.baseline, coefs),
  nrow = 1,
  dimnames = list(NULL, c("R2", "R2_SE", "R2_ADJ", "R2_ADJ_SE", "BETA_PRS",
                          "BETA_PRS_SE", "BETA_PRS_ADJ", "BETA_PRS_ADJ_SE",
                          "BETA_BASELINE", "BETA_BASELINE_SE")))

data.table::fwrite(as.data.frame(res), file = argv$out, quote = FALSE,
                   sep = "\t", na = NA)
