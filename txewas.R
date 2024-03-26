read.plink.file <- function(f, ...) {
  dat <- data.table::fread(f, data.table = FALSE, ...)
  rnames <- dat[, 2]
  cnames <- colnames(dat)[-(1:2)]
  dat <- data.frame(dat[, -(1:2)], stringsAsFactors = FALSE)
  rownames(dat) <- rnames
  colnames(dat) <- cnames
  dat
}


read.sample.file <- function(f, ...) {
  as.character(data.table::fread(f, data.table = FALSE, ...)[, 2])
}


scale.nonbinary <- function(x, center = TRUE, scale = TRUE) {
  mask <- apply(x, 2, function(elem) all(unique(elem) %in% c(0, 1, NA)))
  x[, !mask] <- scale(x[, !mask], center = center, scale = scale)
  x
}


parser <- argparser::arg_parser(
  "Gene-environment interaction test", 
  hide.opts = TRUE)

parser <- argparser::add_argument(
  parser,
  "pheno",
  help = paste0("Tab-delimited file containing the phenotype. The first two ", 
                "columns must contain family and within-family IDs, ", 
                "respectively."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "gene",
  help = paste0("Tab-delimited file containing the genetic factor. ", 
                "The first two columns must contain family and within-family ",
                "IDs, respectively."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "env",
  help = paste0("Tab-delimited file containing the environmental factor. ", 
                "The first two columns must contain family and within-family ",
                "IDs, respectively."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "cov",
  help = paste0("Tab-delimited file containing covariates. The first two ", 
                "columns must contain family and within-family IDs, ", 
                "respectively."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--samples",
  help = paste0("Tab-delimited file of sample IDs selected for the analysis. ",
                "The file must contain two columns with family and ", 
                "within-family IDs, respectively. If not specified, ", 
                "an intersection of samples from the phenotype, genetic, ", 
                "environmental factor and covariates files is used in the ", 
                "analysis."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--cov-index",
  help = paste0("Indices of covariates from the covariate file for which ",
                "interaction terms with genetic and environmental factors ", 
                "will be computed and added to the model. If this option is ",
                "not used, interaction terms for all the covariates will be ",
                "added. Indices can be provided as a (1-based) ", 
                "comma-separated list. Dash can be used to specify a range ", 
                "of indices. E.g. 1,3-5,8."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--logreg",
  help = "Perform the logistic regression-based test.",
  flag = TRUE
)

parser <- argparser::add_argument(
  parser,
  "--no-header",
  help = "Remove header from the output file.",
  flag = TRUE
)

parser <- argparser::add_argument(
  parser,
  "--out",
  help = "Name of the output file.",
  nargs = 1,
  default = "out.tsv"
)

parser <- argparser::add_argument(
  parser,
  "--threads",
  help = "Number of computing threads.",
  nargs = 1,
  default = 1
)

argv <- argparser::parse_args(parser)

RhpcBLASctl::blas_set_num_threads(argv$threads)

pheno <- read.plink.file(argv$pheno)
gene <- read.plink.file(argv$gene)
env <- read.plink.file(argv$env)
cov <- read.plink.file(argv$cov)

if (is.na(argv$cov_index)) {
  inds <- 1:ncol(cov)
} else {
  inds <- unlist(sapply(strsplit(argv$cov_index, ",")[[1]], 
                        function(x) {
                          spl <- stringr::str_split(x, "-")[[1]];
                          if (length(spl) == 2) {
                            if (spl[2] == "") {
                              spl[2] = ncol(cov)
                            }
                            as.integer(spl[1]):as.integer(spl[2])
                          } else {
                            as.integer(spl[1])
                          }
                        })
  )
}

samples <- Reduce(intersect, lapply(list(pheno, gene, env, cov), rownames))

if (!is.na(argv$samples)) {
  samples <- intersect(samples, read.sample.file(argv$samples))
}

pheno <- matrix(pheno[samples,], ncol = 1, 
                dimnames = list(samples, colnames(pheno)))
gene <- matrix(gene[samples,], ncol = 1, 
               dimnames = list(samples, colnames(gene)))
env <- matrix(env[samples,], ncol = 1, dimnames = list(samples, colnames(env)))

cnames <- colnames(cov)

if (is.null(cnames)) {
  cnames <- paste0("C", 1:ncol(cov))
}

cov <- as.matrix(cov[samples,])
rownames(cov) <- samples
colnames(cov) <- cnames
sel.names <- cnames[inds]

csum <- colSums(cov, na.rm = TRUE)
b <- apply(cov, 2, function(x) all(x %in% c(0, 1, NA)))
mask <- b & csum < 20
sel.names <- setdiff(sel.names, colnames(cov)[mask])
cov <- cov[, !mask]
cov <- scale.nonbinary(cov)

cxenv <- sweep(cov[, sel.names], MARGIN = 1, env[,], `*`)
colnames(cxenv) <- paste(colnames(env), sel.names, sep = "x")
mask <- sapply(seq_along(sel.names), 
               function(x) all(cov[, sel.names[x]] == cxenv[, x]))
cxenv <- cxenv[, !mask]
cov <- do.call("cbind", list(cov, env, cxenv))

gname <- colnames(gene)
gene <- scale(gene)

sel.names <- c(sel.names, colnames(env))
cxgene <- sweep(cov[, sel.names], MARGIN = 1, gene[,], `*`)
colnames(cxgene) <- paste(gname, sel.names, sep = "x")
mask <- sapply(seq_along(sel.names), 
               function(x) all(cov[, sel.names[x]] == cxgene[, x]))
cxgene <- cxgene[, !mask]
C <- do.call("cbind", list(cov, cxgene, gene))

if (argv$logreg) {
  lmod <- glm(pheno ~ C, family = binomial(link = "logit"))
} else {
  lmod <- lm(pheno ~ C)
}
  
sandwich_se <- diag(sandwich::vcovHC(lmod, type = "HC"))^0.5
sandwich_t <- coef(summary(lmod))[, 1]/sandwich_se
sandwich_p <- pchisq(sandwich_t^2, 1, lower.tail = FALSE)
  
coefs <- coef(summary(lmod))
coefs[, 2] <- sandwich_se
coefs[, 3] <- sandwich_t
coefs[, 4] <- sandwich_p
  
name <- c(gname, colnames(env), paste(gname, colnames(env), sep = "x"))
coefs <- coefs[paste0("C", name),]
coefs <- data.frame(
  rep(gname, 3),
  c("ADD", colnames(env), paste("ADD", colnames(env), sep = "x")), 
  coefs, 
  rep(length(fitted(lmod)), 3), 
  check.names = FALSE, 
  fix.empty.names = FALSE)
rownames(coefs) <- NULL
colnames(coefs) <- c("GENEID", "TYPE", "BETA", "STDERR", "TVALUE", "PVALUE",
                     "N")

data.table::fwrite(coefs, file = argv$out, quote = FALSE, sep = "\t",
                   row.names = FALSE, col.names = !argv$no_header)
