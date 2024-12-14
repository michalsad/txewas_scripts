path.join <- function(dir.name, file.name) {
  paste(trimws(dir.name, which = "right", whitespace = "/"),
        file.name, sep = "/")
}


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
  "expr",
  help = paste0("List of files containing gene predicted expression (one ", 
                "file per line, no header). The files must be tab-delimited ", 
                "and contain family and within-family IDs as the first two ",
                "columns, respectively. Use the --merged flag to pass one ",
                "file with multiple genes as columns."),
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
  help = paste0("Column numbers of covariates for which interaction terms ",
                "with genetic and environmental factors are computed and ",
                "added to the model. Indices can be provided as a (1-based) ", 
                "comma-separated list. Dash can be used to specify a range ", 
                "of columns. E.g. 1,3-5,8. Default: interaction terms for ",
                "all covariates are added."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--merged",
  help = paste0("File 'expr' contains predicted expression for multiple ",
                "genes. Rows are samples, columns are genes. The first two ",
                "columns contain family and within-family IDs, respectively. ",
                "The file is tab-delimited."),
  flag = TRUE
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
  default = "txewas.stats"
)

argv <- argparser::parse_args(parser)

pheno <- read.plink.file(argv$pheno)
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

samples <- rownames(pheno)

if (!is.na(argv$samples)) {
  samples <- intersect(samples, read.sample.file(argv$samples))
}

pheno <- matrix(pheno[samples,], ncol = 1, 
                dimnames = list(samples, colnames(pheno)))
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

if (argv$merged) {
  expr <- read.plink.file(argv$expr)
  expr.files <- paste0(colnames(expr), ".imp")
} else {
  expr.files <- data.table::fread(argv$expr, header = FALSE, 
                                  data.table = FALSE)[, 1]
  expr <- do.call("cbind", lapply(expr.files, read.plink.file))
}

if (argv$logreg) {
  fit <- function(x) glm(x, family = binomial(link = "logit"))
} else {
  fit <- lm
}

out <- list()

for (i in 1:ncol(expr)) {
  gene <- expr[i]
  gname <- colnames(gene)
  gene <- matrix(gene[samples,], ncol = 1, dimnames = list(samples, gname))
  gene <- scale(gene)
  
  sel.names <- c(sel.names, colnames(env))
  cxgene <- sweep(cov[, sel.names], MARGIN = 1, gene[,], `*`)
  colnames(cxgene) <- paste(gname, sel.names, sep = "x")
  mask <- sapply(seq_along(sel.names), 
                 function(x) all(cov[, sel.names[x]] == cxgene[, x]))
  cxgene <- cxgene[, !mask]
  C <- do.call("cbind", list(cov, cxgene, gene))
  
  lmod <- fit(pheno ~ C)
  
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
  
  out[[i]] <- coefs
}

out <- do.call("rbind", out)
colnames(out) <- c("GENEID", "TYPE", "BETA", "STDERR", "TVALUE", "PVALUE", "N")
data.table::fwrite(out, file = argv$out, quote = FALSE, sep = "\t", 
                   row.names = FALSE, col.names = !argv$no_header)
