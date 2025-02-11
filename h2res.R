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


is.auto <- function(x) {
  if (length(x) == 1) {
    x == "auto"
  } else {
    FALSE
  }
}


save.plink.file <- function(x, out.file, samples = "auto", cnames = "auto") {
  samples <- if (is.auto(samples)) rownames(x) else samples
  cnames <- if (is.auto(cnames)) colnames(x) else cnames
  
  x <- data.frame(samples, samples, x, check.names = FALSE,
                  fix.empty.names = FALSE, stringsAsFactors = FALSE)
  colnames(x) <- c("FID", "IID", cnames)
  data.table::fwrite(x, file = out.file, quote = FALSE, sep = "\t", na = NA)
  out.file
}


save.sample.file <- function(x, out.file) {
  data.table::fwrite(data.frame(FID = x, IID = x), file = out.file,
                     quote = FALSE, sep = "\t", na = NA)
  out.file
}


scale.nonbinary <- function(x, center = TRUE, scale = TRUE) {
  mask <- apply(x, 2, function(elem) all(unique(elem) %in% c(0, 1, NA)))
  x[, !mask] <- scale(x[, !mask], center = center, scale = scale)
  x
}


parser <- argparser::arg_parser(
  "Estimate the heritability of drug response", 
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
  "env",
  help = paste0("Tab-delimited file containing the environmental factor. ", 
                "The first two columns must contain family and within-family ",
                "IDs, respectively."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "geno",
  help = paste0("[prefix] Path to .bed file containing the genotypes ",
                "(without extension '.bed', just the prefix). You need the ",
                "corresponding '.bim' and '.fam' in the same directory."),
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
  "ldak",
  help = paste0("Path to LDAK executable."),
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
                "analysis. Samples with missing values will be discarded."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--binary",
  help = "The phenotype is binary (case/control).",
  flag = TRUE
)

parser <- argparser::add_argument(
  parser,
  "--prevalence",
  help = "The population prevalence for the case/control phenotype.",
  nargs = 1,
  type = "numeric"
)

parser <- argparser::add_argument(
  parser,
  "--out",
  help = "Name of the output file.",
  nargs = 1,
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--tmpdir",
  help = "Directory for storing temporary files.",
  nargs = 1,
  default = "."
)

argv <- argparser::parse_args(parser)


if (argv$binary) {
  tryCatch(error = function(err) stop("You need to specify the prevalence!"), 
           stopifnot(!is.na(argv$prevalence)))
}

tmp.dir <- path.join(
  argv$tmpdir,
  paste("tmp_gxemm", round(stats::runif(1)*1e5), sep = "_"))
dir.create(tmp.dir)

pheno <- read.plink.file(argv$pheno)
env <- read.plink.file(argv$env)
covars <- read.plink.file(argv$cov)

if (is.na(argv$samples)) {
  geno.samples <- read.sample.file(paste(argv$geno, "fam", sep = "."), 
                                   header = FALSE)
  samples <- Reduce(
    intersect, 
    c(lapply(list(pheno, env, covars), rownames), list(geno.samples)))
} else {
  samples <- read.sample.file(argv$samples)
}

pheno <- pheno[samples,]
env <- env[samples,]
covars <- covars[samples,]

na.inds <- which(apply(do.call("cbind", list(pheno, env, covars)), 1,
                       function(x) any(is.na(x))))
warning(paste(length(na.inds), "samples discarded due to missing values!"))

if (length(na.inds) > 0) {
  samples <- samples[-na.inds]
  pheno <- pheno[-na.inds]
  env <- env[-na.inds]
  covars <- covars[-na.inds,]
}

keep.sample.file <- save.sample.file(
  samples,
  out.file = path.join(tmp.dir, "keep.sample"))
grm.file.prefix <- path.join(tmp.dir, "grm_ldak")

cmd <- paste(
  argv$ldak,
  sprintf("--bfile %s", argv$geno),
  sprintf("--calc-kins-direct %s", grm.file.prefix),
  sprintf("--keep %s", keep.sample.file),
  "--ignore-weights YES --power -1 --kinship-raw YES")

tryCatch({
  mes <- system(cmd, intern = TRUE)
  cat(paste(mes, collapse = "\n"))
}, error = function(err) {
  unlink(tmp.dir, recursive = TRUE)
  message(err)
  stop(err)
}, warning = function(war) {
  message(war)
  warning(war)
})

grm.file <- paste0(grm.file.prefix, ".grm.raw")
GRM <- as.matrix(data.table::fread(grm.file, data.table = FALSE))
GRM.sample <- read.sample.file(paste0(grm.file.prefix, ".grm.id"),
                               header = FALSE)
sample.inds <- match(samples, GRM.sample)
GRM <- GRM[sample.inds, sample.inds]

if (argv$binary) {
  mask <- pheno == min(pheno)
  pheno[mask] <- 1
  pheno[!mask] <- 2
} else {
  pheno <- scale(pheno)
}

Z <- cbind(env, 1 - env)
mask <- apply(covars, 2, var) == 0

if (sum(mask) > 0) {
  warning(paste("Covariate(s)", paste(which(mask), collapse = ","),
                "discarded due to being single-valued!"))
}

covars <- scale.nonbinary(covars[, !mask])
covars <- as.matrix(cbind(covars, Z[, 1]))

if (argv$binary) {
  gxemm.out <- GxEMM::GxEMM(y = pheno, Z = Z, X = covars, K = GRM,
                            gtype = "free", etype = "free",
                            ldak_loc = argv$ldak, tmpdir = tmp.dir,
                            binary = TRUE, prev = argv$prevalence)
} else {
  gxemm.out <- GxEMM::GxEMM(y = pheno, Z = Z, X = covars, K = GRM,
                            gtype = "free", etype = "free",
                            ldak_loc = argv$ldak, tmpdir = tmp.dir)
}

v0 <- gxemm.out$sig2g[3]
v1 <- gxemm.out$sig2g[2]
w0 <- gxemm.out$sig2e[2]
w1 <- gxemm.out$sig2e[1]
h2res <- (v0 + v1)/(v0 + v1 + w0 + w1)
h2hom <- gxemm.out$sig2g[1]
pval.res <- c(GxEMM::MVWaldtest(gxemm.out$sig2s[2:3], 
                                gxemm.out$sig2Var[2:3, 2:3]))
pval.hom <- c(GxEMM::Waldtest(gxemm.out$sig2s[1], gxemm.out$sig2Var[1, 1]))

names(h2res) <- NULL
names(h2hom) <- NULL

gxemm.out <- c(
  list(h2hom = h2hom, h2response = h2res, h2hom_pvalue = pval.hom, 
       h2response_pvalue = pval.res),
  gxemm.out)

cat("\nHERITABILITY ESTIMATES:\n")
cat(paste("h2hom:", signif(h2hom, digits = 2), "\n"))
cat(paste("h2hom pvalue:", formatC(pval.hom, format = "e", digits = 2), "\n"))
cat(paste("h2response:", signif(h2res, digits = 2), "\n"))
cat(paste("h2response pvalue:", formatC(pval.res, format = "e", digits = 2),
          "\n"))

if (!is.na(argv$out)) {
  save(gxemm.out, file = argv$out)
}

return(NULL)
