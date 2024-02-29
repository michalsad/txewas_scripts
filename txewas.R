source("covar_prep_utils.R")

args <- commandArgs(trailingOnly = TRUE)

logistic <- as.logical(args[1])
pheno.file <- args[2]
env.file <- args[3]
covars.file <- args[4]
samples.file <- args[5]
rds.file <- args[6]
rm.drugs.file <- args[7]
out.file.prefix <- args[8]
bgen.samples.file <- args[9]
ninter <- args[10]
nthreads <- as.integer(args[11])
stp <- as.integer(args[12])
I <- as.integer(args[13])

RhpcBLASctl::blas_set_num_threads(nthreads)


# if (inter.ind.start == 0 | inter.ind.end == 0) {
#   inter.inds <- NULL
# } else {
#   inter.inds <- inter.ind.start:inter.ind.end
# }

samples <- data.table::fread(samples.file, sep = "\t", header = TRUE,
                             data.table = FALSE)
samples <- as.character(samples[, 2])

bgen.samples <- data.table::fread(bgen.samples.file, sep = " ", header = TRUE,
                                  data.table = FALSE)
bgen.samples <- as.character(bgen.samples[2:nrow(bgen.samples), 2])

samples <- intersect(samples, bgen.samples)

env <- read.plink.file(env.file)

if (toupper(rm.drugs.file) != "NULL") {
  rm.drugs <- read.plink.file(rm.drugs.file)
  rm.drugs <- rm.drugs[, colnames(rm.drugs) != colnames(env)]
  rm.samples <- rownames(rm.drugs)[apply(rm.drugs, 1, function(x) any(x == 1))]
  samples <- setdiff(samples, rm.samples)
}

env <- matrix(env[samples,], ncol = 1, dimnames = list(samples, colnames(env)))

covars <- read.plink.file(covars.file)
covars <- covars[samples,]
csum <- colSums(covars, na.rm = TRUE)
b <- apply(covars, 2, function(x) all(x %in% c(0, 1, NA)))
mask <- b & csum < 20
covars <- covars[, !mask]
covars <- scale.nonbinary(covars)

if (toupper(ninter) == "NULL") {
  ninter <- ncol(covars)
} else {
  ninter <- as.integer(ninter)
}

xenv <- sweep(as.matrix(covars)[, 1:ninter], MARGIN = 1, env[,], `*`)
colnames(xenv) <- paste(colnames(env), colnames(xenv), sep = "x")
ncovars <- ncol(covars)
covars <- do.call("cbind", list(covars, env, xenv))
ncovars <- ncovars + 1

pheno <- read.plink.file(pheno.file)
pheno <- matrix(pheno[samples,], ncol = 1, 
                dimnames = list(samples, colnames(pheno)))

bgen <- bigsnpr::snp_attach(rds.file)

minmax <- data.table::fread(gsub(".rds$", ".minmax", rds.file), sep = "\t",
                            header = TRUE, data.table = FALSE)
geneids <- minmax[, 1]
minmax <- as.matrix(minmax[, 2:3])
rownames(minmax) <- geneids

res <- list()

for (i in I:min(I + stp - 1, nrow(bgen$map))) {
  gname <- bgen$map$rsid[i]
  expr <- matrix(
    descale01(bgen$genotypes[, i], minmax[gname, 1], minmax[gname, 2]),
    ncol = 1, dimnames = list(bgen.samples, gname))
  expr <- matrix(expr[samples,], ncol = 1, dimnames = list(samples, gname))
  expr <- scale(expr)
  
  cxg <- sweep(as.matrix(covars[, 1:ncovars]), MARGIN = 1, expr, `*`)
  colnames(cxg) <- paste(gname, colnames(covars)[1:ncovars], sep = "x")
  C <- do.call("cbind", list(as.matrix(covars), cxg, expr))
  
  if (logistic) {
    g <- glm(pheno ~ C, family = binomial(link = "logit"))
  } else {
    g <- lm(pheno ~ C)
  }
  
  sandwich_se <- diag(sandwich::vcovHC(g, type = "HC"))^0.5
  sandwich_t <- coef(summary(g))[, 1]/sandwich_se
  sandwich_p <- pchisq(sandwich_t^2, 1, lower.tail = FALSE)
  
  coefs <- coef(summary(g))
  coefs[, 2] <- sandwich_se
  coefs[, 3] <- sandwich_t
  coefs[, 4] <- sandwich_p
  
  name <- c(gname, colnames(env), paste(gname, colnames(env), sep = "x"))
  coefs <- coefs[paste0("C", name),]
  coefs <- data.frame(
    rep(gname, 3),
    c("ADD", colnames(env), paste("ADD", colnames(env), sep = "x")), 
    coefs, 
    rep(length(fitted(g)), 3), 
    check.names = FALSE, 
    fix.empty.names = FALSE)
  rownames(coefs) <- NULL
  colnames(coefs) <- c("GENEID", "TYPE", "BETA", "STDERR", "TVALUE", "PVALUE",
                       "N")

  res[[i]] <- coefs
}

res <- do.call("rbind", res)

out.file <- sprintf(
  "%s.%s.on.%s.%s.%d.glmsve", 
  out.file.prefix, 
  colnames(env), 
  colnames(pheno), 
  gsub("ukb.imp.expr.(.*).rds", "\\1", 
       tail(strsplit(rds.file, split = "/")[[1]], 1)), 
  I)
data.table::fwrite(res, file = out.file, quote = FALSE, sep = "\t",
                   row.names = FALSE, col.names = FALSE)
