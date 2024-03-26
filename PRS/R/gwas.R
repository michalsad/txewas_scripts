#' @export
plr <- function(pheno.file, covar.file, rds.file, rds.sample.file,
                train.sample.file, binary = FALSE, threads = 1) {
  rds <- bigsnpr::snp_attach(rds.file)
  G <- rds$genotypes

  rds.sample <- as.character(
    data.table::fread(rds.sample.file, data.table = FALSE)[, 2])
  train.sample <- as.character(
    data.table::fread(train.sample.file, data.table = FALSE)[, 2])

  y <- read.plink.file(pheno.file)
  y.na <- rownames(y)[is.na(y)]

  covar <- read.plink.file(covar.file)
  covar.na <- rownames(covar)[rowSums(is.na(covar)) > 0]

  nas <- union(y.na, covar.na)
  n.tmp <- length(train.sample)
  s <- intersect(rownames(y), rownames(covar))
  train.sample <- setdiff(intersect(s, train.sample), nas)
  train.ind <- which(rds.sample %in% train.sample)

  if (length(nas) > 0) {
    print(sprintf(
      "WARNING: %d samples had to be dropped. %d samples left.",
      n.tmp - length(train.ind), length(train.ind)))
  }

  y <- y[rds.sample[train.ind],]
  covar <- bigstatsr::covar_from_df(covar[rds.sample[train.ind],])

  print("Running penalized regression...")

  if (binary) {
    st <- system.time(
      mod <- bigstatsr::big_spLogReg(G, y, ind.train = train.ind, K = 10,
                                     covar.train = covar,
                                     pf.covar = rep(0, ncol(covar)),
                                     power_scale = c(0, 0.5, 1),
                                     power_adaptive = c(0, 0.5, 1.5),
                                     lambda.min.ratio = 1e-6, nlam.min = 30,
                                     n.abort = 2, dfmax = 200e3,
                                     ncores = threads)
    )
  } else {
    st <- system.time(
      mod <- bigstatsr::big_spLinReg(G, y, ind.train = train.ind, K = 10,
                                     covar.train = covar,
                                     pf.covar = rep(0, ncol(covar)),
                                     power_scale = c(0, 0.5, 1),
                                     power_adaptive = c(0, 0.5, 1.5),
                                     lambda.min.ratio = 1e-6, nlam.min = 30,
                                     n.abort = 2, dfmax = 200e3,
                                     ncores = threads)
    ) # 27h for 50% training set and SNP INFO score 1
  }

  print("Done.")
  print(st)

  mod
}
