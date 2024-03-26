#' @export
calcPRS <- function(geno.file, geno.sample.file, mod.file, test.sample.file,
                    threads = 1) {
  geno <- bigsnpr::snp_attach(geno.file)
  G <- geno$genotypes
  geno.sample <- read.sample.file(geno.sample.file)

  test.sample <- read.sample.file(test.sample.file)
  test.inds <- which(geno.sample %in% test.sample)

  mod <- readRDS(mod.file)
  summ <- summary(mod)
  ncovar <- length(summ$beta[[1]]) - ncol(G)

  scores <- stats::predict(mod, G, ind.row = test.inds, ncores = threads,
                           covar.row = matrix(0, length(test.inds), ncovar))
  samp <- geno.sample[test.inds]
  data.frame(score = scores[match(test.sample, samp)], row.names = test.sample)
}


#' #' @export
#' calcPRS2 <- function(geno.file, mod.file, snp.file, test.sample.file = NULL,
#'                      threads = 1) {
#'   geno <- bigsnpr::snp_attach(geno.file)
#'   G <- geno$genotypes
#'
#'   if (!is.null(test.sample.file)) {
#'     test.sample <- read.sample.file(test.sample.file)
#'     test.inds <- which(as.character(geno$fam$sample.ID) %in% test.sample)
#'   }
#'
#'   mod <- readRDS(mod.file)
#'   summ <- summary(mod)
#'   ncovar <- length(summ$beta[[1]]) - ncol(G)
#'
#'   scores <- stats::predict(mod, G, ind.row = test.inds, ncores = threads,
#'                            covar.row = matrix(0, length(test.inds), ncovar))
#'   samp <- geno.sample[test.inds]
#'   data.frame(score = scores[match(test.sample, samp)], row.names = test.sample)
#' }
