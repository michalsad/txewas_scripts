#' @export
sample.split.simple <- function(samples.file, out.dir, test.frac = 0.5) {
  samples <- read.sample.file(samples.file)
  ntest <- round(test.frac*length(samples))
  samples.test <- sample(samples, size = ntest)
  samples.train <- setdiff(samples, samples.test)

  out.file <- path.join(
    out.dir,
    add.suffix(basename(samples.file),
               suffix = sprintf("simple_test%.1f_test", test.frac)))
  save.sample.file(samples.test, out.file)
  out.file <- path.join(
    out.dir,
    add.suffix(basename(samples.file),
               suffix = sprintf("simple_test%.1f_train", test.frac)))
  save.sample.file(samples.train, out.file)
}


#' @export
sample.split <- function(binary.factor.file, keep.sample, out.dir,
                         test.frac = 0.5) {
  fac <- read.plink.file(binary.factor.file)
  mask <- !is.na(fac)
  fac <- matrix(fac[mask,], ncol = 1, dimnames = list(rownames(fac)[mask]))

  ksample <- data.table::fread(keep.sample, data.table = FALSE)
  ksample <- as.character(ksample[, 2])
  csample <- intersect(rownames(fac), ksample)

  if (length(setdiff(ksample, rownames(fac))) > 0) {
    print(sprintf(
      "WARNING: %d samples have NA values or not found in the factor file. Using %d common samples",
      length(setdiff(ksample, rownames(fac))), length(csample)))
  }

  fac <- matrix(fac[csample,], ncol = 1, dimnames = list(csample))

  exposed <- rownames(fac)[fac == 1]
  nonexposed <- rownames(fac)[fac == 0]

  ##########################
  set.seed(1)
  ##########################

  exposed.test.ind <- sample.int(length(exposed),
                                 size = round(test.frac*length(exposed)))
  non.exposed.test.ind <- sample.int(length(nonexposed),
                                     size = round(test.frac*length(nonexposed)))

  exposed.train <- exposed[-exposed.test.ind]
  nonexposed.train <- nonexposed[-non.exposed.test.ind]
  exposed.test <- exposed[exposed.test.ind]
  nonexposed.test <- nonexposed[non.exposed.test.ind]
  train <- c(exposed.train, nonexposed.train)
  test <- c(exposed.test, nonexposed.test)

  dat <- list(exposed.train, nonexposed.train, exposed.test, nonexposed.test,
              train, test)
  name <- c("train_exposed", "train_nonexposed", "test_exposed",
            "test_nonexposed", "train", "test")

  for (i in seq_along(dat)) {
    out <- data.frame("FID" = dat[[i]], "IID" = dat[[i]])
    out.file <- path.join(
      out.dir,
      add.suffix(basename(keep.sample),
                 suffix = sprintf("test%.1f_%s", test.frac, name[i])))
    data.table::fwrite(out, file = out.file, quote = FALSE, sep = "\t", na = NA)
  }
}
