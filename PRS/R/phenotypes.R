add.samples.dataframe <- function(x, samples){
  miss.samples <- setdiff(samples, rownames(x))
  cnames <- colnames(x)
  m <- data.frame(matrix(NA, nrow = length(miss.samples), ncol = ncol(x),
                         dimnames = list(miss.samples, colnames(x))))
  x <- rbind(x, m)
  x <- data.frame(x[samples, ], row.names = samples)
  colnames(x) <- cnames
  x
}


add.samples.matrix <- function(x, samples){
  miss.samples <- setdiff(samples, rownames(x))
  m <- matrix(NA, nrow = length(miss.samples), ncol = ncol(x),
              dimnames = list(miss.samples, colnames(x)))
  x <- rbind(x, m)
  matrix(x[samples, ], ncol = ncol(x), dimnames = list(samples, colnames(x)))
}


add.samples <- function(x, samples) {
  if (is.data.frame(x)) {
    add.samples.dataframe(x, samples)
  } else {
    add.samples.matrix(x, samples)
  }
}


#' @export
read.phenotype <- function(phenotype.file, phenotype.codes, samples = NULL){
  header <- colnames(data.table::fread(phenotype.file, nrows = 0))
  inds <- union(1, which(header %in% phenotype.codes))
  dat <- data.table::fread(phenotype.file, select = inds, quote = "",
                           showProgress = TRUE, data.table = FALSE)
  phenotype.codes <- phenotype.codes[phenotype.codes != "f.eid"]
  rows <- dat[, "f.eid"]
  dat <- data.frame(dat[, phenotype.codes])
  rownames(dat) <- rows
  colnames(dat) <- phenotype.codes

  if (!is.null(samples)) {
    inds <- samples[samples %in% rownames(dat)]
    notfound <- setdiff(samples, inds)

    if (length(notfound) > 0) {
      print(sprintf("WARNING: %d samples not in data set. Setting their phenotype values to NA", length(notfound)))
    }

    dat <- data.frame(dat[inds, ])
    rownames(dat) <- inds
    colnames(dat) <- phenotype.codes
    add.samples(dat, samples)
  } else {
    dat
  }
}


phenotypes.to.file <- function(phenotype.file, phenotype.codes, samples.file,
                               out.file) {
  samples <- data.table::fread(phenotype.file, data.table = FALSE)
  samples <- as.character(samples[, 2])
  dat <- read.phenotype(phenotype.file, phenotype.codes, samples)
  save.plink.file(dat, out.file)
}


#' @export
birthdate <- function(x) {
  bdate <- (x$f.34.0.0 - 1900) + (x$f.52.0.0 - 0.5) / 12
  ind1 <- match("f.34.0.0", colnames(x))
  ind2 <- match("f.52.0.0", colnames(x))
  x[, ind1] <- bdate
  colnames(x)[ind1] <- "bdate"
  rnames <- rownames(x)
  cnames <- colnames(x)[-ind2]
  x <- as.data.frame(x[, -ind2])
  rownames(x) <- rnames
  colnames(x) <- cnames
  x
}


#' @export
remove.outliers <- function(x, nsd) {
  mu <- mean(x, na.rm = TRUE)
  sigma <- stats::sd(x, na.rm = TRUE)
  x[x < mu - nsd*sigma | x > mu + nsd*sigma] <- NA
  x
}


#' @export
invnorm <- function(x, c = 0.5) {
  invn <- stats::qnorm((rank(x, na.last = "keep") - c)/
                         (sum(!is.na(x)) - 2*c + 1))
  matrix(invn, ncol = 1, dimnames = list(rownames(x), colnames(x)))
}


#' @export
transformf <- function(x, trans) {
  d <- list("scale" = scale, "log" = log, "invnorm" = invnorm)
  d[[trans]](as.matrix(x))
}
