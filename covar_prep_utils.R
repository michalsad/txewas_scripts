path.join <- function(dir, base) {
  paste(trimws(dir, which = "right", whitespace = "/"), base, sep = "/")
}


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
    add.samples.dataframe(dat, samples)
  } else {
    dat
  }
}


add.samples.dataframe <- function(dat, samples){
  miss.samples <- setdiff(samples, rownames(dat))
  cnames <- colnames(dat)
  m <- data.frame(matrix(NA, nrow = length(miss.samples), ncol = ncol(dat),
                         dimnames = list(miss.samples, colnames(dat))))
  dat <- rbind(dat, m)
  dat <- data.frame(dat[samples, ], row.names = samples)
  colnames(dat) <- cnames
  dat
}


add.samples.matrix <- function(mat, samples){
  miss.samples <- setdiff(samples, rownames(mat))
  m <- matrix(NA, nrow = length(miss.samples), ncol = ncol(mat),
              dimnames = list(miss.samples, colnames(mat)))
  mat <- rbind(mat, m)
  matrix(mat[samples, ], ncol = ncol(mat),
         dimnames = list(samples, colnames(mat)))
}


add.samples <- function(x, samples) {
  if (is.data.frame(x)) {
    add.samples.dataframe(x, samples)
  } else {
    add.samples.matrix(x, samples)
  }
}


remove.outliers <- function(x, nsd) {
  mu <- mean(x, na.rm = TRUE)
  sigma <- sd(x, na.rm = TRUE)
  x[x < mu - nsd*sigma | x > mu + nsd*sigma] <- NA
  x
}


qnormalize <- function(x) {
  xqnorm <- preprocessCore::normalize.quantiles.use.target(
    x,
    target = rnorm(length(x)),
    subset = !is.na(x))
  rownames(xqnorm) <- rownames(x)
  colnames(xqnorm) <- colnames(x)
  xqnorm
}


invnorm <- function(x, c = 0.5) {
  invn <- stats::qnorm((rank(x, na.last = "keep") - c)/
                         (sum(!is.na(x)) - 2*c + 1))
  matrix(invn, ncol = 1, dimnames = list(rownames(x), colnames(x)))
}


scale.nonbinary <- function(x, center = TRUE, scale = TRUE) {
  mask <- apply(x, 2, function(elem) all(unique(elem) %in% c(0, 1, NA)))
  x[, !mask] <- scale(x[, !mask], center = center, scale = scale)
  x
}


is.auto <- function(x) {
  if (length(x) == 1) {
    x == "auto"
  } else {
    FALSE
  }
}


tosave <- function(x, samples = "auto", cnames = "auto") {
  samples <- if (is.auto(samples)) rownames(x) else samples
  cnames <- if (is.auto(cnames)) colnames(x) else cnames
  
  x <- data.frame(samples, samples, x, check.names = FALSE, 
                  fix.empty.names = FALSE, stringsAsFactors = FALSE)
  colnames(x) <- c("#FID", "IID", cnames)
  x
}


descale01 <- function(x, mmin, mmax) {
  x * (mmax - mmin) + mmin
}


read.plink.file <- function(f) {
  dat <- data.table::fread(f, header = TRUE, data.table = FALSE)
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


file.add.to.name <- function(x, add, sep) {
  name.split <- strsplit(x, split = "\\.")[[1]]
  paste(paste(name.split[-length(name.split)], collapse = "."), 
        paste(add, name.split[length(name.split)], sep = "."),
        sep = sep)
}


get.drug <- function(meds.tab, drug.name, drug.codes, exclude.codes) {
  res <- matrix(0, ncol = 1, nrow = nrow(meds.tab), 
                dimnames = list(rownames(meds.tab), drug.name))
  res[apply(meds.tab, 1, function(x) any(x %in% drug.codes, na.rm = TRUE)),
      ] <- 1
  
  if (!is.null(exclude.codes)) {
    res[apply(meds.tab, 1, function(x) any(x %in% exclude.codes, na.rm = TRUE)),
        ] <- NA
  }
  
  res
}


get.drugs <- function(pheno.file, drug.names, drug.codes, exclude.codes, 
                      samples = NULL, records1 = 0:0, records2 = 0:47) {
  codes <- c(
    sapply(records1, function(x) sprintf("f.20003.%d.%d", x, records2)))
  dat <- read.phenotype(pheno.file, codes, samples)
  
  lapply(seq_along(drug.names), 
         function(x) get.drug(dat, drug.names[[x]], drug.codes[[x]], 
                              exclude.codes[[x]]))
}


save.plink.file <- function(x, out.file, samples = "auto", cnames = "auto") {
  data.table::fwrite(tosave(x, samples, cnames), file = out.file, quote = FALSE, 
                     sep = "\t", na = "NA")
}


save.sample.file <- function(x, out.file) {
  data.table::fwrite(data.frame(FID = x, IID = x), file = out.file,
                     quote = FALSE, sep = "\t", na = NA)
  out.file
}


save.drug.files <- function(pheno.file, samples, drug.names, drug.codes, 
                            exclude.codes, out.dir, records = 0:47) {
  res <- get.drugs(pheno.file, samples, drug.names, drug.codes, exclude.codes, 
                   records)
  
  for (i in seq_along(res)) {
    save.plink.file(res[[i]], 
                    path.join(out.dir, paste(drug.names[[i]], "cov", sep = ".")))
  }
}


split.cathegorical.var <- function(x, samples = "auto", cname = "auto",
                                   rm.largest = FALSE) {
  samples <- if (is.auto(samples)) rownames(x) else samples
  cname <- if (is.auto(cname)) colnames(x) else cname
  
  vals <- unique(as.vector(x))
  vals <- vals[!is.na(vals)]
  res <- matrix(0, nrow = nrow(x), ncol = length(vals), 
                dimnames = list(samples, paste0(cname, seq_along(vals))))
  
  for (i in seq_along(vals)) {
    res[(!is.na(x)) & x == vals[i], i] <- 1
  }
  
  res[is.na(x),] <- NA
  
  if (rm.largest) {
    ncases <- colSums(res, na.rm = TRUE)
    max.ind <- which(ncases == max(ncases))[1]
    res <- res[, -max.ind]
  }
  
  res
}
