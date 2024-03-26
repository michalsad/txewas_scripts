path.join <- function(dir.name, file.name) {
  paste(trimws(dir.name, which = "right", whitespace = "/"),
        file.name, sep = "/")
}


split.ext <- function(file.name) {
  spl <- strsplit(file.name, split = "\\.")[[1]]
  c(paste(utils::head(spl, -1), collapse = "."), utils::tail(spl, 1))
}


remove.ext <- function(file.name) {
  split.ext(file.name)[1]
}


sub.file.ext <- function(file.name, new.ext) {
  paste(remove.ext(file.name), new.ext, sep = ".")
}


add.suffix <- function(file.name, suffix, sep = "_") {
  spl <- split.ext(file.name)
  paste(paste(spl[1], suffix, sep = sep), spl[2], sep = ".")
}


#' @export
read.plink.file <- function(f, ...) {
  dat <- data.table::fread(f, data.table = FALSE, ...)
  rnames <- dat[, 2]
  cnames <- colnames(dat)[-(1:2)]
  dat <- data.frame(dat[, -(1:2)], stringsAsFactors = FALSE)
  rownames(dat) <- rnames
  colnames(dat) <- cnames
  dat
}


is.auto <- function(x) {
  if (length(x) == 1) {
    x == "auto"
  } else {
    FALSE
  }
}


#' @export
save.plink.file <- function(x, out.file, samples = "auto", cnames = "auto") {
  samples <- if (is.auto(samples)) rownames(x) else samples
  cnames <- if (is.auto(cnames)) colnames(x) else cnames

  x <- data.frame(samples, samples, x, check.names = FALSE,
                  fix.empty.names = FALSE, stringsAsFactors = FALSE)
  colnames(x) <- c("#FID", "IID", cnames)
  data.table::fwrite(x, file = out.file, quote = FALSE, sep = "\t", na = NA)
}


#' @export
my.integer <- function(x) {
  out <- tryCatch({
    as.integer(x)
  }, error = function(cond) {
    stop(cond)
  }, warning = function(cond) {
    return(x)
  })
  return(out)
}


#' @export
read.sample.file <- function(f, ...) {
  as.character(data.table::fread(f, data.table = FALSE, ...)[, 2])
}


#' @export
save.sample.file <- function(x, out.file) {
  data.table::fwrite(data.frame(FID = x, IID = x), file = out.file,
                     quote = FALSE, sep = "\t", na = NA)
  out.file
}
