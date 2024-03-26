#' @export
incR2 <- function(y, x, z) {
  if (typeof(z) == "list") {
    z <- bigstatsr::covar_from_df(z)
  }

  if (typeof(x) == "list") {
    x <- bigstatsr::covar_from_df(x)
  }

  mod <- stats::lm(y ~ z)
  r2 <- summary(mod)$adj.r.squared
  mod <- stats::lm(y ~ z + x)
  summary(mod)$adj.r.squared - r2
}


#' @export
incR2b <- function(y, x, z) {
  if (typeof(z) == "list") {
    z <- bigstatsr::covar_from_df(z)
  }

  if (typeof(x) == "list") {
    x <- bigstatsr::covar_from_df(x)
  }

  mod <- rms::lrm(y ~ z)
  r2 <- mod$stats["R2"]
  mod <- rms::lrm(y ~ z + x)
  mod$stats["R2"] - r2
}


sel <- function(x, inds) {
  if (is.null(dim(x))) {
    x[inds]
  } else {
    x[inds,]
  }
}


bootR2 <- function(y, x, z, n=1e3) {
  if (is.null(dim(x))) {
    l <- length(x)
  } else {
    l <- nrow(x)
  }

  r2boot <- list()

  for (i in 1:n) {
    s <- sample(1:l, replace = TRUE)
    xboot <- sel(x, s)
    yboot <- sel(y, s)
    zboot <- z[s,]
    r2boot[[i]] <- incR2(yboot, xboot, zboot)
  }

  unlist(r2boot)
}


bootR2b <- function(y, x, z, n=1e3) {
  if (is.null(dim(x))) {
    l <- length(x)
  } else {
    l <- nrow(x)
  }

  r2boot <- list()

  for (i in 1:n) {
    s <- sample(1:l, replace = TRUE)
    xboot <- sel(x, s)
    yboot <- sel(y, s)
    zboot <- z[s,]
    r2boot[[i]] <- incR2b(yboot, xboot, zboot)
  }

  unlist(r2boot)
}


#' @export
incR2.bootSE <- function(y, x, z, n=1e3) {
  stats::sd(bootR2(y, x, z, n=n))
}


#' @export
incR2b.bootSE <- function(y, x, z, n=1e3) {
  stats::sd(bootR2b(y, x, z, n=n))
}
