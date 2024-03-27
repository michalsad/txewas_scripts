#!/usr/bin/env Rscript

parser <- argparser::arg_parser("Prepare phenotype file", hide.opts = TRUE)

parser <- argparser::add_argument(
  parser,
  "file",
  help = "UKB .tab file with phenotypes.",
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--field",
  help = "Sequence of fields to extract from --file; in the form f.d+.d+.d+",
  type = "character",
  nargs = Inf,
  short = "-c"
)

parser <- argparser::add_argument(
  parser,
  "--cnames",
  help = "Sequence of strings, which should be used as column names.",
  type = "character",
  nargs = Inf,
  short = "-k"
)

parser <- argparser::add_argument(
  parser,
  "--samples",
  help = paste0("Tab-delimited file of sample IDs. The file must contain two ",
                "columns with family and within-family IDs, respectively."),
  type = "character",
  nargs = 1
)

setClass("tuple")
setMethod("coerce", c(from = "ANY", to = "tuple"),
          function(from, to) {
            lapply(from,
                   function(x) {tup <- strsplit(x, split = ":")[[1]];
                   list(as.integer(strsplit(tup[1], split = "-")[[1]]),
                        PRS::my.integer(tup[2]))}
            )
          }
)

parser <- argparser::add_argument(
  parser,
  "--nsd",
  help = paste0("Sequence of tuples of the form: 'i:n'; where i is an ",
                "(1-based) index of a --field, and n is the number ",
                "of standard deviations from the mean above which ",
                "measurements from field i will be discarded (set to NA). ",
                "For fields with 'i:n' unspecified, all measurements ",
                "will be kept. Can use 'i-j:n' to apply the same operation to ",
                "fields i to j."),
  type = "tuple",
  nargs = Inf
)

parser <- argparser::add_argument(
  parser,
  "--trans",
  help = paste0("Sequence of tuples of the form: 'i:s'; where i is an ",
                "(1-based) index of a --field, and s is the symbol of ",
                "transformation to be applied to field i. ",
                "Fields with 'i:s' unspecified will be left untransformed. ",
                "Available transformations are: log-tranformation ('log'), ",
                "quantile-normalization ('qnorm'), ",
                "inverse normal transfomation ('invnorm') and standarization ",
                "('scale'). Can use 'i-j:n' to apply the same operation to ",
                "fields i to j."),
  type = "tuple",
  nargs = Inf
)

parser <- argparser::add_argument(
  parser,
  "--ifactor",
  help = paste0("Interaction factor. Name of a covariate for which ",
                "interactions with other covariates will be calculated"),
  type = "character",
  nargs = 1
)

parser <- argparser::add_argument(
  parser,
  "--out",
  help = "Name of the output file.",
  type = "character",
  nargs = 1,
  default = "covars.tsv"
)

parser <- argparser::add_argument(
  parser,
  "--rm-NA",
  help = "Remove rows (samples) with NA values.",
  flag = TRUE
)

parser <- argparser::add_argument(
  parser,
  "--add1",
  help = "Add 1 to binary (0/1) covariates.",
  flag = TRUE
)

parser <- argparser::add_argument(
  parser,
  "--bdate",
  help = paste0("Substitute fields f.34.0.0 and f.52.0.0 for a birth date ",
                "variable, calculated as (f.34.0.0 - 1900) + ",
                "(f.52.0.0 - 0.5)/12"),
  flag = TRUE
)


argv <- argparser::parse_args(parser)

if (!is.na(argv$samples)) {
  samples <- data.table::fread(argv$samples, data.table = FALSE)
  samples <- as.character(samples[, 2])
} else {
  samples <- NULL
}

dat <- PRS::read.phenotype(argv$file, argv$field, samples)

mask <- apply(dat, 2, function(x) sum(!is.na(unique(x)))) < 2

if (any(mask)) {
  print(sprintf(
    "WARNING: covariates: %s have only one level - not included in the output file.",
    paste(colnames(dat)[mask], collapse = ",")))
}

if (!is.na(argv$nsd)) {
  for (item in argv$nsd) {
    if (length(item[[1]]) == 1) {
      print(c("nsd", colnames(dat)[item[[1]]]))
      dat[, item[[1]]] <- PRS::remove.outliers(dat[, item[[1]]], item[[2]])
    } else if (length(item[[1]]) == 2) {
      print(c("nsd", colnames(dat)[item[[1]][1]:item[[1]][2]]))
      dat[, item[[1]][1]:item[[1]][2]] <- apply(
        dat[, item[[1]][1]:item[[1]][2]],
        MARGIN = 2,
        function(x) PRS::remove.outliers(x, item[[2]]))
    } else {
      stop("Wrong number of indices")
    }
  }
}

if (!is.na(argv$trans)) {
  for (item in argv$trans) {
    if (length(item[[1]]) == 1) {
      print(c("trans", colnames(dat)[item[[1]]]))
      dat[, item[[1]]] <- PRS::transformf(dat[, item[[1]]], item[[2]])
    } else if (length(item[[1]]) == 2) {
      print(c("trans", colnames(dat)[item[[1]][1]:item[[1]][2]]))
      dat[, item[[1]][1]:item[[1]][2]] <- apply(
        dat[, item[[1]][1]:item[[1]][2]],
        MARGIN = 2,
        function(x) PRS::transformf(x, item[[2]]))
    } else {
      stop("Wrong number of indices")
    }
  }
}

if (argv$add1) {
  for (i in 1:ncol(dat)) {
    if (all(dat[, i] %in% c(0, 1, NA))) {
      dat[, i] <- dat[, i] + 1
    }
  }
}

rnames <- rownames(dat)
cnames <- colnames(dat)
dat <- as.data.frame(dat[, !mask])
rownames(dat) <- rnames
colnames(dat) <- cnames[!mask]

if (argv$bdate) {
  ind.rm <- which(colnames(dat) == "f.52.0.0")
  dat <- PRS::birthdate(dat)
  mask <- mask[-ind.rm]
}

if (!is.na(argv$cnames)) {
  colnames(dat) <- argv$cnames[!mask]
}

if (!is.na(argv$ifactor)) {
  ifactor.index <- which(colnames(dat) == argv$ifactor)
  inters <- dat[, argv$ifactor] * dat[, -ifactor.index]
  colnames(inters) <- paste(argv$ifactor, colnames(dat)[-ifactor.index],
                            sep = "x")
  dat <- cbind(dat, inters)
}

if (argv$bdate) {
  mask <- apply(dat, 1, function(x) any(is.na(x)))
  rnames <- rownames(dat)
  cnames <- colnames(dat)
  dat <- as.data.frame(dat[!mask,])
  colnames(dat) <- cnames
  rownames(dat) <- rnames[!mask]
}

PRS::save.plink.file(dat, out.file = argv$out)
