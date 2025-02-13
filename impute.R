path.join <- function(dir.name, file.name) {
  paste(trimws(dir.name, which = "right", whitespace = "/"),
        file.name, sep = "/")
}


sub.ext <- function(x, new.ext = NULL) {
  f <- basename(x)
  spl <- strsplit(f, split = "\\.")[[1]]
  newf <- paste(c(utils::head(spl, -1), new.ext), collapse = ".")
  
  if (f == x) {
    return(newf)
  } else {
    return(path.join(dirname(x), newf))
  }
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


sub.alleles <- function(x) sub("([a-zA-Z0-9:_]+)_([a-zA-Z]+)_([a-zA-Z]+)",
                               "\\1_\\3_\\2", x)


vprod <- function(genotypes, wgts){
  wgts.snpid <- rownames(wgts)
  geno.snpidx <- which(genotypes$map$marker.ID %in%
                         c(wgts.snpid, sub.alleles(wgts.snpid)))
  snpid <- genotypes$map$marker.ID[geno.snpidx]
  mask <- snpid %in% wgts.snpid
  snpid[!mask] <- sub.alleles(snpid[!mask])
  
  if (any(mask)) {
    w <- wgts[snpid,]
    names(w) <- snpid
    w[snpid[mask]] <- -w[snpid[mask]]
    imp.expr <- bigstatsr::big_prodVec(genotypes$genotypes, w,
                                       ind.col = geno.snpidx)
    2*sum(wgts[snpid[mask],]) + imp.expr
  } else {
    bigstatsr::big_prodVec(genotypes$genotypes, wgts[snpid,],
                           ind.col = geno.snpidx)
  }
}


parser <- argparser::arg_parser(
  "Impute gene expression based on genotypes", 
  hide.opts = TRUE)

parser <- argparser::add_argument(
  parser,
  "rds",
  help = paste0("File with extension '.rds' that stores genotypes in a ", 
                "bigSNP object. The corresponding '.bk' and '.sample' files ", 
                "must be in the same directory."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "models",
  help = paste0("File listing paths (one path per line, no header) to ", 
                "expression prediction models. Prediction models must be in ",
                "the .Rdata format. If the name of a model file contains the ",
                "'ENSG' ID of a gene, the ID is used in the output file; ",
                "otherwise, the gene is identified by its position on the ", 
                "list (e.g. GENE4)."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--threshold",
  help = paste0("Generate predictions only for models with cross-validation ", 
                "p-value < threshold."),
  nargs = 1,
  default = 1.0
)

parser <- argparser::add_argument(
  parser,
  "--outdir",
  help = paste0("Directory where the output files will be saved. Output ",
                "files use model file names with '.imp' extension added."),
  nargs = 1,
  default = "./"
)

argv <- argparser::parse_args(parser)

models <- data.table::fread(argv$models, header = FALSE, 
                            data.table = FALSE)[, 1]
samples <- read.sample.file(sub.ext(argv$rds, new.ext = "rds.sample"))
geno <- bigsnpr::snp_attach(argv$rds)

for (i in seq_along(models)) {
  model <- models[i]
  
  if (grepl(".*(ENSG\\d+\\.\\d+).*", basename(model))) {
    gene.name <- sub(".*(ENSG\\d+\\.\\d+).*", "\\1", basename(model))
  } else {
    gene.name <- paste0("GENE", i)
  }
  
  load(model, wgt <- new.env())
  
  if (wgt$cv.performance[2, 1] < argv$threshold){
    rownames(wgt$wgt.matrix) <- stringi::stri_c(wgt$snps[, 1], wgt$snps[, 4],
                                                wgt$snps[, 5], wgt$snps[, 6],
                                                sep = "_")
    
    imp.expr <- vprod(geno, wgt$wgt.matrix)
    
    out.file <- path.join(argv$outdir, paste0(basename(model), ".imp"))
    save.plink.file(imp.expr, out.file = out.file, samples = samples,
                    cnames = gene.name)
  } else {
    warning(sprintf(
      "Model %s did not pass the p-value threshold. No file saved.", 
      basename(model)))
  }
}
