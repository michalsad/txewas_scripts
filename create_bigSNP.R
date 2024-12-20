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


save.sample.file <- function(x, out.file) {
  data.table::fwrite(data.frame(FID = x, IID = x), file = out.file,
                     quote = FALSE, sep = "\t", na = NA)
  out.file
}


sub.alleles <- function(x) sub("([a-zA-Z0-9:_]+)_([a-zA-Z]+)_([a-zA-Z]+)",
                               "\\1_\\3_\\2", x)


parser <- argparser::arg_parser(
  "Create an easy-to-load bigSNP object from a BGEN file", 
  hide.opts = TRUE)

parser <- argparser::add_argument(
  parser,
  "bgen",
  help = paste0("BGEN genotype file with extension '.bgen'. The ", 
                "corresponding '.bgen.bgi' index file and '.bgen.mfi' ",
                "variant information file must exist."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--bgi-dir",
  help = paste0("Directory where the bgen index file is stored. Default is ",
                "the directory of the bgen file."),
  nargs = 1,
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--mfi-dir",
  help = paste0("Directory where the bgen varint information file is stored. ",
                "Default is the directory of the bgen file."),
  nargs = 1,
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--snps",
  help = paste0("File listing SNP IDs (one per line) in the form ",
                "'<chr>_<pos>_<a1>_<a2>' (e.g. '1_88169_C_T') to be kept in ",
                "the bigSNP. This function assumes that these IDs are ",
                "uniquely identifying variants. If not specified, all SNPs ",
                "from the bgen file will be kept in the bigSNP."),
  nargs = 1,
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--samples",
  help = paste0("Two files containing sample IDs. The first one should list ", 
                "samples from the bgen file. The second one should contain ", 
                "samples to be kept in the bigSNP. The files must contain two ",
                "tab-delimited columns with family and within-family IDs, ",
                "respectively. If not specified, all samples from the bgen ",
                "file will be kept in the bigSNP."),
  nargs = 2,
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--out",
  help = paste0("Path (without extension) for the output files ('.bk' and ",
                "'.rds') that are created for storing the bigSNP object. If ",
                "not specified, path to the bgen file is used."),
  nargs = 1,
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--threads",
  help = "Number of computing threads.",
  nargs = 1,
  default = 1
)

argv <- argparser::parse_args(parser)

if (is.na(argv$bgi_dir)) {
  bgi.dir <- dirname(argv$bgen)
} else {
  bgi.dir <- argv$bgi_dir
}

if (is.na(argv$out)) {
  backing.file <- sub.ext(argv$bgen)
} else {
  backing.file <- argv$out
}

if (length(argv$samples) == 2) {
  bgen.sample <- read.sample.file(argv$sample[1])
  keep.sample <- read.sample.file(argv$sample[2])
  keep.idx <- which(bgen.sample %in% keep.sample)
  
  save.sample.file(bgen.sample[keep.idx], 
                   out.file = paste0(backing.file, ".sample"))
  
  keep.unmatched <- keep.sample[!keep.sample %in% bgen.sample]
  
  if (length(keep.sample) > 0){
    out.unmatched <- paste0(backing.file, ".unmatched.samples")
    data.table::fwrite(list(keep.unmatched), file = out.unmatched, 
                       quote = FALSE, row.names = FALSE, col.names = FALSE)
    warning(sprintf("%d samples not matched! Saved to %s", 
                    length(keep.unmatched), out.unmatched))
  }
} else {
  keep.idx <- NULL
}

if (is.na(argv$mfi_dir)) {
  mfi.file <- paste0(argv$bgen, ".mfi")
} else {
  mfi.file <- path.join(argv$mfi_dir, paste0(basename(argv$bgen), ".mfi"))
}

snps <- data.table::fread(mfi.file, select = 1, data.table = FALSE,
                          showProgress = FALSE)[, 1]

if (grepl("^\\d+:", snps[1])) {
  snps <- stringi::stri_replace_first_fixed(snps, ":", "_")
}

if (!is.na(argv$snps)) {
  snps.keep <- data.table::fread(argv$snps, header = FALSE, data.table = FALSE,
                                 showProgress = FALSE)[, 1]
  matched.alleles <- snps[snps %in% snps.keep]
  switched.alleles <- snps[snps %in% sub.alleles(snps.keep)]
  unmatched <- snps.keep[!snps.keep %in% 
                      c(matched.alleles, sub.alleles(switched.alleles))]

  if (length(unmatched) > 0){
    out.unmatched <- paste0(backing.file, ".unmatched.snps")
    data.table::fwrite(list(unmatched), file = out.unmatched, quote = FALSE, 
                       row.names = FALSE, col.names = FALSE)
    warning(sprintf("%d SNPs not matched! Saved to %s", length(unmatched),
                    out.unmatched))
  }
  
  snps <- c(matched.alleles, switched.alleles)
}

cat("Loading the BGEN file...\n")
bigSNP.file <- bigsnpr::snp_readBGEN(bgenfiles = argv$bgen,
                                     backingfile = backing.file,
                                     list_snp_id = list(snps),
                                     ind_row = keep.idx,
                                     bgi_dir = bgi.dir,
                                     read_as = "dosage",
                                     ncores = argv$threads)
cat("Done.\n")

bigSNP <- bigsnpr::snp_attach(bigSNP.file)
bigSNP$map$chromosome <- as.numeric(bigSNP$map$chromosome)
bigSNP$map$marker.ID <- stringi::stri_c(bigSNP$map$chromosome, 
                                        bigSNP$map$physical.pos, 
                                        bigSNP$map$allele1, 
                                        bigSNP$map$allele2,
                                        sep = "_")
bigSNP <- bigsnpr::snp_save(bigSNP)
cat(paste("Output saved to", bigSNP.file, "\n"))
