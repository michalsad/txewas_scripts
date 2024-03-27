#!/usr/bin/env Rscript

parser <- argparser::arg_parser(
  "Create an RDS file from the input BED file",
  hide.opts = TRUE)

parser <- argparser::add_argument(
  parser,
  "bedfile",
  help = paste0("Path to file with extension '.bed' to read. You need the ",
                "corresponding '.bim' and '.fam' in the same directory."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--out",
  help = paste0("Path (without extension) to the output file. Default takes ",
                "the bedfile without the '.bed' extension."),
  nargs = 1
)

argv <- argparser::parse_args(parser)

if (is.na(argv$out)) {
  PRS::read.bed(argv$bedfile)
} else {
  PRS::read.bed(argv$bedfile, backingfile = argv$out)
}
