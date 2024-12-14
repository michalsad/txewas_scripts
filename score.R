path.join <- function(dir.name, file.name) {
  paste(trimws(dir.name, which = "right", whitespace = "/"),
        file.name, sep = "/")
}


parser <- argparser::arg_parser(
  "Generate a score file for individual level prediction", 
  hide.opts = TRUE)

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
  "--names",
  help = paste0("File where each row contains a name of the model to be ",
                "selected from the corresponding file listed under 'models'. ",
                "Useful when more than one approach was used in 'models' ",
                "(default: model with the most significant cross-validation ",
                "p-value is selected)."),
  nargs = 1,
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--outdir",
  help = paste0("Directory where the output files will be saved. Output ",
                "files use model file names with '.score' extension added."),
  nargs = 1,
  default = "./"
)

argv <- argparser::parse_args(parser)

models <- data.table::fread(argv$models, header = FALSE, 
                            data.table = FALSE)[, 1]

if (is.na(argv$names)) {
  name.list <- rep("best", length(models))
} else {
  name.list <- data.table::fread(argv$names, header = FALSE, 
                                 data.table = FALSE)[, 1]
}

for (i in seq_along(models)) {
  model <- models[i]
  name <- name.list[i]
  load(model)
  
  if (name == "best") {
    model.idx <- which.min(cv.performance[2,])
    name <- colnames(cv.performance)[model.idx]
  } else {
    if (name %in% colnames(cv.performance)) {
      model.idx <- which(colnames(cv.performance) == name)
    } else {
      stop(sprintf(
        "There is no model named %s in %s. The available models are: %s.", 
        name,
        model,
        paste(colnames(cv.performance), collapse = ", ")))
    }
  }
  
  if (name == "top1") {
    keep <- which.max(wgt.matrix[, model.idx]^2)
  } else {
    keep <- wgt.matrix[, model.idx] != 0
  }
  
  if (grepl(".*(ENSG\\d+\\.\\d+).*", basename(model))) {
    gene.name <- sub(".*(ENSG\\d+\\.\\d+).*", "\\1", basename(model))
  } else {
    gene.name <- paste0("GENE", i)
  }
  
  out <- format(cbind((snps[, c(2, 5, 6)]),  wgt.matrix[, model.idx])[keep,], 
                digits = 3)
  colnames(out) <- c("ID", "A1", "A2", gene.name)
  out.file <- path.join(argv$outdir, paste0(basename(model), ".score"))
  data.table::fwrite(out, file = out.file, quote = FALSE, sep = "\t", na = NA)
}
