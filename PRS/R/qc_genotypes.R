#' @export
qc.snp.list.submit <- function(path.to.plink, bgen.file.pattern, bgen.sample,
                               mfi.file.pattern, keep.sample, out.file.prefix,
                               logfile.prefix, info.score.min = 0.8,
                               maf.min = 0.01, geno.max = 0.1, hwe.pval = 1e-10,
                               memory.mb = 8000, threads = 1, time.hours = 10) {
  memory.gb <- memory.mb/1000
  bgen.file.pattern <- sub("\\*", "${SGE_TASK_ID}", bgen.file.pattern)
  mfi.file.pattern <- sub("\\*", "${SGE_TASK_ID}", mfi.file.pattern)
  out.file.prefix <- paste(out.file.prefix, "chr${SGE_TASK_ID}", sep = "_")

  template <- readr::read_file(path.join(
    fs::path_package("PRS", "tools"), "submit_template_qc_snp_list.txt"))
  out <- stringr::str_glue(template)
  submit.file <- "submit_qc_snplist.sh"
  readr::write_file(out, file = submit.file)
  cmd <- paste("qsub", submit.file)

  tryCatch({
    mes <- system(cmd, intern = TRUE, timeout = 30)
    Sys.sleep(1)
    file.remove(submit.file)
  }, error = function(err) {
    message(err)
  }, warning = function(war) {
    message(war)
  })
}


#' @export
qc.genotypes <- function(snplist.file.pattern, bgen.file.pattern, bgen.sample,
                         bgi.dir, mfi.file.pattern, keep.sample,
                         out.file.prefix, threads = 1,
                         rm.snplist.files = FALSE) {
  snp.files <- sprintf(sub("\\*", "%d", snplist.file.pattern), 1:22)
  bgen.files <- sprintf(sub("\\*", "%d", bgen.file.pattern), 1:22)
  mfi.files <- sprintf(sub("\\*", "%d", mfi.file.pattern), 1:22)
  snp.list <- list()

  print("Reading snp lists...")

  for (i in seq_along(snp.files)) {
    snps <- unlist(data.table::fread(snp.files[i], header = FALSE,
                                     data.table = FALSE),
                   use.names = FALSE)
    mfi <- data.table::fread(mfi.files[i], data.table = FALSE)
    rownames(mfi) <- mfi[, 2]
    snps <- sub(":", "_", mfi[snps, 1])
    snp.list[[i]] <- snps
  }

  if (rm.snplist.files) {
    unlink(dirname(snplist.file.pattern), recursive = TRUE)
  }

  print("Done.")

  bsample <- as.character(
    data.table::fread(bgen.sample, skip = 2, data.table = FALSE)[, 2])
  ksample <- as.character(
    data.table::fread(keep.sample, data.table = FALSE)[, 2])
  rinds <- which(bsample %in% ksample)

  if (length(setdiff(ksample, bsample)) > 0) {
    print(sprintf(
      "WARNING: %d samples not found in the genotype file. Using %d common samples",
      length(setdiff(ksample, bsample)), length(rinds)))
  }

  print("Reading BGEN files...")

  st <- system.time(
    rds <- bigsnpr::snp_readBGEN(bgenfiles = bgen.files, list_snp_id = snp.list,
                                 backingfile = out.file.prefix, ind_row = rinds,
                                 bgi_dir = bgi.dir, ncores = threads)
  )

  print("Done.")
  print(st)

  data.table::fwrite(data.frame("FID" = bsample[rinds], "IID" = bsample[rinds]),
                     file = paste(out.file.prefix, "sample", sep = "."),
                     quote = FALSE, sep = "\t", na = NA)
}


#' @export
save.to.bed <- function(rds.file, rds.sample, out.dir) {
  bgen <- bigsnpr::snp_attach(rds.file)
  bgen$map <- tibble::add_column(bgen$map, .after = 3, genetic.dist = 0)
  sample <- data.table::fread(rds.sample, data.table = FALSE)[, 2]
  bgen$fam <- data.frame(family.ID = sample, sample.ID = sample,
                         paternal.ID = 0, maternal.ID = 0, sex = -9,
                         affection = -9)
  out.file <- path.join(out.dir,
                        sub.file.ext(basename(rds.file), new.ext = "bed"))
  system.time(
    bed <- bigsnpr::snp_writeBed(bgen, bedfile = out.file)
  ) # 40 min
}


qc.snp.list <- function(path.to.plink, bgen, bgen.sample, mfi, keep.sample,
                        info.score.min = 0.8, maf.min = 0.01, geno.max = 0.1,
                        hwe.pval = 1e-10, memory = 8000, threads = 1,
                        out.file) {
  cmd <- paste(
    path.to.plink,
    sprintf("--memory %d", memory),
    sprintf("--threads %d", threads),
    sprintf("--bgen %s ref-first", bgen),
    sprintf("--sample %s", bgen.sample),
    sprintf("--keep %s", keep.sample),
    sprintf("--maf %.6f", maf.min),
    sprintf("--extract-col-cond %s 8 2", mfi),
    sprintf("--extract-col-cond-min %.4f", info.score.min),
    sprintf("--geno %.4f", geno.max),
    sprintf("--hwe %.e midp", hwe.pval),
    "--write-snplist",
    sprintf("--out %s", out.file))

  tryCatch({
    mes <- system(cmd, intern = TRUE)
  }, error = function(err) {
    message(err)
  }, warning = function(war) {
    message(war)
  })
}


#' @export
qc.snp.list.parallel <- function(path.to.plink, bgen.files, bgen.sample,
                                 mfi.files, keep.sample, info.score.min = 0.8,
                                 maf.min = 0.01, geno.max = 0.1,
                                 hwe.pval = 1e-10, memory.per.thread = 8000,
                                 threads = 1, temp.dir = "tmp") {
  `%dopar%` <- foreach::`%dopar%`
  n <- length(bgen.files)
  threads.per.chr <- ifelse(threads %/% n == 0, 1, threads %/% n)
  dir.create(temp.dir, showWarnings = FALSE)
  out.files <- sprintf(path.join(temp.dir, "plink%d"), 1:n)

  ncl <- ifelse(threads > n, n, threads)

  print(sprintf("Running %d concurrent job(s); using %d thread(s) per job.",
                ncl, threads.per.chr))

  cl <- parallel::makeForkCluster(ncl)
  doParallel::registerDoParallel(cl)
  system.time({
    foreach::foreach(i = 1:n) %dopar% {
      qc.snp.list(path.to.plink, bgen.files[i], bgen.sample, mfi.files[i],
                  keep.sample, info.score.min, maf.min, geno.max, hwe.pval,
                  memory.per.thread, threads.per.chr, out.files[i])
    }
  })
  parallel::stopCluster(cl)

  out.files <- paste(out.files, "snplist", sep = ".")
  qcd.snps <- unlist(lapply(
    out.files,
    function(x) data.table::fread(x, header = FALSE, data.table = FALSE)))
  names(qcd.snps) <- NULL
  unlink(temp.dir, recursive = TRUE)
  qcd.snps
}


#' @export
qc.genotypes.parallel <- function(path.to.plink, backingfile, bgen.pattern,
                                  bgen.sample, mfi.pattern, keep.sample,
                                  info.score.min = 0.8, maf.min = 0.01,
                                  geno.max = 0.1, hwe.pval = 1e-10,
                                  memory.per.thread = 8000, threads = 1) {
  bgen.files <- sprintf(sub("\\*", "%d", bgen.pattern), 1:22)
  mfi.files <- sprintf(sub("\\*", "%d", mfi.pattern), 1:22)

  print("Calculating snp list...")
  snps <- qc.snp.list.parallel(path.to.plink, bgen.files, bgen.sample,
                               mfi.files, keep.sample, info.score.min, maf.min,
                               geno.max, hwe.pval, memory.per.thread, threads,
                               dirname(backingfile))
  print("Done.")

  bsample <- as.character(
    data.table::fread(bgen.sample, skip = 2, data.table = FALSE)[, 2])
  ksample <- as.character(
    data.table::fread(keep.sample, data.table = FALSE)[, 2])
  rinds <- which(bsample %in% ksample)

  print("Preparing genotype file...")

  if (length(setdiff(ksample, bsample)) > 0) {
    print(sprintf(
      "WARNING: %d samples not found in the genotype file. Using %d common samples",
      length(setdiff(ksample, bsample)), length(rinds)))
  }

  system.time(
    rds <- bigsnpr::snp_readBGEN(bgenfiles = bgen.files, list_snp_id = snps,
                                 backingfile = backingfile, ind_row = rinds,
                                 ncores = threads)
  )
  bgen <- bigsnpr::snp_attach(rds)
  system.time(
    bed <- bigsnpr::snp_writeBed(bgen, bedfile = backingfile)
  )

  print("Done.")
}


#' @export
read.bed <- function(bedfile, ...) {
  print("Reading the file...")
  rds <- bigsnpr::snp_readBed(bedfile, ...)
  print("Done.")

  dat <- bigsnpr::snp_attach(rds)

  print("Imputing missing values...")
  dat.imp <- bigsnpr::snp_fastImputeSimple(dat$genotypes, method = "random")
  print("Done.")

  dat$genotypes <- dat.imp
  save.sample.file(dat$fam$sample.ID,
                   out.file = sub.file.ext(rds, new.ext = "sample"))
  bigsnpr::snp_save(dat)
}
