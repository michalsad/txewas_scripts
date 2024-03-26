suppressMessages(library(dplyr))

treeQTL <- function(pval.tab, tmp.dir, mult.test.adj = "BH", thr = 0.1) {
  multi_trait_pvals <- data.frame(SNP = pval.tab$gene, trait = pval.tab$tissue, 
                                  t_stat = NA, beta = NA, 
                                  "p-value" = pval.tab$pvalue, FDR = NA)
  
  # Threshold and sort appropriately
  multi_trait_pvals <- multi_trait_pvals[multi_trait_pvals$p.value <= thr, ]
  multi_trait_pvals <- multi_trait_pvals[order(multi_trait_pvals$p.value), ]
  
  # Write to disk
  names(multi_trait_pvals)[5] <- "p-value"
  write.table(multi_trait_pvals, 
              paste(trimws(tmp.dir, which = "right", whitespace = "/"), 
                    "multi_trait_pvals.txt", sep = "/"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Get mumber of tests per gene, i.e. number of tissues it was tested on
  n_tests_per_gene=pval.tab %>% group_by(gene) %>% 
    summarise(n_tests=sum(!is.na(pvalue)))
  colnames(n_tests_per_gene)=c("family", "n_tests")
  
  eqtl.out <- paste(trimws(tmp.dir, which = "right", whitespace = "/"), 
                    "multi_trait_pvals.txt", sep = "/")
  assoc.tmp.out <- paste(trimws(tmp.dir, which = "right", whitespace = "/"), 
                         "eAssoc_multi_trait.txt", sep = "/")
  
  # Get list of eGENEs (i.e. GENEs which affect any of the traits)
  eGenes <- TreeQTL::get_eSNPs(n_tests_per_SNP = n_tests_per_gene,
                               m_eqtl_out = eqtl.out, level1 = thr, 
                               level2 = thr, method = mult.test.adj)
  
  eAssociations <- TreeQTL::get_eAssociations(eDiscoveries = eGenes,
                                              n_tests = n_tests_per_gene,
                                              m_eqtl_out = eqtl.out,
                                              out_file = assoc.tmp.out,
                                              by_snp = TRUE)
  
  eAssociations <- select(eAssociations, SNP, gene, p.value, BBFDR)
  colnames(eAssociations) <- c("GENE_ID", "TISSUE", "PVALUE", "BBFDR")
  eAssociations <- data.table::setorder(eAssociations, GENE_ID)
  eAssociations <- eAssociations[eAssociations$GENE_ID != "testGene",]
  colnames(eGenes) <- c("GENE_ID", "FAMILY_PVALUE", "NTISSUES")
  eGenes <- eGenes[eGenes$GENE_ID != "testGene",]
  
  list(eGenes, eAssociations)
}


parser <- argparser::arg_parser(
  "Hierarchical error control with TreeQTL", 
  hide.opts = TRUE)

parser <- argparser::add_argument(
  parser,
  "files",
  help = paste0("File listing paths to files containing TxEWAS associations. ",
                "There should be one file for each tissue considered in the ",
                "analysis"),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "tissues",
  help = paste0("File listing tissue names. The order must match the order ",
                "of files containing TxEWAS associations for these tissues."),
  type = "character"
)

parser <- argparser::add_argument(
  parser,
  "--tmp-dir",
  help = paste0("Directory for storing temporary files."),
  type = "character",
  default = "tmpdir"
)

parser <- argparser::add_argument(
  parser,
  "--threshold",
  help = "Target error rate.",
  nargs = 1,
  default = 0.1
)

parser <- argparser::add_argument(
  parser,
  "--out",
  help = "Name of the output file without the extension.",
  nargs = 1,
  default = "out"
)

argv <- argparser::parse_args(parser)

dir.create(argv$tmp_dir, showWarnings = FALSE)

files <- unlist(read.delim(argv$files, header = FALSE))
tissues <- unlist(read.delim(argv$tissues, header = FALSE))

pval.tab <- list()

for (i in seq_along(files)) {
  dat <- data.table::fread(files[i], sep = "\t", data.table = FALSE)
  mask <- grepl("ADDx.+", dat$TYPE)
  dat <- data.frame(gene = dat$GENEID[mask], tissue = tissues[i], 
                    pvalue = dat$PVALUE[mask], rsq = NA, 
                    stringsAsFactors = FALSE, row.names = NULL)
  dat <- dat[!duplicated(dat$gene),]
  dat <- rbind(dat, data.frame(gene = "testGene", tissue = tissues[i], 
                               pvalue = 1e-12, rsq = NA))
  pval.tab[[i]] <- dat
}

pval.tab <- do.call("rbind", pval.tab)
eGenes <- treeQTL(pval.tab, tmp.dir = argv$tmp_dir, mult.test.adj = "BH", 
                  thr = argv$threshold)

data.table::fwrite(eGenes[[1]], file = paste(argv$out, "egenes", sep = "."),
                   quote = FALSE, sep = "\t", na = NA)
data.table::fwrite(eGenes[[2]], file = paste(argv$out, "assoc", sep = "."),
                   quote = FALSE, sep = "\t", na = NA)

unlink(argv$tmp_dir, recursive = TRUE)
