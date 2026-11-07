#!/usr/bin/env Rscript
# =============================================================================
# download_gse40279.R  --  step 1 of 2 for the real-data run (download + cache)
#
# Fetches GSE40279 (Hannum 2013, PMID 23177740; 656 whole-blood 450K samples,
# ages 19-101) ONCE and caches it as R objects, so the processing step
# (preprocess_gse40279.R) never touches the network and can be re-run freely.
#
# Writes into --datadir (default ../data/raw/gse40279, git-ignored):
#   beta.rds         CpGs x samples beta matrix (full cohort, ~2.5 GB in RAM)
#   pheno.rds        data.frame(sample_id, age, sex)
#   beta_test.rds    small age-stratified subset (--test-n samples, default 48)
#   pheno_test.rds   matching phenotype rows
#   GSE40279_series_matrix.txt.gz   GEOquery's own download cache
#
# Re-running is cheap: if beta.rds + pheno.rds exist it only (re)builds the test
# subset; pass --force to download again.
#
# Usage (from src/):
#   Rscript download_gse40279.R [--datadir ../data/raw/gse40279] [--test-n 48] [--force]
# Needs: GEOquery, Biobase, data.table  (~16 GB RAM for the full matrix)
# =============================================================================

suppressMessages({
  need <- c("GEOquery", "Biobase", "data.table")
  miss <- need[!sapply(need, requireNamespace, quietly = TRUE)]
  if (length(miss)) stop(
    "Missing packages: ", paste(miss, collapse = ", "),
    "\nInstall with:\n  install.packages(c('BiocManager','data.table'))\n",
    "  BiocManager::install(c('GEOquery','Biobase'))")
  library(GEOquery); library(Biobase); library(data.table)
})

args    <- commandArgs(trailingOnly = TRUE)
getarg  <- function(flag, default) { i <- which(args == flag); if (length(i)) args[i + 1] else default }
datadir <- getarg("--datadir", file.path("..", "data", "raw", "gse40279"))
test_n  <- as.integer(getarg("--test-n", "48"))
force   <- "--force" %in% args
dir.create(datadir, showWarnings = FALSE, recursive = TRUE)
f_beta  <- file.path(datadir, "beta.rds")
f_pheno <- file.path(datadir, "pheno.rds")

# ---- phenotype parsing ------------------------------------------------------
# GEO characteristics look like "age (y): 67", "gender: F". Find the column whose
# values START with the key (anchored, so e.g. "stage" never matches "age").
charcol <- function(ph, pattern) {
  for (cn in grep("^characteristics", names(ph), value = TRUE))
    if (any(grepl(pattern, ph[[cn]], ignore.case = TRUE))) return(as.character(ph[[cn]]))
  NULL
}
kv <- function(x) trimws(sub("^[^:]*:", "", x))

parse_pheno <- function(ph, sample_ids) {
  a <- charcol(ph, "^\\s*age")
  s <- charcol(ph, "^\\s*(gender|sex)")
  if (is.null(a)) stop("could not find an 'age' characteristics column in the GEO phenotype table")
  age <- suppressWarnings(as.numeric(kv(a)))
  sex <- if (is.null(s)) rep(NA_character_, length(age)) else toupper(substr(kv(s), 1, 1))
  data.frame(sample_id = sample_ids, age = age, sex = sex, stringsAsFactors = FALSE)
}

# ---- 1. download (or reuse cache) --------------------------------------------
if (!force && file.exists(f_beta) && file.exists(f_pheno)) {
  message("[1/2] Cached matrix found in ", datadir, " -> skipping download (use --force to refetch)")
  beta  <- readRDS(f_beta)
  pheno <- readRDS(f_pheno)
} else {
  message("[1/2] Downloading GSE40279 series matrix into ", datadir, " (large; GEOquery caches it)...")
  gse  <- getGEO("GSE40279", GSEMatrix = TRUE, getGPL = FALSE, destdir = datadir)[[1]]
  beta <- as.matrix(exprs(gse))
  if (nrow(beta) == 0) {
    # some GEO series keep the betas only as a supplementary file
    message("    series matrix has no values -> fetching the supplementary beta table")
    sf <- getGEOSuppFiles("GSE40279", baseDir = datadir, makeDirectory = FALSE)
    bf <- grep("beta", rownames(sf), ignore.case = TRUE, value = TRUE)[1]
    if (is.na(bf)) stop("no beta file among supplementary files: ", paste(basename(rownames(sf)), collapse = ", "))
    dt   <- fread(bf)
    beta <- as.matrix(dt[, -1]); rownames(beta) <- dt[[1]]
    # map supplementary column names to GSM ids when they differ
    if (!all(colnames(beta) %in% colnames(gse)) && ncol(beta) == ncol(gse)) colnames(beta) <- colnames(gse)
  }
  pheno <- parse_pheno(pData(gse), colnames(gse))
  pheno <- pheno[match(colnames(beta), pheno$sample_id), ]
  saveRDS(beta, f_beta); saveRDS(pheno, f_pheno)
}
message(sprintf("    beta: %d CpGs x %d samples | age %.0f-%.0f (%d missing) | sex %s",
                nrow(beta), ncol(beta), min(pheno$age, na.rm = TRUE), max(pheno$age, na.rm = TRUE),
                sum(is.na(pheno$age)), paste(names(table(pheno$sex)), table(pheno$sex), sep = "=", collapse = " ")))

# ---- 2. age-stratified test subset --------------------------------------------
# Evenly spaced over the age-sorted cohort, so the quick run still spans 19-101.
message("[2/2] Writing age-stratified test subset (n = ", test_n, ")...")
ok  <- which(!is.na(pheno$age))
ord <- ok[order(pheno$age[ok])]
idx <- ord[unique(round(seq(1, length(ord), length.out = min(test_n, length(ord)))))]
saveRDS(beta[, idx, drop = FALSE], file.path(datadir, "beta_test.rds"))
saveRDS(pheno[idx, ],             file.path(datadir, "pheno_test.rds"))
message(sprintf("    test: %d samples, age %.0f-%.0f", length(idx), min(pheno$age[idx]), max(pheno$age[idx])))

message("\nDONE. Next (from src/):\n",
        "  Rscript preprocess_gse40279.R --test   # quick check on the subset (minutes, low RAM)\n",
        "  Rscript preprocess_gse40279.R          # full cohort")
