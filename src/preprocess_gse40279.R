#!/usr/bin/env Rscript
# =============================================================================
# preprocess_gse40279.R  --  step 2 of 2 for the real-data run (processing)
# Cached GSE40279 beta matrix  ->  model feature matrix (features.csv)
#
# Reads the cache written by download_gse40279.R (no network access here) and
# turns the whole-blood 450K matrix (Hannum 2013, PMID 23177740) into the exact
# per-sample features the Prevention & Longevity Map methylation layer consumes:
# 7 model-input features + age + an EXTERNAL clock age + two validation-only
# EpiScores. Each model feature is a COHORT PERCENTILE (0-1), which is what
# src/methylation_layer.run_patient expects.
#
# TWO MODES, SAME CODE PATH
#   --test   uses the small age-stratified subset (beta_test.rds, ~48 samples):
#            minutes and a few GB RAM -> features_test.csv. Use it to check the
#            R stack, the weight tables and the CpG overlaps before the full run.
#            Percentiles are relative to the subset, so numbers are not final.
#   (none)   full cohort (beta.rds, 656 samples, ~16 GB RAM) -> features.csv
#
# DESIGN NOTE (non-circularity): the epigenetic clock is computed for EXTERNAL
# comparison ONLY. It is never fed back into the model as an input. Likewise the
# CRP/IL-6 and GDF15 EpiScores are written as val_* columns used only to check
# the model's Aging_inflammation / Aging_secretory sub-scores against an
# independent molecular signal.
#
# WHAT RUNS OUT OF THE BOX vs WHAT YOU SUPPLY
#   - Cell composition (EpiDISH, built-in blood reference) ....... out of the box
#   - Epigenetic clock (methylclock: Horvath/Hannum) ............. out of the box
#   - Protein / exposure EpiScores (CRP, IL6, GDF15, adiponectin,
#     leptin, liver, smoking) ................................... need published
#       CpG weight tables in episcore_weights/ (see episcore_weights/README.md).
#   Any EpiScore whose weight table is absent is written as NA and imputed to the
#   cohort median downstream (geo_validation.py reports which were imputed), so a
#   PARTIAL run (cell composition + clock + whatever scores you have) still works.
#
# Usage (from src/, after download_gse40279.R):
#   Rscript preprocess_gse40279.R --test        # quick check on the subset
#   Rscript preprocess_gse40279.R               # full cohort
# Options:
#   --datadir ../data/raw/gse40279   --weights episcore_weights   --outdir .
# Output:
#   features.csv  (features_test.csv with --test); columns documented at the bottom
# =============================================================================

suppressMessages({
  need <- c("EpiDISH", "methylclock", "data.table")
  miss <- need[!sapply(need, requireNamespace, quietly = TRUE)]
  if (length(miss)) stop(
    "Missing packages: ", paste(miss, collapse = ", "),
    "\nInstall with:\n  install.packages(c('BiocManager','data.table'))\n",
    "  BiocManager::install(c('EpiDISH','methylclock'))")
  library(EpiDISH); library(methylclock); library(data.table)
})

# ---- args -------------------------------------------------------------------
args    <- commandArgs(trailingOnly = TRUE)
getarg  <- function(flag, default) { i <- which(args == flag); if (length(i)) args[i + 1] else default }
test    <- "--test" %in% args
# legacy: also feed the CRP/IL6 EpiScores into the `inflammation` input. Off by default because the
# same EpiScores form val_inflammation_episcore, which would make that check circular.
infl_crp <- "--inflammation-with-crp" %in% args
datadir <- getarg("--datadir", file.path("..", "data", "raw", "gse40279"))
wdir    <- getarg("--weights", "episcore_weights")   # = src/episcore_weights when run from src/
outdir  <- getarg("--outdir", ".")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
sfx     <- if (test) "_test" else ""
outfile <- file.path(outdir, paste0("features", sfx, ".csv"))
pctl    <- function(x) { r <- rank(x, na.last = "keep", ties.method = "average"); (r - 1) / (sum(!is.na(x)) - 1) }

# =============================================================================
# 1. Load cached beta matrix + phenotype (from download_gse40279.R)
# =============================================================================
f_beta  <- file.path(datadir, paste0("beta", sfx, ".rds"))
f_pheno <- file.path(datadir, paste0("pheno", sfx, ".rds"))
if (!file.exists(f_beta) || !file.exists(f_pheno))
  stop("No cached data at ", f_beta, "\nRun first:  Rscript download_gse40279.R --datadir ", datadir)
message(sprintf("[1/5] Loading %s cohort from %s ...", if (test) "TEST" else "FULL", datadir))
beta  <- readRDS(f_beta)
pheno <- readRDS(f_pheno)
stopifnot(identical(colnames(beta), pheno$sample_id))
message(sprintf("    beta: %d CpGs x %d samples | age %.0f-%.0f", nrow(beta), ncol(beta),
                min(pheno$age, na.rm = TRUE), max(pheno$age, na.rm = TRUE)))
if (test) message("    TEST MODE: percentiles are relative to this subset; use for checking, not results")

# =============================================================================
# 2. Cell composition (EpiDISH, built-in blood reference) -> inflammation/immune
# =============================================================================
message("[2/5] Estimating leukocyte composition (EpiDISH RPC)...")
data(centDHSbloodDMC.m, package = "EpiDISH")      # 6-cell blood reference (ships with EpiDISH)
cell <- as.data.frame(epidish(beta.m = beta, ref.m = centDHSbloodDMC.m, method = "RPC")$estF)  # CD4T,CD8T,NK,B,Mono,Neutro
nlr  <- cell$Neutro / (cell$CD4T + cell$CD8T + cell$B + cell$NK)   # neutrophil:lymphocyte ratio
cd48 <- cell$CD4T / pmax(cell$CD8T, 1e-6)                          # CD4:CD8 (falls with immune aging)

# =============================================================================
# 3. Protein / exposure EpiScores (user-supplied CpG weight tables)
#    Each table: CSV with columns  cpg,weight  (optionally an intercept row cpg="(Intercept)").
#    Sources: Gadd et al. 2022 eLife (protein EpiScores: CRP, IL6, GDF15, ADIPOQ, LEP, liver);
#             McCartney et al. 2018 (smoking EpiScore). MethylDetectR (Hillary & Marioni) bundles
#             several of these weight tables ready to drop in.
# =============================================================================
message("[3/5] Computing EpiScores from supplied weight tables in '", wdir, "'...")
if (!dir.exists(wdir)) message("    (weights folder not found -> all EpiScores NA; cell composition + clock still run)")
episcore <- function(beta, file) {
  fp <- file.path(wdir, file)
  if (!file.exists(fp)) { message("    - ", file, " MISSING -> NA"); return(rep(NA_real_, ncol(beta))) }
  w   <- fread(fp); setnames(w, 1:2, c("cpg", "weight"))
  interc <- w[cpg %in% c("(Intercept)", "Intercept"), weight]; interc <- if (length(interc)) interc[1] else 0
  w   <- w[!cpg %in% c("(Intercept)", "Intercept")]
  common <- intersect(w$cpg, rownames(beta))
  if (!length(common)) { message("    - ", file, ": 0 CpGs overlap -> NA"); return(rep(NA_real_, ncol(beta))) }
  message(sprintf("    - %-22s %d/%d CpGs", file, length(common), nrow(w)))
  wv <- setNames(w$weight, w$cpg)[common]
  as.numeric(interc + crossprod(beta[common, , drop = FALSE], wv))
}
es_crp    <- episcore(beta, "CRP.csv")
es_il6    <- episcore(beta, "IL6.csv")              # optional: no IL6 EpiScore in MethylDetectR/Gadd -> normally NA
es_gdf15  <- episcore(beta, "GDF15.csv")
es_adipoq <- episcore(beta, "ADIPOQ.csv")
es_leptin <- episcore(beta, "LEP.csv")
es_liver  <- episcore(beta, "liver.csv")           # ALT/GGT composite, or a single hepatic EpiScore
es_smoke  <- episcore(beta, "smoking.csv")
es_nad    <- episcore(beta, "NAD_CD38.csv")        # optional; no standard table exists

# =============================================================================
# 4. Epigenetic clock  -- EXTERNAL comparison ONLY, never a model input
# =============================================================================
message("[4/5] Computing epigenetic clock (methylclock)...")
clock_age <- tryCatch({
  cl <- DNAmAge(beta, clocks = c("Horvath", "Hannum"))
  if ("Hannum" %in% names(cl)) as.numeric(cl$Hannum) else as.numeric(cl[[2]])
}, error = function(e) { message("    clock failed (", conditionMessage(e), ") -> NA"); rep(NA_real_, ncol(beta)) })

# =============================================================================
# 5. Assemble model features as COHORT PERCENTILES (0-1) and write features.csv
#    Feature construction (documented, blood-aging directions):
#      inflammation  = percentile of neutrophil:lymphocyte (EpiDISH)          (higher = more inflamed)
#                      (--inflammation-with-crp: mean percentile of CRP, IL6, N:L -- legacy, makes the
#                       inflammation check below circular)
#      immune_aging  = percentile of LOW CD4:CD8 (ratio falls with age)     (higher = more aged)
#      nad_depletion = NA unless a CD38/NAD surrogate score is supplied      (imputed median downstream)
#      hepatic       = percentile of the liver EpiScore
#      adiposity_lep = percentile of the leptin EpiScore
#      adiponectin   = percentile of the adiponectin EpiScore (protective; model treats as ADIPOQ input)
#      damage_smoke  = percentile of the smoking EpiScore
#    Validation-only (NOT model inputs):
#      val_inflammation_episcore = mean percentile of CRP (+IL6 if a table exists) -> checks
#                                  Aging_inflammation; independent of the inputs unless
#                                  --inflammation-with-crp
#      inflammation_def          = which inflammation definition was used ("NLR" / "NLR+CRP")
#      val_gdf15_episcore        = percentile of GDF15          -> checks Aging_secretory
# =============================================================================
message("[5/5] Building percentile features + writing ", basename(outfile), "...")
if (infl_crp) message("    --inflammation-with-crp: inflammation includes CRP/IL6, so the ",
                      "Aging_inflammation vs CRP check is NOT independent")
meanp <- function(...) {
  v <- rowMeans(apply(cbind(...), 2, pctl), na.rm = TRUE)
  v[is.nan(v)] <- NA                                # all inputs missing -> NA (not NaN)
  v
}

feat <- data.table(
  sample_id     = pheno$sample_id,
  age           = pheno$age,
  sex           = pheno$sex,
  # ---- 7 model-input features (cohort percentiles 0-1) ----
  inflammation  = if (infl_crp) meanp(es_crp, es_il6, nlr) else pctl(nlr),
  immune_aging  = pctl(-cd48),
  nad_depletion = pctl(es_nad),
  hepatic       = pctl(es_liver),
  adiposity_lep = pctl(es_leptin),
  adiponectin   = pctl(es_adipoq),
  damage_smoke  = pctl(es_smoke),
  # ---- external / validation (never model inputs) ----
  clock_age                 = clock_age,
  val_inflammation_episcore = meanp(es_crp, es_il6),
  val_gdf15_episcore        = pctl(es_gdf15),
  inflammation_def          = if (infl_crp) "NLR+CRP" else "NLR"
)
fwrite(feat, outfile)
na_cols <- names(feat)[sapply(feat, function(c) all(is.na(c)))]
message(sprintf("\nDONE -> %s  (%d samples)", outfile, nrow(feat)))
if (length(na_cols)) message("  all-NA columns (weight table not supplied; imputed downstream): ",
                             paste(na_cols, collapse = ", "))
message("  Next:  python3 geo_validation.py --features ", outfile)

# -----------------------------------------------------------------------------
# features.csv columns (+ inflammation_def: "NLR" or "NLR+CRP", provenance for the inflammation check):
#   sample_id, age, sex,
#   inflammation, immune_aging, nad_depletion, hepatic, adiposity_lep,
#   adiponectin, damage_smoke        <- the 7 model inputs (percentiles 0-1)
#   clock_age                        <- EXTERNAL clock (years), comparison only
#   val_inflammation_episcore, val_gdf15_episcore  <- validation percentiles only
# -----------------------------------------------------------------------------
