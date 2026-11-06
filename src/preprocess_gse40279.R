#!/usr/bin/env Rscript
# =============================================================================
# preprocess_gse40279.R
# Raw GEO 450K methylation  ->  model feature matrix (features.csv)
#
# Turns the GSE40279 whole-blood 450K beta matrix (Hannum 2013, PMID 23177740,
# 656 samples, ages 19-101) into the exact per-sample features the Prevention &
# Longevity Map methylation layer consumes:  7 model-input features + age +
# an EXTERNAL clock age + two validation-only EpiScores. Each model feature is a
# COHORT PERCENTILE (0-1), which is what src/methylation_layer.run_patient expects.
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
#       CpG weight tables dropped into  episcore_weights/  (see section 3).
#   Any EpiScore whose weight table is absent is written as NA and imputed to the
#   cohort median downstream (geo_validation.py reports which were imputed), so a
#   PARTIAL run (cell composition + clock + whatever scores you have) still works.
#
# Heavy step: this is the only part of the project that needs R/Bioconductor,
# the raw matrix (~GB) and ~16 GB RAM. Run it where the data lives; everything
# after features.csv is light Python that runs on a laptop.
#
# Usage:
#   Rscript preprocess_gse40279.R [--outdir .] [--weights episcore_weights]
# Output:
#   features.csv   (one row per sample; columns documented at the bottom)
# =============================================================================

suppressMessages({
  need <- c("GEOquery", "EpiDISH", "methylclock", "data.table")
  miss <- need[!sapply(need, requireNamespace, quietly = TRUE)]
  if (length(miss)) stop(
    "Missing packages: ", paste(miss, collapse = ", "),
    "\nInstall with:\n  install.packages('BiocManager')\n",
    "  BiocManager::install(c('GEOquery','EpiDISH','methylclock','data.table'))")
  library(GEOquery); library(EpiDISH); library(methylclock); library(data.table)
})

# ---- args -------------------------------------------------------------------
args    <- commandArgs(trailingOnly = TRUE)
getarg  <- function(flag, default) { i <- which(args == flag); if (length(i)) args[i + 1] else default }
outdir  <- getarg("--outdir", ".")
wdir    <- getarg("--weights", "episcore_weights")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
pctl    <- function(x) { r <- rank(x, na.last = "keep", ties.method = "average"); (r - 1) / (sum(!is.na(x)) - 1) }

# =============================================================================
# 1. Download + assemble beta matrix and phenotype
# =============================================================================
message("[1/5] Downloading GSE40279 (this is the big step)...")
gse  <- getGEO("GSE40279", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
beta <- as.matrix(exprs(gse))                      # CpGs x 656 (beta values 0-1)
ph   <- pData(gse)

# age + sex live in the characteristics columns; parse defensively.
grab <- function(ph, key) {
  col <- names(ph)[sapply(ph, function(c) any(grepl(key, c, ignore.case = TRUE)))][1]
  sub(paste0(".*", key, "[: ]*"), "", grep(key, ph[[col]], ignore.case = TRUE, value = FALSE))
}
age  <- as.numeric(sub(".*age \\(y\\): ",   "", ph[[grep("age",    sapply(ph, function(c) any(grepl("age \\(y\\)", c))))[1]]]))
if (all(is.na(age))) {                             # fallback parser
  ac  <- ph[[ names(ph)[grepl("characteristics", names(ph))][
            which(sapply(ph[grepl("characteristics", names(ph))], function(c) any(grepl("age", c, TRUE))))[1]] ]]
  age <- as.numeric(gsub("[^0-9.]", "", ac))
}
sexc <- ph[[ names(ph)[grepl("characteristics", names(ph))][
          which(sapply(ph[grepl("characteristics", names(ph))], function(c) any(grepl("gender|sex", c, TRUE))))[1]] ]]
sex  <- ifelse(grepl("F", sexc, ignore.case = TRUE), "F", "M")
message(sprintf("    beta: %d CpGs x %d samples | age %.0f-%.0f", nrow(beta), ncol(beta), min(age, na.rm=T), max(age, na.rm=T)))

# =============================================================================
# 2. Cell composition (EpiDISH, built-in blood reference) -> inflammation/immune
# =============================================================================
message("[2/5] Estimating leukocyte composition (EpiDISH RPC)...")
data(centDHSbloodDMC.m)                             # 6-cell blood reference (ships with EpiDISH)
cell <- epidish(beta.m = beta, ref.m = centDHSbloodDMC.m, method = "RPC")$estF  # CD4T,CD8T,NK,B,Mono,Neu
cell <- as.data.frame(cell)
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
es_crp   <- episcore(beta, "CRP.csv")
es_il6   <- episcore(beta, "IL6.csv")
es_gdf15 <- episcore(beta, "GDF15.csv")
es_adipoq<- episcore(beta, "ADIPOQ.csv")
es_leptin<- episcore(beta, "LEP.csv")
es_liver <- episcore(beta, "liver.csv")            # ALT/GGT composite, or a single hepatic EpiScore
es_smoke <- episcore(beta, "smoking.csv")

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
#      inflammation  = mean percentile of CRP, IL6, neutrophil:lymphocyte   (higher = more inflamed)
#      immune_aging  = percentile of LOW CD4:CD8 (ratio falls with age)     (higher = more aged)
#      nad_depletion = NA unless a CD38/NAD surrogate score is supplied      (imputed median downstream)
#      hepatic       = percentile of the liver EpiScore
#      adiposity_lep = percentile of the leptin EpiScore
#      adiponectin   = percentile of the adiponectin EpiScore (protective; model treats as ADIPOQ input)
#      damage_smoke  = percentile of the smoking EpiScore
#    Validation-only (NOT model inputs):
#      val_inflammation_episcore = percentile of mean(CRP,IL6)  -> checks Aging_inflammation
#      val_gdf15_episcore        = percentile of GDF15          -> checks Aging_secretory
# =============================================================================
message("[5/5] Building percentile features + writing features.csv...")
meanp <- function(...) { m <- cbind(...); rowMeans(apply(m, 2, pctl), na.rm = TRUE) }

feat <- data.table(
  sample_id     = colnames(beta),
  age           = age,
  sex           = sex,
  # ---- 7 model-input features (cohort percentiles 0-1) ----
  inflammation  = meanp(es_crp, es_il6, nlr),
  immune_aging  = pctl(-cd48),
  nad_depletion = episcore(beta, "NAD_CD38.csv") |> pctl(),   # NA unless supplied
  hepatic       = pctl(es_liver),
  adiposity_lep = pctl(es_leptin),
  adiponectin   = pctl(es_adipoq),
  damage_smoke  = pctl(es_smoke),
  # ---- external / validation (never model inputs) ----
  clock_age                 = clock_age,
  val_inflammation_episcore = meanp(es_crp, es_il6),
  val_gdf15_episcore        = pctl(es_gdf15)
)
fwrite(feat, file.path(outdir, "features.csv"))
na_cols <- names(feat)[sapply(feat, function(c) all(is.na(c)))]
message(sprintf("\nDONE -> %s  (%d samples)", file.path(outdir, "features.csv"), nrow(feat)))
if (length(na_cols)) message("  all-NA columns (weight table not supplied; imputed downstream): ",
                             paste(na_cols, collapse = ", "))
message("  Next:  python3 geo_validation.py --features ", file.path(outdir, "features.csv"))

# -----------------------------------------------------------------------------
# features.csv columns:
#   sample_id, age, sex,
#   inflammation, immune_aging, nad_depletion, hepatic, adiposity_lep,
#   adiponectin, damage_smoke        <- the 7 model inputs (percentiles 0-1)
#   clock_age                        <- EXTERNAL clock (years), comparison only
#   val_inflammation_episcore, val_gdf15_episcore  <- validation percentiles only
# -----------------------------------------------------------------------------
