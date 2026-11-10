# Setup & reproduction — mirror of the working environment

This documents **exactly what I used**, how to reproduce it on your own PC, and the step-by-step for
the real GSE40279 run. Two environments: a light **Python model environment** (everything in this
repo) and a heavier **R methylation environment** (only for turning a real methylation matrix into
model features).

---

## 1. What I used (mirror this)

| Component | Version | Role |
|---|---|---|
| OS | Linux x86-64 | any OS works (Win/Mac/Linux) |
| Python | **3.13.17** | model + analysis |
| cmaboss | **1.0.0b34** | MaBoSS engine; ships `cmaboss_128n` (≤128 nodes), `_256n`, `_512n`, `_1024n` |
| numpy | 2.5.3 | arrays |
| matplotlib | 3.11.2 | figures |
| openpyxl | 3.1.5 | Excel workbooks |
| playwright | 1.56 (optional) | only to re-render the schematic figures |

The coupled model imports **`cmaboss_128n`** (not plain `cmaboss`) because it has 98 nodes and the
default `cmaboss` build caps at 64. `pip install cmaboss` ships all the larger builds in one package.

## 2. Python environment (runs the whole model on your PC)

```bash
# option A — uv (recommended; fast, and fetches Python 3.13 itself if missing)
uv venv --python 3.13                  # creates .venv/ in the repo root (git-ignored)
source .venv/bin/activate              # Windows: .venv\Scripts\activate
uv pip install -r requirements.txt

# option B — plain venv
python3 -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install -r requirements.txt

# option C — conda
conda env create -f environment.yml && conda activate plmap
```

With uv you can also skip activating: from `src/`, run e.g.
`uv run --python ../.venv python coupled_model.py`.

Check it works (no data needed):

```bash
cd src
python3 coupled_model.py        # -> "98 nodes; 22 inputs", writes MUS_P1_coupled.bnd/.cfg
python3 methylation_layer.py    # 4 synthetic patients -> profiles + recommendations
python3 geo_validation.py       # calibrated cohort -> aging index vs age (r ≈ 0.73) + figure
```

If `import cmaboss_128n` fails, your `cmaboss` wheel is too old — `uv pip install -U cmaboss` (or `pip install -U cmaboss`)
(need ≥ 1.0.0b34). Nothing here needs a GPU or large RAM; a laptop is fine.

## 3. What runs where

- **On your personal PC, immediately:** the whole model — building it, scenarios, robustness,
  characterization, the aging layer, the methylation layer on synthetic or already-derived features,
  and the calibrated-cohort validation. Seconds to ~2 min each.
- **Needs the R methylation stack (one dataset, ~16 GB RAM):** turning a *raw* GEO methylation matrix
  into the features the model consumes (cell composition, EpiScores, a clock). Do this once per
  dataset; then the per-sample model runs are light again.

---

## 4. Real GSE40279 run — step by step

GSE40279 = Hannum 2013, 656 whole-blood 450K samples, ages 19–101 (PMID 23177740). Full rationale,
the feature mapping, the `features.csv` schema and the coverage caveats are in
`validation/geo_validation_pipeline.md`; this is the operational recipe. Three scripts, run in order:

| step | script | what | cost |
|---|---|---|---|
| download | `src/download_gse40279.R` | fetch GSE40279 **once**, cache as `.rds` + a small age-stratified test subset | network, ~16 GB RAM, once |
| process | `src/preprocess_gse40279.R [--test]` | cache → `features.csv` (`--test`: subset → `features_test.csv`) | test: minutes / few GB; full: ~16 GB |
| model | `src/geo_validation.py --features …` | features → correlations + figure | light Python, laptop |

Download and processing are separate so the big download happens once and processing can be re-run
offline as often as needed (new weight tables, tweaks). The test subset runs the **same** processing
code, so a passing `--test` run means the full run will work too.

### 4a. R environment for preprocessing (once)

```r
install.packages(c("BiocManager", "data.table"))
BiocManager::install(c("GEOquery", "Biobase", "EpiDISH", "methylclock"))   # minfi/sesame optional (raw IDATs)
# download step needs GEOquery + Biobase; processing step needs EpiDISH + methylclock
```

### 4b. EpiScore weight tables (once) — populate `src/episcore_weights/`

Cell composition and the clock run out of the box. The protein/exposure **EpiScores need published CpG
<<<<<<< HEAD
weight tables** — drop one CSV per score into `episcore_weights/` at the repo root (columns `cpg,weight`, optional
`(Intercept)` row):
=======
weight tables**, one CSV per score (`cpg,weight`, optional `(Intercept)`) in **`src/episcore_weights/`**
(that path — not the repo root — because `preprocess_gse40279.R` is run from `src/`). The full recipe,
mapping, and sources are in `src/episcore_weights/README.md`. The fast path:
>>>>>>> 9ae43ce (Add EpiScore weights converter + fix folder path convention)

```bash
cd src
curl -L -o Predictors_Shiny_by_Groups.csv \
  "https://zenodo.org/records/15548022/files/Predictors_Shiny_by_Groups.csv?download=1"
# in R once:  writeLines(rownames(beta), "beta_cpgs.txt")
python3 make_episcore_weights.py Predictors_Shiny_by_Groups.csv \
    --outdir episcore_weights --beta-cpgs beta_cpgs.txt \
    --map CRP.csv="CRP" ADIPOQ.csv="Adiponectin" LEP.csv="Body Fat %" liver.csv="HGF" smoking.csv="Smoking"
```

That writes `CRP/ADIPOQ/LEP/liver/smoking.csv` (→ inflammation, adiponectin, adiposity_lep, hepatic,
damage_smoke). `GDF15`/`IL6`/true leptin aren't in the bundle (GDF15 is validation-only; leptin is
covered by Body Fat %). Any missing table → that feature is `NA` and imputed downstream, so a partial
run still works.

### 4c. Download once (R)

```bash
cd src
Rscript download_gse40279.R              # -> ../data/raw/gse40279/{beta,pheno}.rds + *_test.rds
```

Caches the series matrix and the parsed phenotype (age, sex) under `data/raw/gse40279/` (git-ignored),
plus a 48-sample subset spread evenly over the age range (`--test-n` to change). Re-running skips the
download if the cache exists (`--force` to refetch).

### 4d. Produce the feature matrix (R) — test first, then full

```bash
Rscript preprocess_gse40279.R --test     # subset -> features_test.csv  (minutes, a few GB RAM)
Rscript preprocess_gse40279.R            # full   -> features.csv       (~16 GB RAM)
#  one row per sample; 7 model features as cohort percentiles (0–1)
#  + age, sex, clock_age (external), val_inflammation_episcore, val_gdf15_episcore
```

Runs EpiDISH cell composition, projects the EpiScore weight tables (default `../episcore_weights`)
onto the beta matrix, computes the Horvath/Hannum clock (external only), and writes the percentile
feature matrix. No network access. The `--test` run checks the R stack, the weight tables and their
CpG overlaps; its percentiles are relative to the subset, so treat its numbers as a check, not a result.

### 4e. Run the cohort through the model (light, Python — runs on your PC)

```bash
python3 geo_validation.py --features features_test.csv   # -> geo_validation_real_test.json + fig_…_test.png
python3 geo_validation.py --features features.csv        # -> geo_validation_real.json + fig_geo_validation_real.png
```

This reads `features.csv`, imputes any absent feature to the cohort median (reported), runs every
sample through `methylation_layer.run_patient`, and writes the correlations + figure. Without
`--features` the same script runs the literature-calibrated simulated cohort (the pipeline proof).

### 4f. Read out (printed + in the JSON)

- `Aging_index` + 7 sub-scores vs chronological age (Pearson + Spearman + slope/decade).
- `Aging_index` vs the clock (external), and age-acceleration: residual(clock~age) vs residual(index~age).
- **Non-circular** checks: `Aging_inflammation` vs the CRP/IL-6 EpiScore; `Aging_secretory` vs the GDF15
  EpiScore.
- Note the coverage boundary (see the pipeline doc): the blood-feature set drives the
  inflammation/NAD/nutrient/mito/secretory arms, not genomic-senescence/proteostasis/Sarcopenia — those
  need additional inputs.

### Honest split
4c–4d are the only heavy, R-based part. 4e–4f reuse this repo unchanged. Run 4c–4d where the data and
R live; 4e–4f run anywhere Python does — including your PC.

---

## 5. Notes for moving to GitHub
- `src/` is self-contained and path-free (run scripts from inside `src/`).
- Never commit raw methylation matrices (see `.gitignore`); keep `data/` for small derived JSONs only.
- The two-environment split (Python model / R preprocessing) is deliberate — keep them separate so the
  model stays laptop-friendly.
