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
# option A — venv
python3 -m venv plmap && source plmap/bin/activate     # Windows: plmap\Scripts\activate
pip install -r requirements.txt

# option B — conda
conda env create -f environment.yml && conda activate plmap
```

Check it works (no data needed):

```bash
cd src
python3 coupled_model.py        # -> "98 nodes; 22 inputs", writes MUS_P1_coupled.bnd/.cfg
python3 methylation_layer.py    # 4 synthetic patients -> profiles + recommendations
python3 geo_validation.py       # calibrated cohort -> aging index vs age (r ≈ 0.73) + figure
```

If `import cmaboss_128n` fails, your `cmaboss` wheel is too old — `pip install -U cmaboss`
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
`validation/geo_validation_pipeline.md`; this is the operational recipe. **Both steps are now
implemented** — `src/preprocess_gse40279.R` (heavy, R) and `src/geo_validation.py --features` (light,
Python).

### 4a. R environment for preprocessing (once)

```r
install.packages(c("BiocManager", "data.table"))
BiocManager::install(c("GEOquery", "EpiDISH", "methylclock"))   # minfi/sesame optional (raw IDATs)
```

### 4b. EpiScore weight tables (once, user-supplied)

Cell composition and the clock run out of the box. The protein/exposure **EpiScores need published CpG
weight tables** — drop one CSV per score into `episcore_weights/` (columns `cpg,weight`, optional
`(Intercept)` row):

```
episcore_weights/CRP.csv  IL6.csv  GDF15.csv  ADIPOQ.csv  LEP.csv  liver.csv  smoking.csv  [NAD_CD38.csv]
```

Sources: Gadd et al. 2022 *eLife* (protein EpiScores), McCartney et al. 2018 (smoking), or the
MethylDetectR bundle (Hillary & Marioni). Any missing table → that feature is `NA` and imputed to the
cohort median downstream (the Python step reports which), so a **partial run still works**.

### 4c. Produce the feature matrix (heavy, R)

```bash
cd src
Rscript preprocess_gse40279.R --outdir . --weights episcore_weights
#  -> features.csv : one row per sample; 7 model features as cohort percentiles (0–1)
#     + age, sex, clock_age (external), val_inflammation_episcore, val_gdf15_episcore
```

The script downloads GSE40279 (GEOquery), runs EpiDISH cell composition, projects the EpiScore weight
tables onto the beta matrix, computes the Horvath/Hannum clock (external only), and writes the
percentile feature matrix. ~16 GB RAM for the 450K matrix.

### 4d. Run the cohort through the model (light, Python — runs on your PC)

```bash
python3 geo_validation.py --features features.csv
#  -> geo_validation_real.json + fig_geo_validation_real.png
```

This reads `features.csv`, imputes any absent feature to the cohort median (reported), runs every
sample through `methylation_layer.run_patient`, and writes the correlations + figure. Without
`--features` the same script runs the literature-calibrated simulated cohort (the pipeline proof).

### 4e. Read out (printed + in the JSON)

- `Aging_index` + 7 sub-scores vs chronological age (Pearson + Spearman + slope/decade).
- `Aging_index` vs the clock (external), and age-acceleration: residual(clock~age) vs residual(index~age).
- **Non-circular** checks: `Aging_inflammation` vs the CRP/IL-6 EpiScore; `Aging_secretory` vs the GDF15
  EpiScore.
- Note the coverage boundary (see the pipeline doc): the blood-feature set drives the
  inflammation/NAD/nutrient/mito/secretory arms, not genomic-senescence/proteostasis/Sarcopenia — those
  need additional inputs.

### Honest split
4c is the only heavy, R-based part. 4d–4e reuse this repo unchanged. Run 4c where the data and R live;
4d–4e run anywhere Python does — including your PC.

---

## 5. Notes for moving to GitHub
- `src/` is self-contained and path-free (run scripts from inside `src/`).
- Never commit raw methylation matrices (see `.gitignore`); keep `data/` for small derived JSONs only.
- The two-environment split (Python model / R preprocessing) is deliberate — keep them separate so the
  model stays laptop-friendly.
