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

GSE40279 = Hannum 2013, 656 whole-blood 450K samples, ages 19–101 (PMID 23177740). Full rationale and
the mapping are in `validation/geo_validation_pipeline.md`; this is the operational recipe.

### 4a. R environment for preprocessing (once)

```r
install.packages("BiocManager")
BiocManager::install(c("minfi", "EpiDISH", "GEOquery", "sesame"))
# EpiScore weights: MethylDetectR (Hillary & Marioni) or the published EpiScore CpG tables
# Clocks: methylclock (Bioconductor) or dnaMethyAge
BiocManager::install(c("methylclock"))
```

### 4b. Get the data

```r
library(GEOquery)
gse <- getGEO("GSE40279", GSEMatrix = TRUE)[[1]]
beta <- exprs(gse)                      # CpGs × 656 samples (beta values)
age  <- as.numeric(gsub("age \\(y\\): ", "", pData(gse)$characteristics_ch1.1))
# (sex, etc. are in pData too)
```

### 4c. Derive the model features (per sample, all standard)

```r
library(EpiDISH)
cell  <- epidish(beta, ref.m = centDHSbloodDMC.m, method = "RPC")$estF   # CD4T,CD8T,NK,B,Mono,Neu
# features: neutrophil:lymphocyte, CD4:CD8  -> inflammation, immune_aging
# EpiScores (CRP, IL-6, GDF15, adiponectin, leptin, liver enzymes, smoking): apply the
#   published CpG weight tables to beta -> inflammation, hepatic, adiposity_lep, adiponectin,
#   damage_smoke, and a GDF15 EpiScore to *validate* the secretory sub-score
library(methylclock)
clock <- DNAmAge(beta)                  # Horvath / Hannum -> EXTERNAL comparison only, never an input
```

Write one CSV, `features.csv`: one row per sample, columns =
`age, inflammation, immune_aging, nad_depletion, hepatic, adiposity_lep, adiponectin, damage_smoke, clock_age`.
Each feature is its **cohort percentile (0–1)** — that is what the model layer expects.

### 4d. Run the cohort through the model (Python, light)

Point the existing `src/geo_validation.py` at `features.csv` instead of the simulated cohort: replace
the synthetic `feat` dict with each row's columns (already 0–1), then it runs each sample through
`methylation_layer.run_patient`, computes the aging index + sub-scores, and correlates them with `age`
and `clock_age`. (A ~15-line edit — the `for age in ages:` loop becomes `for row in features.csv`.)

### 4e. Read out

- `Aging_index` vs chronological age and vs the clock (Pearson/Spearman + scatter).
- Sub-score-specific, **non-circular** checks: does `Aging_NAD` track age, `Aging_inflammation` track
  the CRP/IL-6 EpiScore, `Aging_secretory` track the **GDF15 EpiScore**?
- Age-acceleration: residual(clock~age) vs residual(index~age) — do people the clock calls "old for
  age" score high on specific mechanistic hallmarks?

### Honest split
4a–4c (download + EpiDISH + EpiScores + clock) are the only heavy, R-based part. 4d–4e reuse this repo
unchanged. Run 4a–4c where the data and R live; 4d–4e run anywhere Python does — including your PC.

---

## 5. Notes for moving to GitHub
- `src/` is self-contained and path-free (run scripts from inside `src/`).
- Never commit raw methylation matrices (see `.gitignore`); keep `data/` for small derived JSONs only.
- The two-environment split (Python model / R preprocessing) is deliberate — keep them separate so the
  model stays laptop-friendly.
