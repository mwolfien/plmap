# GEO validation — real-data pipeline (GSE40279)

**Goal.** Test whether the coupled model's aging index / hallmark sub-scores, parameterized from
*real* blood-methylation features, track chronological age and an epigenetic clock.

**Why GSE40279** (Hannum 2013, PMID 23177740, [doi:10.1016/j.molcel.2012.10.016](https://doi.org/10.1016/j.molcel.2012.10.016)):
656 whole-blood Illumina 450K samples, ages 19–101, public, widely used for clocks — the cleanest
public match for a blood-readout aging model. Alternatives: GSE87571 (wide age range), GSE55763
(large), GSE132203 (with phenotypes).

**Status.** The machinery is proven on a literature-calibrated simulated cohort
(`python3 geo_validation.py`, `fig_geo_validation.png`): the aging index tracks simulated age
(r ≈ 0.73), NAD-decline and inflammaging sub-scores strongest. The real-data path is now **wired**:
`download_gse40279.R` fetches and caches the data once, `preprocess_gse40279.R` produces the feature
matrix (`--test` on a small age-stratified subset first), and `geo_validation.py --features features.csv`
consumes it. Preprocessing needs the beta matrix + an R methylation stack + ~16 GB RAM, so run it
where the data lives; everything after `features.csv` is light Python.

## Commands (from `src/`)

```bash
# 0. ONCE, R, needs network: GEO -> cached ../data/raw/gse40279/{beta,pheno}.rds (+ *_test.rds subset)
Rscript download_gse40279.R

# 1. HEAVY, R, offline: cached matrix -> features.csv   (--test: 48-sample subset -> features_test.csv)
Rscript preprocess_gse40279.R --test
Rscript preprocess_gse40279.R

# 2. LIGHT, Python, anywhere (laptop): features.csv -> correlations + figure
python3 geo_validation.py --features features.csv
```

## Step 0 — `download_gse40279.R` (once)

Downloads the GSE40279 series matrix with GEOquery, parses age and sex from the GEO characteristics
columns, and caches the beta matrix + phenotype as `.rds` under `data/raw/gse40279/` (git-ignored). It
also writes a small subset (default 48 samples, evenly spaced over the age range) for the `--test` run.
Re-running reuses the cache; `--force` refetches.

## Step 1 — `preprocess_gse40279.R` (the heavy R job)

Reads the cache (no network), then derives, per sample:

- **Cell composition** — `EpiDISH` RPC with the built-in `centDHSbloodDMC.m` blood reference
  (CD4T, CD8T, NK, B, Mono, Neu). Derived: neutrophil:lymphocyte ratio → `inflammation`;
  CD4:CD8 (falls with age) → `immune_aging`. *Runs out of the box.*
- **Protein / exposure EpiScores** — a projection of published CpG weight tables onto beta:
  CRP, IL-6, GDF15, adiponectin, leptin, liver, smoking. **You supply the weight tables** as
  `episcore_weights/<NAME>.csv` (columns `cpg,weight`, optional `(Intercept)` row). Sources:
  Gadd et al. 2022 *eLife* (protein EpiScores), McCartney et al. 2018 (smoking), or the
  MethylDetectR bundle (Hillary & Marioni). Any table that is absent → that feature is written
  `NA` and imputed to the cohort median downstream (reported), so a **partial run still works**.
- **Epigenetic clock** — `methylclock::DNAmAge` (Horvath/Hannum). **External readout only — never a
  model input** (this is what keeps the validation non-circular).

Each model feature is written as a **cohort percentile (0–1)** — exactly what
`methylation_layer.run_patient` expects.

**`features.csv` schema** (one row per sample):

| column | kind | → model |
|---|---|---|
| `sample_id`, `age`, `sex` | metadata | — |
| `inflammation`, `immune_aging`, `nad_depletion`, `hepatic`, `adiposity_lep`, `adiponectin`, `damage_smoke` | 7 model inputs (percentile 0–1) | `Inflammaging`, `LEP`/`Nutrients`, `ADIPOQ`, `DNA_damage`, `hyperammonemia`/`FGL1`, … |
| `clock_age` | external (years) | comparison only |
| `val_inflammation_episcore`, `val_gdf15_episcore` | validation (percentile) | non-circular checks |

> `nad_depletion` has no standard blood EpiScore; supply a CD38/NAD surrogate table as
> `episcore_weights/NAD_CD38.csv` or it is imputed to the cohort median.

## Step 2 — `geo_validation.py --features features.csv` (light Python)

Runs every sample through the coupled model and reports:

- `Aging_index` + 7 sub-scores vs **chronological age** (Pearson + Spearman + slope/decade).
- `Aging_index` vs the **clock** (external, non-circular), plus an **age-acceleration** test:
  residual(clock~age) vs residual(index~age) — do people the clock calls "old for age" also score
  high mechanistically?
- **Non-circular molecular checks**: `Aging_inflammation` vs the CRP/IL-6 EpiScore;
  `Aging_secretory` vs the GDF15 EpiScore.
- Writes `geo_validation_real.json` + `fig_geo_validation_real.png`.

## What the blood-methylation feature set can and cannot drive (important)

The wiring test (synthetic cohort in the real schema) confirms the layer strongly parameterizes the
**inflammation / NAD / nutrient-sensing / mitochondrial / secretory** arms — those sub-scores rise
with age as designed. It does **not**, on its own, drive **genomic-senescence**, **proteostasis**, or
the **muscle/Sarcopenia** arm, because the 7 blood features don't supply those inputs. Expect those
sub-scores to be flat vs age until the layer is extended with:

- a damage/senescence EpiScore (e.g. an EpiTOC2 mitotic-age or a p16/senescence surrogate) → genomic-senescence;
- a proteostasis-relevant signal → proteostasis;
- muscle-relevant methylation features or clinical/DXA inputs → Sarcopenia (the muscle module was
  designed to be parameterized from more than blood).

This is an honest coverage boundary, not a failure — it tells the student exactly which inputs to add
next, and the correlations on the covered arms are the meaningful, non-circular validation.

## Feasibility note
Steps 0–1 (download, then EpiDISH + EpiScores + clock) are a standard but heavy R job (GEOquery/EpiDISH/
methylclock + reference panels + a large matrix, ~16 GB RAM; the `--test` run needs only a few GB). Step 2 reuses this repo unchanged and
runs anywhere Python does. Run step 1 where the data + R live; step 2 on a laptop.
