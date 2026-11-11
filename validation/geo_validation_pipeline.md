# GEO validation — real-data pipeline (GSE40279)

**Goal.** Test whether the coupled model's aging index / hallmark sub-scores, parameterized from
*real* blood-methylation features, track chronological age and an epigenetic clock.

**Why GSE40279** (Hannum 2013, PMID 23177740, [doi:10.1016/j.molcel.2012.10.016](https://doi.org/10.1016/j.molcel.2012.10.016)):
656 whole-blood Illumina 450K samples, ages 19–101, public, widely used for clocks — the cleanest
public match for a blood-readout aging model. Alternatives: GSE87571 (wide age range), GSE55763
(large), GSE132203 (with phenotypes).

**Status.** The machinery is proven on a literature-calibrated simulated cohort
(`python3 geo_validation.py`, `fig_geo_validation.png`): the aging index tracks simulated age
(r ≈ 0.63; NAD 0.80, inflammation 0.77 — after NAD became emergent), NAD-decline and inflammaging sub-scores strongest. The real-data path is now **wired**:
`download_geo.R` fetches and caches the data once, `preprocess_geo.R` produces the feature
matrix (`--test` on a small age-stratified subset first), and `geo_validation.py --features features_<gse>.csv`
consumes it. Two cohorts are kept strictly apart: **GSE40279 is the development cohort** (all tuning)
and **GSE87571 the held-out validation cohort** (run once, with the settings frozen on GSE40279). Preprocessing needs the beta matrix + an R methylation stack + ~16 GB RAM, so run it
where the data lives; everything after the features file is light Python.

## Commands (from `src/`)

```bash
# 0. ONCE per cohort, R, needs network: GEO -> cached ../data/raw/<gse>/{beta,pheno}.rds (+ *_test.rds)
Rscript download_geo.R --gse40279
Rscript download_geo.R --gse87571

# 1. HEAVY, R, offline: cached matrix -> features_<gse>.csv   (--test: subset -> features_<gse>_test.csv)
Rscript preprocess_geo.R --gse40279 --test
Rscript preprocess_geo.R --gse40279 --by-sex
Rscript preprocess_geo.R --gse87571 --min-age 19 --by-sex      # validation: same options as dev

# 2. LIGHT, Python, anywhere (laptop): features -> correlations + figure (+ per-sex table)
python3 geo_validation.py --features features_gse40279.csv --by-sex
python3 geo_validation.py --features features_gse87571.csv --by-sex   # once, frozen settings
```

## Development vs validation, and sex

- **Tune only on GSE40279**: feature subset (`sweep_configs.py`), index aggregation, pooled vs
  within-sex percentiles, clock. `sweep_configs.py` refuses to run on GSE87571.
- **Then run GSE87571 once** with exactly those settings and report every pre-specified metric.
- **Sex:** body fat, CRP and leukocyte composition differ by sex, so pooled percentiles partly rank
  sex rather than biology. `--by-sex` ranks each sample against same-sex peers, and
  `geo_validation.py --by-sex` reports the results for women and men separately. Compare pooled vs
  within-sex on GSE40279 and freeze the choice before touching GSE87571. Separate model *weights*
  per sex would be a larger step; within-sex percentiles plus per-sex reporting come first.

## Step 0 — `download_geo.R` (once)

Downloads the named cohort's series matrix with GEOquery, parses age and sex from the GEO characteristics
columns, and caches the beta matrix + phenotype as `.rds` under `data/raw/gse40279/` (git-ignored). It
also writes a small subset (default 48 samples, evenly spaced over the age range) for the `--test` run.
Re-running reuses the cache; `--force` refetches.

## Step 1 — `preprocess_geo.R` (the heavy R job)

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

**`features_<gse>.csv` schema** (one row per sample):

| column | kind | → model |
|---|---|---|
| `sample_id`, `age`, `sex` | metadata | — |
| `inflammation`, `immune_aging`, `nad_depletion`, `hepatic`, `adiposity_lep`, `adiponectin`, `damage_smoke` | 7 model inputs (percentile 0–1) | `Inflammaging`, `LEP`/`Nutrients`, `ADIPOQ`, `DNA_damage`, `hyperammonemia`/`FGL1`, … |
| `clock_age` (+ `clock_name`, `clock_hannum`, `clock_horvath`) | external (years) | comparison only |
| `val_inflammation_episcore`, `val_gdf15_episcore` | validation (percentile) | sub-score checks |
| `inflammation_def`, `cohort`, `percentile_scope` | provenance | how the file was built |

> `nad_depletion` has no standard blood EpiScore; supply a CD38/NAD surrogate table as
> `episcore_weights/NAD_CD38.csv` or it is imputed to the cohort median.

## Step 2 — `geo_validation.py --features features_<gse>.csv` (light Python)

Runs every sample through the coupled model and reports:

- `Aging_index` + 7 sub-scores vs **chronological age** (Pearson + Spearman + slope/decade).
- `Aging_index` vs the **clock** (external, non-circular), plus an **age-acceleration** test:
  residual(clock~age) vs residual(index~age) — do people the clock calls "old for age" also score
  high mechanistically?
- **Molecular checks**: `Aging_inflammation` vs the CRP EpiScore (independent: since the N:L-only change
  CRP is no longer part of the `inflammation` input; `inflammation_def` in the features file records which
  definition was used), and `Aging_secretory` vs the GDF15 EpiScore (needs a GDF15 table).
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

## Troubleshooting weak correlations

Run the preflight gate first — it checks the feature matrix in seconds and exits non-zero if inputs are
empty, so you catch it before the (slow) model run:

```bash
python3 check_features.py features_gse40279.csv # per-feature coverage + r vs age; FAIL on empty inputs
```

If `Aging_index` vs age comes back near zero (r ≈ 0.1), check the feature matrix **before** the model:

1. **How many features actually have data?** `geo_validation.py` now prints an `EXCLUDED` line and a
   `WARNING` when a column is entirely empty. Empty EpiScore columns mean `episcore_weights/` had no
   table, or the table's CpG IDs did not overlap the beta-matrix rownames (the R script writes `NA` on
   zero overlap). This is the usual cause — the model then runs on cell-composition features only, whose
   age ceiling is ~0.18.
2. **Do the present features track age at all?** Correlate each feature column against `age` directly. In
   GSE40279, `inflammation` (neutrophil:lymphocyte) tracks age (r ≈ 0.19) but `immune_aging` (CD4:CD8)
   does not (r ≈ −0.07) — a poorly-estimated feature dilutes the signal, so validate per feature before
   trusting the index.
3. **Missing columns are left OFF, not imputed to a constant.** A fully-empty column is excluded so its
   inputs stay off; combined inputs use a weighted blend (not a capped sum) so shared hubs like
   Inflammaging cannot saturate. These fixes restore honest dynamic range but cannot create signal that
   is absent from the inputs — the lift comes from supplying the EpiScore tables.

## Findings (GSE40279 real run, n=656)

A systematic sweep (`src/sweep_configs.py`: feature subsets x index aggregations) over the real cohort:

- **OLS ceiling ≈ 0.40** (in-sample) — the maximum age signal *linearly extractable* from the blood features. This is
  the number to beat; a mechanistic model that merely re-predicts chronological age adds nothing over an
  epigenetic clock (which hits r≈0.95; here `clock_age` vs age = 0.946).
- **Default setup (all features + the model's built-in `Aging_index` node): r ≈ 0.00 vs age** — the signal
  is cancelled by (a) wrong-signed/flat features and (b) an unstable aggregation node.
  - Feature-level truth: `inflammation` +0.31 and `adiponectin` +0.18 track age; `adiposity_lep` is
    **−0.20** (the Body Fat % leptin-substitute falls with age here), `immune_aging` −0.07 and
    `damage_smoke` −0.06 are flat. The bad ones drag the index to zero.
- **Two fixes recover it to r ≈ 0.24 (≈60% of the ceiling)** — *in-sample*: the subset was chosen because its
  features track age in this same cohort, so 0.24 is optimistic until confirmed on held-out samples or a
  second dataset (e.g. GSE87571):
  1. **Feature selection** — use `--use-features inflammation,adiponectin,hepatic` ("core3"); drop the
     wrong-signed/flat features.
  2. **Aggregation** — use `Aging_index_mean` (mean of the 7 sub-scores), which is robust; the model's
     built-in node gave +0.24 for core3 but −0.02 for core2, while the mean gave ~0.23 for both.
- **`Aging_inflammation` vs the CRP/IL6 EpiScore = 0.86–0.88 — a transduction check, NOT independent
  validation.** The CRP EpiScore is part of the `inflammation` input (`preprocess_geo.R`:
  `inflammation = meanp(CRP, IL6, NLR)`, `val_inflammation_episcore = meanp(CRP, IL6)`; with no IL6 table
  the check column is just CRP, half of the input). With pure-noise CRP and NLR the same check already
  gives r ≈ 0.6, and it rises as other drivers of Inflammaging are dropped, which is why it "strengthens"
  for core3. It shows the input reaches the sub-score intact, nothing more. An independent version needs
  the check variable kept out of the inputs, e.g. inflammation from NLR only, validated against CRP, or the
  GDF15 → `Aging_secretory` check (needs a GDF15 weight table). **Done since:** `inflammation` is now
  built from N:L only by default and CRP stays on the validation side, so the check is independent
  (`--inflammation-with-crp` restores the old, overlapping definition for comparison). The 0.86–0.88
  above was measured with the old definition and must be re-run. Caveat: CRP EpiScores partly reflect
  leukocyte composition, so some shared variance with N:L is expected even when independent by
  construction.
- **Age-acceleration ≈ 0.35** (index residual vs clock residual, in-sample; same selection caveat). This is
  the genuinely non-circular result: nothing from the clock enters the model.

**Recommendation.** Chronological age is a weak target for this model; report **age-acceleration** as
primary, add a truly independent sub-score check (see above), confirm the chosen subset out-of-sample, and
benchmark against health outcomes next. Per dataset, run
`sweep_configs.py` to pick the feature subset, then `geo_validation.py --use-features <subset>`.

## Feasibility note
Steps 0–1 (download, then EpiDISH + EpiScores + clock) are a standard but heavy R job (GEOquery/EpiDISH/
methylclock + reference panels + a large matrix, ~16 GB RAM; the `--test` run needs only a few GB). Step 2 reuses this repo unchanged and
runs anywhere Python does. Run step 1 where the data + R live; step 2 on a laptop.
