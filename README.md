# Prevention & Longevity Map

A mechanistic, executable model of aging that turns a blood DNA-methylation profile into an
interpretable, patient-specific readout and an intervention ranking. A primary-literature-curated
Boolean network (CaSQ → MaBoSS) couples a **nutrient-sensing / NAD⁺ core (P1)** to a
**muscle / myokine module (MUS)**, linked through a shared inflammation hub and a hepatic
(hyperammonemia / FGL1) bridge. On top sits an **aging layer** — seven hallmark sub-scores and an
integrated index, built from pathway activity and metabolites (not from clocks, so the clocks stay an
independent external check). A **methylation input layer** parameterizes the model per patient.

Status: isolated validation done; methylation layer prototype done; GEO validation pipeline proven
in-silico (aging index vs age r ≈ 0.63) and **wired for real data** (`download_geo.R` → `preprocess_geo.R` →
`geo_validation.py --features`); organ modules (liver, stem-cell/CHIP) scaffolded.

## Repository layout

```
src/            runnable model + analysis (run scripts from here)
  p1_model.py            P1 nutrient-sensing/NAD core (node/edge gates)
  coupled_model.py       P1 + MUS + recalibrated phenotype layer + aging layer  -> writes MUS_P1_coupled.bnd/.cfg
  recal.py               recalibration harness (anchor scenarios before/after)
  aging_layer.py         hallmark sub-scores across archetypes + interventions
  methylation_layer.py   methylation features -> input fractions -> per-patient readout
  run_coupled.py         intervention scenarios on the coupled model
  robustness_coupled.py  ±30% weight-noise sign-stability check
  characterize.py        systematic behaviour battery (inputs × readouts, dose, interactions)
  make_figs.py           behaviour-map + dose-response figures
  geo_validation.py      cohort validation: simulated (default) or real (--features features_<gse>.csv)
  download_geo.R         R: fetch ONE named GEO cohort once (--gse40279 / --gse87571) -> cached .rds
  preprocess_geo.R       R: cached matrix -> features_<gse>.csv (EpiDISH + EpiScores + clock; --by-sex)
  cohorts.R              the known cohorts and their roles (development / validation)
  check_features.py      preflight gate: per-feature coverage + age-correlation on a features file
  sweep_configs.py       feature-subset x aggregation sweep (development cohort only)
  make_episcore_weights.py  split a MethylDetectR/Gadd weights table -> per-score episcore_weights/*.csv
  episcore_weights/      EpiScore CpG weight tables (see its README; weights gitignored)
  MUS_module_v3.json     curated muscle edges (191; 149 in executable model)
  P1_curation.json       curated P1 edges
data/           result JSONs (characterization, scenarios, robustness, aging, methylation, geo)
modules/        liver + stem-cell/CHIP module scaffolds (collaboration starters)
validation/     geo_validation_pipeline.md — real-data recipe (GSE40279 dev, GSE87571 validation)
figures/        PNG/SVG figures and schematics
docs/           schematics (SVG)
```

## Quickstart

```bash
uv venv --python 3.13 && source .venv/bin/activate   # or plain venv / conda, see SETUP.md
uv pip install -r requirements.txt
cd src
python3 coupled_model.py                 # build the 98-node model -> MUS_P1_coupled.bnd/.cfg
python3 run_coupled.py                   # intervention scenarios
python3 methylation_layer.py             # per-patient demo (4 synthetic profiles)
python3 geo_validation.py                # simulated-cohort validation (aging index vs age)

```

### Real data: two cohorts, kept apart

Run from `src/`, where the matrix and ~16 GB RAM live. The cohort is always named explicitly, and each
gets its own data folder (`data/raw/<gse>/`) and its own `features_<gse>.csv`, so they never mix.
EpiScore weights ship in `src/episcore_weights/` (to rebuild them, see its README).

**1. Development cohort: GSE40279** (Hannum 2013, n = 656). All tuning happens here: feature
subset, index aggregation, sex handling.

```bash
Rscript download_geo.R --gse40279                     # once -> ../data/raw/gse40279/
Rscript preprocess_geo.R --gse40279 --test            # quick check on a 48-sample subset
Rscript preprocess_geo.R --gse40279 --by-sex          # -> features_gse40279.csv (within-sex percentiles)
python3 check_features.py features_gse40279.csv       # preflight: fails fast on empty inputs
python3 sweep_configs.py --features features_gse40279.csv --subsample 160   # pick subset / aggregation
python3 geo_validation.py --features features_gse40279.csv --by-sex --use-features <subset>
```

**2. Validation cohort: GSE87571** (Johansson 2013, ~730, ages ~14–94). Held out: run it **once**
with the settings frozen on GSE40279, and never tune on it (`sweep_configs.py` refuses this cohort).

```bash
Rscript download_geo.R --gse87571                     # once -> ../data/raw/gse87571/
Rscript preprocess_geo.R --gse87571 --min-age 19 --by-sex   # same options as the dev run, + age range match
python3 check_features.py features_gse87571.csv
python3 geo_validation.py --features features_gse87571.csv --by-sex --use-features <same subset>
```

Clock note: the Hannum clock was trained on GSE40279, so its age acceleration there is in-sample.
`preprocess_geo.R --clock Horvath` writes Horvath to `clock_age` instead (both clocks are always kept as
`clock_hannum` / `clock_horvath`); verify the clock you pick was not trained on the cohort you analyse.

Everything model-side is lightweight (pure Python + `cmaboss`); a full run is seconds to a couple of
minutes on a laptop. The only heavy, R-based step is turning a **real** methylation matrix into model
features — see `SETUP.md` and `validation/geo_validation_pipeline.md`.

Evidence throughout is PubMed-curated (primary or canonical); DOIs are in the module JSONs and the
analysis workbooks.
