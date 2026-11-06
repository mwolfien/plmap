# Prevention & Longevity Map

A mechanistic, executable model of aging that turns a blood DNA-methylation profile into an
interpretable, patient-specific readout and an intervention ranking. A primary-literature-curated
Boolean network (CaSQ → MaBoSS) couples a **nutrient-sensing / NAD⁺ core (P1)** to a
**muscle / myokine module (MUS)**, linked through a shared inflammation hub and a hepatic
(hyperammonemia / FGL1) bridge. On top sits an **aging layer** — seven hallmark sub-scores and an
integrated index, built from pathway activity and metabolites (not from clocks, so the clocks stay an
independent external check). A **methylation input layer** parameterizes the model per patient.

Status: isolated validation done; methylation layer prototype done; GEO validation pipeline proven
in-silico (aging index vs age r ≈ 0.73) and **wired for real data** (`preprocess_gse40279.R` →
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
  geo_validation.py      cohort validation: simulated (default) or real (--features features.csv)
  preprocess_gse40279.R  heavy R step: raw GEO 450K -> features.csv (EpiDISH + EpiScores + clock)
  MUS_module_v3.json     curated muscle edges (191; 149 in executable model)
  P1_curation.json       curated P1 edges
data/           result JSONs (characterization, scenarios, robustness, aging, methylation, geo)
modules/        liver + stem-cell/CHIP module scaffolds (collaboration starters)
validation/     geo_validation_pipeline.md — real-data (GSE40279) recipe
figures/        PNG/SVG figures and schematics
docs/           schematics (SVG)
```

## Quickstart

```bash
pip install -r requirements.txt          # or: conda env create -f environment.yml
cd src
python3 coupled_model.py                 # build the 98-node model -> MUS_P1_coupled.bnd/.cfg
python3 run_coupled.py                   # intervention scenarios
python3 methylation_layer.py             # per-patient demo (4 synthetic profiles)
python3 geo_validation.py                # simulated-cohort validation (aging index vs age)

# real data (heavy R step first, where the matrix + ~16GB RAM live):
Rscript preprocess_gse40279.R --weights episcore_weights   # -> features.csv
python3 geo_validation.py --features features.csv           # -> fig_geo_validation_real.png
```

Everything model-side is lightweight (pure Python + `cmaboss`); a full run is seconds to a couple of
minutes on a laptop. The only heavy, R-based step is turning a **real** methylation matrix into model
features — see `SETUP.md` and `validation/geo_validation_pipeline.md`.

Evidence throughout is PubMed-curated (primary or canonical); DOIs are in the module JSONs and the
analysis workbooks.
