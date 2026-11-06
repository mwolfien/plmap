# GEO validation — real-data pipeline (GSE40279)

**Goal.** Test whether the coupled model's aging index / hallmark sub-scores, parameterized from
*real* blood-methylation features, track chronological age and an epigenetic clock.

**Why GSE40279** (Hannum 2013, PMID 23177740, [doi:10.1016/j.molcel.2012.10.016](https://doi.org/10.1016/j.molcel.2012.10.016)):
656 whole-blood Illumina 450K samples, ages 19–101, public, widely used for clocks — the cleanest
public match for a blood-readout aging model. Alternatives: GSE87571 (wide age range), GSE55763
(large), GSE132203 (with phenotypes).

**Status.** The machinery is proven on a literature-calibrated simulated cohort (`geo_validation.py`,
`fig_geo_validation.png`): the aging index tracks simulated age (r ≈ 0.73), NAD-decline and
inflammaging sub-scores strongest. The step below swaps the simulated feature vectors for real ones.
It needs the beta matrix + an R methylation stack, so run it where the data lives (or stage it here).

## Steps

1. **Get the data.** Download the GSE40279 processed beta matrix (GEO supplementary / series matrix)
   and the sample table (age, sex). ~656 × ~473k CpGs.

2. **Derive the methylation features** the layer needs (all standard, blood-validated):
   - **Cell composition** — `EpiDISH` (RPC, `centDHSbloodDMC.m` reference): CD4T, CD8T, NK, B, Mono, Neu.
     Derived features: neutrophil:lymphocyte ratio, CD4:CD8, naive:memory proxy → `inflammation`,
     `immune_aging`.
   - **Protein EpiScores** — Gadd 2022 / Hillary `MethylDetectR` weights: CRP, IL-6, GDF15, adiponectin,
     leptin, liver enzymes (GGT/ALT) → `inflammation`, `hepatic`, `adiposity_lep`, `adiponectin`,
     and the GDF15 EpiScore to *validate* the secretory sub-score.
   - **Smoking EpiScore** (e.g. Elliott/`EpiSmokEr`) → `damage_smoke`.
   - **Epigenetic clock** — Horvath and/or Hannum (and TFMethyl once available from Patel/Hoffmann) as
     the **external** readout to compare against — never an input.

3. **Normalize to percentiles.** Convert each feature to its cohort percentile (0–1); that is the
   value `methylation_layer.features_to_fractions` expects. (Keep the reference cohort fixed so a new
   patient is placed against it.)

4. **Run each sample through the model.** For each subject, `features_to_fractions(feat)` → input
   `istate` fractions → `cmaboss_128n` run → read `Aging_index`, the 7 sub-scores, `Sarcopenia`,
   `Insulin_resistance`, `MSTN`, `IL6` (exactly as `geo_validation.py` does, replacing the simulated
   `feat` with the real one).

5. **Validate.**
   - Correlate model `Aging_index` with chronological age and with the clock (Pearson/Spearman).
   - Per sub-score: does `Aging_NAD` track age, `Aging_inflammation` track the CRP/IL-6 EpiScore,
     `Aging_secretory` track the GDF15 EpiScore (a direct, non-circular check of that sub-score)?
   - Age-acceleration test: residual of clock-on-age vs residual of index-on-age — do the people the
     clock calls "old for their age" also score high on specific mechanistic hallmarks?
   - Report correlations + the scatter (as in `fig_geo_validation.png`).

6. **Interpret.** Strong index–age correlation = the mechanistic layer recovers aging from real
   methylation. Sub-score-specific hits (GDF15 EpiScore → secretory sub-score; CRP → inflammaging)
   are the stronger, non-circular validations. Divergence between clock age and index flags where the
   model is incomplete (e.g. telomere/epigenetic-drift hallmarks not yet in the network).

## Feasibility note
Steps 1–2 (download + EpiDISH + EpiScores + clock) are a standard but heavy R job (minfi/EpiDISH +
reference panels, a large matrix). Steps 3–5 reuse `methylation_layer.py` unchanged. The honest
split: run 1–2 where the data and R stack live; steps 3–5 run anywhere the model does.
