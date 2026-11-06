# EpiScore weight tables (user-supplied)

`preprocess_gse40279.R` turns the blood 450K matrix into model features. Cell composition (EpiDISH)
and the epigenetic clock (methylclock) run **out of the box**. The protein / exposure **EpiScores**
need published CpG weight tables, which are not redistributed here (licensing / size varies by source).
Drop one CSV per score into this folder; any table that is absent → that feature is written `NA` and
imputed to the cohort median downstream (the Python step reports which).

## File format

One CSV per score, named exactly as below. Two columns, header required:

```csv
cpg,weight
cg00000029,0.1234
cg00000108,-0.0456
(Intercept),2.7100        # optional; omit if the source gives no intercept
```

The score for a sample is `intercept + Σ_cpg weight[cpg] · beta[cpg, sample]`, over the CpGs present in
both the table and the matrix (the script reports the overlap count). Use the published weights as-is;
do not re-scale — percentile normalization happens afterwards in the script.

## Which table → which file → which model feature

| file | EpiScore | feeds model feature | source |
|---|---|---|---|
| `CRP.csv`     | C-reactive protein     | `inflammation` (with IL6 + N:L ratio) | Gadd 2022 |
| `IL6.csv`     | Interleukin-6          | `inflammation` | Gadd 2022 |
| `GDF15.csv`   | GDF15                  | **validation only** (`val_gdf15_episcore` → checks `Aging_secretory`) | Gadd 2022 |
| `ADIPOQ.csv`  | Adiponectin            | `adiponectin` (protective; model `ADIPOQ` input) | Gadd 2022 |
| `LEP.csv`     | Leptin                 | `adiposity_lep` | Gadd 2022 |
| `liver.csv`   | Liver enzyme / hepatic protein (e.g. a GGT- or ALT-associated EpiScore) | `hepatic` | Gadd 2022 |
| `smoking.csv` | Smoking score          | `damage_smoke` | McCartney 2018 |
| `NAD_CD38.csv`| CD38 / NAD surrogate (optional) | `nad_depletion` | none standard — see note |

CRP and IL6 together also form `val_inflammation_episcore` (a non-circular check of `Aging_inflammation`).

## Sources

According to PubMed:

- **Protein EpiScores** — Gadd DA, Hillary RF, McCartney DL, et al. *Epigenetic scores for the
  circulating proteome as tools for disease prediction.* eLife 2022;11:e71802.
  [DOI](https://doi.org/10.7554/eLife.71802) · PMID 35023833. 109 protein EpiScores (CRP, IL6, GDF15,
  adiponectin/ADIPOQ, leptin/LEP, and liver-associated proteins among them); CpG weights are in the
  supplementary. The **MethylDetectR** resource (Hillary & Marioni) bundles several of these
  ready-to-use.
- **Smoking** — McCartney DL, Hillary RF, Stevenson AJ, et al. *Epigenetic prediction of complex
  traits and death.* Genome Biol 2018;19:136. [DOI](https://doi.org/10.1186/s13059-018-1514-1) ·
  PMID 30257690. (DNAm predictors for smoking, BMI, WHR, alcohol — smoking → `damage_smoke`;
  BMI/WHR can substitute for `adiposity_lep` if you prefer an anthropometric score to leptin.)
- **Clock (external, not an input)** — computed in-script via the `methylclock` package
  (Horvath / Hannum). Hannum G, et al. Mol Cell 2013;49:359-367.
  [DOI](https://doi.org/10.1016/j.molcel.2012.10.016) · PMID 23177740.

## Note on `nad_depletion`

There is no established, widely-validated blood **NAD/CD38 EpiScore**. Options, in order of preference:
1. leave `NAD_CD38.csv` absent → `nad_depletion` is imputed to the cohort median (neutral);
2. supply a CD38-promoter / NAD-pathway surrogate table if your group has one;
3. train a small EpiScore against measured NAD⁺/CD38 in a reference set, then drop its weights here.

The model's NAD arm is still exercised through the inflammation/CD38 coupling even when this feature is
imputed; a dedicated table only sharpens it.
