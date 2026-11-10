# EpiScore weight tables — `src/episcore_weights/`

This folder holds the CpG weight tables that turn the blood 450K matrix into the model's
protein/exposure features. Cell composition (EpiDISH) and the clock (methylclock) run without any of
these; the five feature files below are what light up `inflammation`, `adiponectin`, `adiposity_lep`,
`hepatic`, and `damage_smoke`.

**Location matters.** `preprocess_gse40279.R` lives in `src/` and is run from `src/`, so it looks for
`episcore_weights/` *relative to `src/`* — i.e. **this folder, `src/episcore_weights/`**. (Putting it at
the repo root is the classic "weights folder not found → all EpiScores NA" trap.)

## File format

One CSV per score, two columns, header required. The score for a sample is
`intercept + Σ weight·beta[cpg]` over CpGs present in both the table and the matrix.

```csv
cpg,weight
cg25124965,-0.010712146
cg06629775,0.012668712
(Intercept),2.7100        # optional
```

## Fastest way to produce them — MethylDetectR + the converter

The MethylDetectR bundle is the best single source: it is 450K-oriented and contains the protein
EpiScores **plus** lifestyle/adiposity/smoking predictors, in exactly the long format the converter
reads (`CpG_Site,Coefficient,Mean_Beta_Value,Predictor`).

```bash
# 1. download (≈675 kB)
curl -L -o Predictors_Shiny_by_Groups.csv \
  "https://zenodo.org/records/15548022/files/Predictors_Shiny_by_Groups.csv?download=1"

# 2. (recommended) in R, dump your array's CpGs once:  writeLines(rownames(beta), "beta_cpgs.txt")

# 3. generate the five files INTO this folder, checking 450K overlap as it goes:
cd src
python3 make_episcore_weights.py Predictors_Shiny_by_Groups.csv \
    --outdir episcore_weights --beta-cpgs beta_cpgs.txt \
    --map CRP.csv="CRP" ADIPOQ.csv="Adiponectin" LEP.csv="Body Fat %" \
          liver.csv="HGF" smoking.csv="Smoking"
```

### Recommended mapping (confirmed present in MethylDetectR)

| file | predictor | CpGs | model feature | note |
|---|---|---|---|---|
| `CRP.csv`     | `CRP`         | 100 | `inflammation` | combined with the neutrophil:lymphocyte ratio |
| `ADIPOQ.csv`  | `Adiponectin` | 39  | `adiponectin`  | |
| `LEP.csv`     | `Body Fat %`  | 968 | `adiposity_lep`| **substitute** — true leptin isn't in the bundle; body-fat % drives the same LEP/Nutrients inputs |
| `liver.csv`   | `HGF`         | 72  | `hepatic`      | HGF rises with liver injury (correct sign); avoid SHBG, which is *inversely* related to liver fat |
| `smoking.csv` | `Smoking`     | 233 | `damage_smoke` | no separate EpiSmokEr needed |

Not available in MethylDetectR / the Gadd files: **GDF15, IL6, true leptin**. Consequences:
`IL6` doesn't matter (inflammation runs on CRP + N:L); leptin is covered by Body Fat %; **GDF15 is
validation-only** — the `Aging_secretory` vs GDF15 non-circular check simply won't run unless you add a
`GDF15.csv` from another source. `NAD_CD38.csv` is optional (NAD depletion is an emergent model state).

This takes the feature matrix from 2/7 to **6/7 live** (all required), so `check_features.py` passes.

## Alternative sources

- **Gadd 2022 raw protein weights** — `LBC_weights_finalised_120121.csv` (450K-native but only 25
  proteins, no CRP/adiponectin → too thin) and `KORA_weights_finalised_120121.csv` (EPIC/SomaScan, 84
  proteins incl. CRP/Adiponectin/SHBG — usable, but verify 450K overlap and it has no lifestyle scores).
  The converter reads both (`CpG.site`/`Coefficient`/`Predictor` are auto-detected). From
  [DanniGadd/EpiScores-for-protein-levels](https://github.com/DanniGadd/EpiScores-for-protein-levels)
  `Weights_files/` or eLife Supplementary file 1E.
- **Smoking** — [EpiSmokEr](https://github.com/sailalithabollepalli/EpiSmokEr) or McCartney 2018.

## Sources & attribution (according to PubMed)

- Protein EpiScores: Gadd et al. 2022, *eLife* [doi:10.7554/eLife.71802](https://doi.org/10.7554/eLife.71802) (PMID 35023833).
- Lifestyle / BMI / smoking DNAm scores: McCartney et al. 2018, *Genome Biol* [doi:10.1186/s13059-018-1514-1](https://doi.org/10.1186/s13059-018-1514-1) (PMID 30257690).
- MethylDetectR bundle: Zenodo record [10.5281/zenodo.15548022](https://zenodo.org/records/15548022) (CC-BY-4.0).
- Clock (external only): Hannum et al. 2013, *Mol Cell* [doi:10.1016/j.molcel.2012.10.016](https://doi.org/10.1016/j.molcel.2012.10.016) (PMID 23177740).

## Version control

The weight CSVs are **gitignored by default** (third-party derived data; licensing varies). Regenerate
them with the converter as above. If you prefer to vendor them in the repo for turnkey reproducibility
(MethylDetectR is CC-BY-4.0 — keep this attribution), force-add them:

```bash
git add -f src/episcore_weights/CRP.csv src/episcore_weights/ADIPOQ.csv src/episcore_weights/LEP.csv \
           src/episcore_weights/liver.csv src/episcore_weights/smoking.csv
```
