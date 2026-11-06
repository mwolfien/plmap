# Standards & interoperability (SBML / SBGN / SED-ML / COMBINE)

For publication and FAIR reuse the model is provided in community standards alongside the exact MaBoSS
engine files. The guiding split: **standard formats carry the structure and logic; the MaBoSS files
carry the exact weighted dynamics.** Both ship together in the COMBINE archive.

| File | Standard | Carries | Honest caveat |
|---|---|---|---|
| `PLM_model.sbml` | **SBML L3 + qual** | 98 nodes, 188 **signed** regulatory edges, canonical sign-consistent Boolean function | the exact weighted transition *rates* are MaBoSS-specific and are **not** representable in SBML-qual — they live in the `.bnd`/`.cfg` |
| `MUS_P1_coupled.bnd` | MaBoSS network | exact executable dynamics (weighted `rate_up`/`rate_down`) | MaBoSS-native format |
| `MUS_P1_coupled.cfg` | MaBoSS config | **the authoritative simulation setup** (max_time, time_tick, sample_count, initial states) | MaBoSS-native (the native analogue of SED-ML) |
| `simulation.sedml` | **SED-ML L1V4** | time-course description mirroring the `.cfg` | MaBoSS↔SED-ML is not push-button; run via BioSimulators; **confirm the BKMC KISAO term** before submission |
| `PLM_model.omex` | **COMBINE / OMEX** | all of the above + `metadata.rdf` manifest | the single-file FAIR artifact |
| `PLM_network.sif`, `PLM_edges.tsv`, `PLM_nodes.tsv` | SIF / tables | signed influence graph for quick import | structure only |

## SBML-qual
Generated with libSBML 5.21.2 (qual package). It is the interoperable structural model: the
GINsim / bioLQM / MINERVA / CaSQ / Cell Collective ecosystem all read SBML-qual. Validated
(`checkConsistency` → 0 errors). Each transition lists its regulators as signed inputs
(`positive`/`negative`) with a canonical Boolean reading: a node is active when an activator is active
and no inhibitor is active; basal-driven nodes are active unless repressed. **This reproduces the
model's wiring and sign logic, not its exact weighted kinetics** — for the latter use the MaBoSS files.

## SBGN
A Boolean regulatory model maps naturally to **SBGN Activity Flow (SBGN-AF)** — nodes are biological
activities, arcs are positive/negative influence. Two easy routes:
1. Import `PLM_model.sbml` into **GINsim** (export → SBGN-AF) or **MINERVA** (renders qual models in
   SBGN), or **Cell Collective**.
2. Import `PLM_network.sif` + `PLM_edges.tsv` into **Cytoscape** (→ SBGN via `cytoscape-cd2sbgnml` or
   the Newt editor) for a quick signed-network figure.

The upstream NSCLC disease map this project extends is already in **CellDesigner / SBGN-PD**, so the
detailed molecular layer and this executable layer are both SBGN-expressible (PD for the map, AF for
the Boolean model).

## SED-ML & COMBINE
`simulation.sedml` describes the time course (initialTime 0, outputEndTime 40, 80 steps — mirroring the
`.cfg`). MaBoSS's stochastic trajectory count (`sample_count`) has no SED-ML equivalent and stays in the
`.cfg`. Execution through SED-ML goes via **BioSimulators** (which packages MaBoSS as a container with a
registered KISAO term); the template's KISAO id is flagged for verification. `PLM_model.omex` bundles
SBML-qual + MaBoSS files + SED-ML + metadata with an OMEX manifest — open it with any COMBINE-aware tool
(or just unzip it).

## Reproduce the exports
```bash
cd src
python3 export_standards.py     # -> PLM_model.sbml, PLM_network.sif, PLM_edges.tsv, PLM_nodes.tsv
python3 build_omex.py           # -> simulation.sedml, manifest.xml, metadata.rdf, PLM_model.omex
```
Regenerate these whenever the model changes (add a module, recalibrate) so the standard exports stay in
sync with `coupled_model.py`.
