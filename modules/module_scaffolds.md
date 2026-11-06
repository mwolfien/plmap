# Organ-module scaffolds — collaboration starters

Two lightweight module scaffolds that attach to the existing validated core through the shared
inflammation hub / secretome. **Anchor edges are PubMed-verified primary papers; proposed edges are
marked and left for the collaborator to curate** — the point is a concrete frame to co-develop, not a
finished module. Full edge lists: `liver_module_scaffold.json`, `chip_module_scaffold.json`.

## Liver / metabolic module (LIV) — for Carolin Schneider

Attaches to: Insulin_resistance, AMPK, NF-κB, MSTN, GDF15, hyperammonemia, FGL1, LECT2.

Verified anchor edges:

- **Fetuin-A (AHSG) → TLR4 → insulin resistance** — Fetuin-A is the endogenous TLR4 ligand for lipid-induced IR (Pal 2012, [10.1038/nm.2851](https://doi.org/10.1038/nm.2851)).
- **Selenoprotein P ⊣ AMPK → insulin resistance** — liver-derived SeP inactivates AMPK (Misu 2010, [10.1016/j.cmet.2010.09.015](https://doi.org/10.1016/j.cmet.2010.09.015)).
- **FGL1 (hepassocin) → insulin resistance** — already in the muscle module (Jung 2017, [10.1016/j.mce.2017.10.014](https://doi.org/10.1016/j.mce.2017.10.014)).
- **Hyperammonemia → myostatin** — the hepatic→muscle sarcopenia bridge already in the map (Qiu 2013, [10.1073/pnas.1317049110](https://doi.org/10.1073/pnas.1317049110)).

Proposed (to curate): MASLD → hepatokine secretion; LECT2 → IR (Lan 2014, PMID not retrievable here); IR → DNL → steatosis feed-forward; cirrhosis → urea-cycle failure → hyperammonemia; fibrosis progression; bile-acid/FXR axis (intervention node).

**Co-develop:** complete the hepatokine panel and targets; the ammonia–urea–myostatin sarcopenia axis (your domain); MASLD progression + intervention nodes (FXR agonists, resmetirom).

## Stem-cell / clonal haematopoiesis module (CHIP) — for Ingmar Glauche

Attaches to: IL-6, IL-1β, NF-κB, Inflammaging, DNA_damage, and the Aging_inflammation / Aging_senescence_genomic sub-scores.

Verified anchor edges:

- **CHIP clone (DNMT3A/TET2) expands with age → adverse outcomes** (Jaiswal 2014, [10.1056/NEJMoa1408617](https://doi.org/10.1056/NEJMoa1408617); CVD risk Jaiswal 2017, [10.1056/NEJMoa1701719](https://doi.org/10.1056/NEJMoa1701719)).
- **TET2 loss → IL-6** — TET2 normally represses IL-6 via HDAC2; its loss de-represses IL-6 (Zhang 2015 Nature, [10.1038/nature15252](https://doi.org/10.1038/nature15252)).
- **TET2 CHIP → NLRP3 → IL-1β** — Tet2-deficient clonal haematopoiesis drives IL-1β/NLRP3 (Sano 2018, [10.1016/j.jacc.2017.12.037](https://doi.org/10.1016/j.jacc.2017.12.037)).
- **TET2 CHIP → inflammaging / atherosclerosis**, IL-6-dependent (Fuster 2017, [10.1126/science.aag1381](https://doi.org/10.1126/science.aag1381)).

Proposed (to curate): DNA damage → somatic mutation → clone; inflammation (IL-6/TNF) → selects mutant clones (feed-forward); clonal expansion → myeloid bias → ↓lymphoid output → immune aging; JAK2 axis.

**Co-develop:** the clonal-dynamics layer (fitness, expansion rate vs age) — your HSC-modelling expertise would make the single "clonal expansion" node quantitative; the inflammation↔clone feed-forward loop; wiring CHIP output to the methylation immune-aging features and aging sub-scores.

---
Evidence from PubMed; DOIs above. Each module plugs into a working, validated core — so each collaborator extends a running system rather than starting from a blank network.
