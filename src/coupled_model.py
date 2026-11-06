"""Coupled P1 (nutrient-sensing/NAD+) + MUS (muscle/sarcopenia) MaBoSS model.
Reuses the P1 core N-dict and INPUTS unchanged; adds a muscle layer whose gates are
hand-curated over the PRIMARY-evidence MUS v2 edges (model=='include'/'use P1').
Weighted-rate formulation (same as P1): rate_up = basal + sum(w*activator); rate_down = basal + sum(w*inhibitor).
Uses cmaboss_128n (P1 core ~60 nodes + muscle layer > 64).
Each muscle gate names the supporting MUS rid + PMID verified this session.
"""
import sys, json
pass  # p1_model co-located
import p1_model as P1

INPUTS = list(P1.INPUTS) + ['Exercise', 'leucine', 'FST', 'hyperammonemia', 'glucocorticoids',
                            'LEP', 'ADIPOQ', 'FGL1', 'Omega3_PUFA']
N = dict(P1.N)  # P1 hub reused verbatim

# ---- muscle layer (labels harmonised to P1 node names; PPARGC1A->PGC1A, NF-kB->NFkB) ----
MUS = {
 # anabolic input arms onto the P1 hub
 'leucine_mTOR': None,  # placeholder removed below
 'Testosterone': ([(0.1, 1), (1.5, 'Exercise')], [(1, 1)]),                       # MUS-163 (Exercise->T)
 'VDR':          ([(0.2, 1), (1, 'Vitamin_D')], [(1, 1)]),                         # MUS-046
 'Vitamin_D':    ([(0.3, 1)], [(1, 1)]),                                           # input-like (diet/sun), kept internal
 # myostatin / TGF-beta atrophy arm
 'FGF21_mito':   None,
 'Mito_dysfunction': ([(0.2, 1), (1, 'ROS'), (1.5, 'hyperammonemia')], [(0.3, 1), (2, 'PGC1A')]),  # MUS-053, MUS-082; MUS-055 (PGC1A-|)
 'MSTN':         ([(0.1, 1), (1.5, 'NFkB')], [(0.3, 1), (2, 'FST'), (1, 'Testosterone')]),         # MUS-023(primary), MUS-105(primary FST), MUS-084 weak
 'ACVR2':        ([(1, 'MSTN')], [(0.5, 1)]),                                      # MUS-022 (primary Lee2001)
 'smad2_3':      ([(0.1, 1), (1.5, 'ACVR2')], [(0.5, 1)]),                         # MSTN/ActRIIB->Smad2/3 (Sartori/Trendelenburg 2009, primary)
 # proteolysis machinery
 'TRIM63':       ([(0.1, 1), (1, 'NFkB'), (1, 'FOXO3')], [(0.5, 1)]),             # MUS-065 (primary Cai2004); FOXO3 atrogene (Sandri2004)
 'UPS':          ([(0.1, 1), (1, 'NFkB'), (1, 'TRIM63'), (1, 'smad2_3'), (1.5, 'glucocorticoids')], [(0.6, 1)]),  # MUS-024(primary), MUS-096(primary), MUS-021
 'Proteolysis':  ([(0.05, 1), (1.5, 'UPS'), (1, 'Autophagy'), (1, 'FOXO3')], [(0.5, 1)]),           # MUS-006, MUS-189(primary), MUS-008(primary)
 # anabolic output: Protein_synthesis is a P1 node (mTOR/S6K1 driven); add muscle activators
 'Protein_synthesis_mus': None,
 # --- v3 myokine / adipokine / hepatokine layer (newly primary-backed regulatory edges) ---
 'TNF':          ([(0.1, 1), (1, 'LEP'), (1.5, 'Inflammaging')], [(0.5, 1)]),      # MUS-174 (LEP->TNF, primary); Inflammaging driver
 'IL6':          ([(0.1, 1), (1, 'LEP'), (1, 'Inflammaging')], [(0.3, 1), (1, 'Omega3_PUFA')]),  # MUS-175 (LEP->IL6), MUS-153 (PUFA-|IL6)
 'FNDC5':        ([(0.1, 1), (1.5, 'Exercise')], [(1, 1)]),                        # irisin, exercise-induced (MUS-116/157)
 'METRNL':       ([(0.1, 1), (1.5, 'Exercise')], [(1, 1)]),                        # MUS-157 (exercise myokine)
 'IL10':         ([(0.1, 1), (1, 'METRNL')], [(1, 1)]),                            # MUS-126 (METRNL->IL10, anti-inflammatory)
 'IL15':         ([(0.1, 1), (1, 'Exercise')], [(1, 1)]),                          # MUS-157
 'DCN':          ([(0.1, 1), (1.5, 'Exercise')], [(1, 1)]),                        # MUS-135 (Exercise->DCN, primary)
 'Insulin_resistance': ([(0.1, 1), (1, 'FGL1'), (1, 'TNF'), (1, 'NFkB')], [(0.3, 1), (1.5, 'ADIPOQ'), (1, 'Insulin_sensitivity')]),  # MUS-049 (FGL1, primary), MUS-088 (ADIPOQ-|, sign fix)
 # phenotype aggregation (structural glue; definitional)
 'Anabolism':    ([(0.1, 1), (1.5, 'Protein_synthesis')], [(0.5, 1)]),            # MUS-155
 'Catabolism':   ([(0.1, 1), (1.5, 'Proteolysis'), (1, 'Mito_dysfunction')], [(0.5, 1)]),  # MUS-001/156
 'Sarcopenia':   ([(0.1, 1), (2, 'Catabolism')], [(0.3, 1), (2, 'Anabolism')]),  # MUS-002, MUS-156
}
del MUS['leucine_mTOR']; del MUS['FGF21_mito']; del MUS['Protein_synthesis_mus']

N.update(MUS)

# ---- rewire a few P1 hub nodes to accept muscle inputs (extend, do not replace, P1 logic) ----
# mTOR: + leucine (MUS-140 primary 11015466/18682538), + BCAAs via leucine proxy
up, dn = N['mTOR']; N['mTOR'] = (up + [(1, 'leucine')], dn)
# AMPK: + Exercise/ATP-shortage (MUS-173 primary), leucine -| AMPK (MUS-141)
up, dn = N['AMPK']; N['AMPK'] = (up + [(1, 'Exercise')], dn + [(0.5, 'leucine')])
# Akt: smad2_3 -| Akt (MUS-020 primary Sartori2009); IGF1/Testosterone arm via P1 IGF1
up, dn = N['Akt']; N['Akt'] = (up, dn + [(1.5, 'smad2_3')])
# FOXO3: + muscle ROS arm (MUS-071 primary Brunet2004) already via SIRT1; AMPK->FOXO3 (MUS-185)
up, dn = N['FOXO3']; N['FOXO3'] = (up + [(0.5, 'AMPK')], dn)
# PGC1A: + Exercise (MUS-113 primary), smad2_3 -| PGC1A (MUS-184)
up, dn = N['PGC1A']; N['PGC1A'] = (up + [(1, 'Exercise')], dn + [(1, 'smad2_3')])
# NFkB: + hyperammonemia (MUS-081 primary Qiu2013)
up, dn = N['NFkB']; N['NFkB'] = (up + [(1.5, 'hyperammonemia')], dn)
# FGF21: + Mito_dysfunction (MUS-139 primary, conf5)
up, dn = N['FGF21']; N['FGF21'] = (up + [(1, 'Mito_dysfunction')], dn)
# IGF1: + Testosterone (MUS-089 primary)
up, dn = N['IGF1']; N['IGF1'] = (up + [(0.7, 'Testosterone')], dn)
# Protein_synthesis: + leucine & VDR (MUS-140/180); + irisin/FNDC5 (MUS-117 reroute, Reza2017)
up, dn = N['Protein_synthesis']; N['Protein_synthesis'] = (up + [(0.5, 'leucine'), (0.3, 'VDR'), (0.5, 'FNDC5')], dn)
# --- v3 myokine/adipokine rewires onto the muscle-mass hub ---
# MSTN: + decorin antagonism (MUS-134, primary)
up, dn = N['MSTN']; N['MSTN'] = (up, dn + [(1, 'DCN')])
# UPS: -| IL-15 (MUS-132, primary Busquets/Quinn)
up, dn = N['UPS']; N['UPS'] = (up, dn + [(1, 'IL15')])
# Proteolysis: -| irisin/FNDC5 (MUS-118 reroute, Reza2017)
up, dn = N['Proteolysis']; N['Proteolysis'] = (up, dn + [(0.8, 'FNDC5')])
# NFkB: + TNF (MUS-099, primary); ROS: + TNF (MUS-033, primary Li&Reid)
up, dn = N['NFkB']; N['NFkB'] = (up + [(1, 'TNF')], dn)
up, dn = N['ROS']; N['ROS'] = (up + [(0.7, 'TNF')], dn)
# AMPK: + adiponectin (MUS-088 sign fix, Yamauchi2002)
up, dn = N['AMPK']; N['AMPK'] = (up + [(0.8, 'ADIPOQ')], dn)
# Insulin_sensitivity: -| insulin resistance arm (FGL1/TNF/ADIPOQ integrated in Insulin_resistance)
up, dn = N['Insulin_sensitivity']; N['Insulin_sensitivity'] = (up, dn + [(1, 'Insulin_resistance')])

# ---- RECALIBRATION of the aggregation stack (bottom-up, anchor-tuned; see recal.py) ----
# Fixes the compressed Sarcopenia/Catabolism readouts: re-anchors healthy resting levels low,
# routes MSTN (high-range catabolic driver) into Catabolism, uses balance/difference formulations.
# Roughly doubles dynamic range (Sarcopenia 0.17->0.31) with healthy < aged < cachectic ordering.
_u, _d = N['FOXO3']; N['FOXO3'] = ([(0.12, 1)] + [(w, e) for w, e in _u if e != 1], _d)  # drop high up-basal 0.5
N['TRIM63'] = ([(0.08, 1), (1, 'NFkB'), (1, 'FOXO3')], [(1.0, 1)])
N['UPS'] = ([(0.08, 1), (1, 'NFkB'), (1, 'TRIM63'), (1, 'smad2_3'), (1.5, 'glucocorticoids')], [(1.3, 1), (1, 'IL15')])
N['Mito_dysfunction'] = ([(0.12, 1), (1, 'ROS'), (1.5, 'hyperammonemia')], [(0.5, 1), (2.5, 'PGC1A')])
N['Proteolysis'] = ([(0.05, 1), (1.5, 'UPS'), (0.5, 'Autophagy'), (1.2, 'FOXO3')], [(2.0, 1), (1.0, 'FNDC5')])
N['Anabolism'] = ([(0.3, 1), (2.4, 'Protein_synthesis')], [(0.3, 1), (0.8, 'Proteolysis')])
N['Catabolism'] = ([(0.05, 1), (1.2, 'Proteolysis'), (1.3, 'MSTN'), (0.9, 'Mito_dysfunction')], [(0.3, 1), (1.6, 'Anabolism')])
N['Sarcopenia'] = ([(0.05, 1), (3.4, 'Catabolism')], [(0.05, 1), (3.4, 'Anabolism')])

# ================= AGING LAYER =================
# Hallmarks-of-aging operationalized as readout sub-scores over mechanistic nodes already in the map
# (López-Otín 2013 PMID 23746838 / 2023 PMID 36599349). Higher = more aged. Built from pathway activity
# and secretory markers, NOT from epigenetic clocks — so clocks stay an independent, external comparison
# (no circularity). Each sub-score is a pure readout; nothing downstream depends on them.
# Directional anchors (PubMed-verified): NAD decline via CD38 — Camacho-Pereira 2016 (27304511);
# inflammaging IL-6+TNF mortality index — Varadhan 2014 (23689826); GDF15 all-cause mortality — Wiklund 2010 (20854422).
N['Aging_nutrient_sensing']   = ([(0.05, 1), (1.0, 'mTOR'), (0.7, 'IGF1'), (0.8, 'Insulin_resistance')], [(0.3, 1), (1.2, 'AMPK')])
N['Aging_mito']               = ([(0.05, 1), (1.2, 'Mito_dysfunction'), (0.8, 'ROS')], [(0.3, 1), (1.3, 'Mito_biogenesis')])
N['Aging_NAD']                = ([(0.05, 1), (1.3, 'CD38')], [(0.3, 1), (1.4, 'NAD')])
N['Aging_inflammation']       = ([(0.05, 1), (1.0, 'NFkB'), (1.0, 'IL6'), (1.0, 'TNF')], [(0.5, 1)])
N['Aging_senescence_genomic'] = ([(0.05, 1), (1.0, 'p53'), (1.0, 'PARP1'), (1.0, 'DNA_damage')], [(0.5, 1)])
N['Aging_proteostasis']       = ([(0.05, 1), (1.3, 'Proteolysis')], [(0.3, 1), (1.0, 'Autophagy'), (0.6, 'Protein_synthesis')])
N['Aging_secretory']          = ([(0.05, 1), (1.3, 'GDF15'), (0.8, 'MSTN')], [(0.5, 1)])
_SUB = ['Aging_nutrient_sensing', 'Aging_mito', 'Aging_NAD', 'Aging_inflammation',
        'Aging_senescence_genomic', 'Aging_proteostasis', 'Aging_secretory']
N['Aging_index'] = ([(0.05, 1)] + [(0.7, s) for s in _SUB], [(2.3, 1)])

READOUTS = ['Aging_index'] + _SUB + ['Sarcopenia', 'Catabolism', 'Anabolism', 'Proteolysis', 'Protein_synthesis',
            'MSTN', 'mTOR', 'AMPK', 'FOXO3', 'Autophagy', 'PGC1A', 'Insulin_sensitivity', 'NAD', 'FGF21',
            'Insulin_resistance', 'TNF', 'IL6']


def write(prefix='MUS_P1_coupled'):
    lines = []
    for n in INPUTS:
        lines.append(f'Node {n} {{ logic = {n}; rate_up = 0.0; rate_down = 0.0; }}')
    for n, (up, dn) in N.items():
        if n in INPUTS:
            continue
        lines.append(f'Node {n} {{ rate_up = {" + ".join(P1.term(w, e) for w, e in up)}; '
                     f'rate_down = {" + ".join(P1.term(w, e) for w, e in dn)}; }}')
    open(prefix + '.bnd', 'w').write('\n'.join(lines) + '\n')
    cfg = ['time_tick = 0.5;', 'max_time = 40;', 'sample_count = 500;', 'discrete_time = 0;',
           'use_physrandgen = 0;', 'seed_pseudorandom = 7;', 'thread_count = 8;', 'statdist_traj_count = 0;']
    for n in INPUTS + [x for x in N if x not in INPUTS]:
        cfg.append(f'{n}.istate = 0;')
    for n in N:
        if n not in INPUTS:
            cfg.append(f'{n}.is_internal = 0;')
    open(prefix + '.cfg', 'w').write('\n'.join(cfg) + '\n')
    return prefix, len(INPUTS) + len([x for x in N if x not in INPUTS])


if __name__ == '__main__':
    pre, ncount = write()
    print(ncount, 'nodes;', len(INPUTS), 'inputs')
    json.dump({'inputs': INPUTS, 'nodes': [x for x in N if x not in INPUTS], 'readouts': READOUTS},
              open('coupled_nodes.json', 'w'), indent=1)
