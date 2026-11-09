"""Methylation input layer (prototype v0).
Maps standard blood-methylation features (EPIC/450K-derived) to model input activation fractions,
so a single patient's methylation profile parameterizes the coupled MUS-P1 model and yields a
personalized aging sub-score profile + sarcopenia readout. Clocks are NOT used as inputs (they
remain external readouts) -> no circularity.

Feature -> node mapping (each feature is a 0-1 percentile within a reference cohort):
  inflammation   (CRP + IL-6 EpiScores, neutrophil:lymphocyte, senescent-CD8 fraction) -> Inflammaging
  immune_aging   (CD4:CD8 ratio, naive:memory, EpiTOC/immune drift)                    -> Inflammaging (adds)
  hepatic        (ALT/GGT/liver EpiScores; ammonia proxy)                              -> hyperammonemia, FGL1
  adiposity_lep  (leptin / BMI / waist EpiScores)                                      -> LEP, Nutrients
  adiponectin    (adiponectin EpiScore)                                                -> ADIPOQ
  damage_smoke   (smoking EpiScore, genotoxic exposure)                                -> DNA_damage
  nad_depletion  (NAD-pathway / CD38 surrogate, if available)                          -> Inflammaging (CD38 arm)
TFMethyl TF-binding-site CpG scores (Patel/Hoffmann) can later clamp internal TF nodes
(NF-kB, FOXO3, p53, PPARG) directly; v0 parameterizes the robust state inputs only.
"""
import cmaboss_128n as cmaboss, json, tempfile, shutil
import coupled_model as CM

INPUTS = CM.INPUTS; NODES = [x for x in CM.N if x not in INPUTS]
SUB = ['Aging_nutrient_sensing', 'Aging_mito', 'Aging_NAD', 'Aging_inflammation',
       'Aging_senescence_genomic', 'Aging_proteostasis', 'Aging_secretory']
READ = ['Aging_index'] + SUB + ['Sarcopenia', 'Anabolism', 'Catabolism', 'MSTN',
                                'Insulin_resistance', 'IL6', 'NAD']
BND = '\n'.join([f'Node {n} {{ logic={n}; rate_up=0.0; rate_down=0.0; }}' for n in INPUTS] +
                [f'Node {n} {{ rate_up = {" + ".join(CM.P1.term(w, e) for w, e in up)}; '
                 f'rate_down = {" + ".join(CM.P1.term(w, e) for w, e in dn)}; }}'
                 for n, (up, dn) in CM.N.items() if n not in INPUTS])

# feature -> {input: weight}; a feature's percentile is distributed to inputs by these weights
FEATURE_MAP = {
    'inflammation':  {'Inflammaging': 1.0},
    'immune_aging':  {'Inflammaging': 0.5},
    'hepatic':       {'hyperammonemia': 0.9, 'FGL1': 0.8},
    'adiposity_lep': {'LEP': 1.0, 'Nutrients': 0.6},
    'adiponectin':   {'ADIPOQ': 1.0},
    'damage_smoke':  {'DNA_damage': 1.0},
    'nad_depletion': {'Inflammaging': 0.3},   # proxy: drives CD38 arm -> NAD decline
}


def features_to_fractions(feat):
    """feat: {feature: percentile 0-1}. Returns {input: activation fraction 0-1}.

    A feature that is missing / None / NaN is skipped (it does NOT drive its inputs),
    and an input with no present contributor is left unset -> OFF (never a constant
    mid-level drive). Each input is a WEIGHTED blend of its present contributors,
    normalized by the summed weight ONLY when that would exceed full drive: this keeps
    the intended (possibly partial) drive of a single low-weight feature, while
    preventing several features from saturating a shared input (e.g. Inflammaging)."""
    num, wsum = {}, {}
    for f, val in feat.items():
        if val is None:
            continue
        try:
            v = float(val)
        except (TypeError, ValueError):
            continue
        if v != v:                                   # NaN
            continue
        for inp, w in FEATURE_MAP.get(f, {}).items():
            num[inp] = num.get(inp, 0.0) + w * v
            wsum[inp] = wsum.get(inp, 0.0) + w
    return {inp: min(1.0, num[inp] / max(1.0, wsum[inp])) for inp in num}


def run_patient(feat, extra_inputs=None, sc=1500):
    frac = features_to_fractions(feat)
    fix = dict(extra_inputs or {})                 # interventions as full ON
    d = tempfile.mkdtemp(); open(d + '/m.bnd', 'w').write(BND)
    cfg = ['time_tick=0.5;', 'max_time=40;', f'sample_count={sc};', 'discrete_time=0;', 'use_physrandgen=0;',
           'seed_pseudorandom=7;', 'thread_count=8;', 'statdist_traj_count=0;']
    for n in INPUTS:
        if n in fix:
            cfg.append(f'{n}.istate={fix[n]};')
        elif n in frac:
            p = round(frac[n], 4); cfg.append(f'[{n}].istate = {p}[1],{round(1-p,4)}[0];')
        else:
            cfg.append(f'{n}.istate=0;')
    for n in NODES: cfg.append(f'{n}.istate=0;')
    for n in NODES: cfg.append(f'{n}.is_internal=0;')
    open(d + '/m.cfg', 'w').write('\n'.join(cfg))
    a, t, nodes, _ = cmaboss.MaBoSSSim(network=d + '/m.bnd', config=d + '/m.cfg').run().get_last_nodes_probtraj()
    shutil.rmtree(d, ignore_errors=True)
    L = {n: float(a[0][i]) for i, n in enumerate(nodes)}
    return {k: round(L.get(k, 0), 3) for k in READ}


# --- synthetic demonstration patients (percentiles 0-1) ---
PATIENTS = {
 'P1 healthy (60y, fit)':        {'inflammation': 0.2, 'immune_aging': 0.25, 'hepatic': 0.1,
                                  'adiposity_lep': 0.25, 'adiponectin': 0.7, 'damage_smoke': 0.1, 'nad_depletion': 0.25},
 'P2 inflammaged (75y)':         {'inflammation': 0.85, 'immune_aging': 0.8, 'hepatic': 0.3,
                                  'adiposity_lep': 0.5, 'adiponectin': 0.3, 'damage_smoke': 0.4, 'nad_depletion': 0.8},
 'P3 metabolic/hepatic (65y)':   {'inflammation': 0.55, 'immune_aging': 0.4, 'hepatic': 0.8,
                                  'adiposity_lep': 0.85, 'adiponectin': 0.15, 'damage_smoke': 0.3, 'nad_depletion': 0.6},
 'P4 ex-smoker, lean (70y)':     {'inflammation': 0.45, 'immune_aging': 0.6, 'hepatic': 0.2,
                                  'adiposity_lep': 0.3, 'adiponectin': 0.55, 'damage_smoke': 0.85, 'nad_depletion': 0.5},
}
IVS = {'Exercise+Leucine': {'Exercise': 1, 'leucine': 1}, 'NMN': {'NMN_supp': 1}, 'Omega-3': {'Omega3_PUFA': 1},
       'Metformin': {'Metformin': 1}, 'Rapamycin': {'Rapamycin': 1}, 'DR': {'DR': 1}}

if __name__ == "__main__":
    out = {'feature_map': FEATURE_MAP, 'patients': {}}
    for name, feat in PATIENTS.items():
        base = run_patient(feat)
        rec = {'features': feat, 'input_fractions': {k: round(v, 3) for k, v in features_to_fractions(feat).items()},
               'baseline': base, 'interventions': {}}
        for ivn, iv in IVS.items():
            r = run_patient(feat, iv)
            rec['interventions'][ivn] = {k: round(r[k] - base[k], 3) for k in READ}
        out['patients'][name] = rec
    json.dump(out, open('methylation_layer.json', 'w'), indent=1)

    # --- personalized recommendation: target the hallmark most ELEVATED vs a healthy reference ---
    # (absolute sub-score levels still carry resting-level offset; deviation-from-reference is the
    #  clinically meaningful signal and the right operationalization until the sub-scores are re-anchored)
    REF = out['patients']['P1 healthy (60y, fit)']['baseline']
    for name, rec in out['patients'].items():
        b = rec['baseline']
        devs = {s: round(b[s] - REF[s], 3) for s in SUB}
        rec['deviation_vs_healthy'] = devs
        dom = max(SUB, key=lambda s: devs[s])                     # most-elevated hallmark vs healthy
        best = min(rec['interventions'].items(), key=lambda kv: kv[1][dom])
        rec['dominant_hallmark'] = dom
        rec['recommended'] = {'intervention': best[0], 'delta_on_dominant': best[1][dom],
                              'delta_sarcopenia': best[1]['Sarcopenia'], 'dominant_elevation': devs[dom]}
    json.dump(out, open('methylation_layer.json', 'w'), indent=1)

    # console summary
    cols = ['Aging_index', 'Aging_inflammation', 'Aging_NAD', 'Aging_nutrient_sensing', 'Aging_secretory', 'Sarcopenia']
    print('Patient baseline aging profiles (from methylation features):')
    print('patient'.ljust(30) + ''.join(c.replace('Aging_', '')[:10].rjust(12) for c in cols))
    for name, rec in out['patients'].items():
        b = rec['baseline']
        print(name.ljust(30) + ''.join(f'{b[c]:12.3f}' for c in cols))
    print('\nPersonalized recommendation (target the dominant elevated hallmark):')
    for name, rec in out['patients'].items():
        r = rec['recommended']
        print(f'  {name:30} dominant: {rec["dominant_hallmark"].replace("Aging_",""):18} '
              f'-> {r["intervention"]:16} (Δ {r["delta_on_dominant"]:+.3f} on it, ΔSarcopenia {r["delta_sarcopenia"]:+.3f})')

    # --- figure: patient × sub-score profile heatmap ---
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, numpy as np
    names = list(out['patients'].keys()); short = [s.replace('Aging_', '') for s in SUB]
    M = np.array([[out['patients'][n]['baseline'][s] for s in SUB] for n in names])
    fig, ax = plt.subplots(figsize=(9.5, 3.8))
    im = ax.imshow(M, cmap='Oranges', vmin=0.2, vmax=0.85, aspect='auto')
    ax.set_xticks(range(len(SUB))); ax.set_xticklabels(short, rotation=35, ha='right', fontsize=8)
    ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=8)
    for yi in range(len(names)):
        for xi in range(len(SUB)):
            ax.text(xi, yi, f'{M[yi,xi]:.2f}', ha='center', va='center', fontsize=7,
                    color='white' if M[yi, xi] > 0.6 else '#222')
    ax.set_title('Methylation-driven aging sub-score profile per patient (higher = more aged)', fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label='sub-score')
    fig.tight_layout(); fig.savefig('fig_methylation_patients.png', dpi=170, bbox_inches='tight')
    print('\nsaved fig_methylation_patients.png')
