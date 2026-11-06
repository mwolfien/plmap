"""Systematic characterization of the coupled MUS-P1 model: isolated input sweeps, dose-response,
pairwise interactions, patient archetypes, and a dynamic-range / dead-lever boundary scan.
Writes coupled_characterization.json. Uses cmaboss_128n."""
import cmaboss_128n as cmaboss
import json, tempfile, shutil, itertools
import coupled_model as CM

INPUTS = CM.INPUTS
NODES = [x for x in CM.N if x not in INPUTS]
READOUTS = ['Aging_index', 'Aging_nutrient_sensing', 'Aging_mito', 'Aging_NAD', 'Aging_inflammation',
            'Aging_senescence_genomic', 'Aging_proteostasis', 'Aging_secretory',
            'Sarcopenia', 'Anabolism', 'Catabolism', 'Proteolysis', 'Protein_synthesis', 'MSTN',
            'mTOR', 'AMPK', 'FOXO3', 'Autophagy', 'PGC1A', 'Mito_biogenesis', 'NAD', 'SIRT1',
            'Insulin_sensitivity', 'Insulin_resistance', 'IL6', 'TNF', 'FGF21', 'ROS', 'NFkB']
BND_LINES = ([f'Node {n} {{ logic = {n}; rate_up = 0.0; rate_down = 0.0; }}' for n in INPUTS] +
             [f'Node {n} {{ rate_up = {" + ".join(CM.P1.term(w, e) for w, e in up)}; '
              f'rate_down = {" + ".join(CM.P1.term(w, e) for w, e in dn)}; }}'
              for n, (up, dn) in CM.N.items() if n not in INPUTS])
BND = '\n'.join(BND_LINES)


def run(fix, frac=None, sc=600):
    """fix: {input:0/1}; frac: {input:p in [0,1]} fractional activation."""
    frac = frac or {}
    d = tempfile.mkdtemp()
    open(d + '/m.bnd', 'w').write(BND)
    cfg = ['time_tick=0.5;', 'max_time=40;', f'sample_count={sc};', 'discrete_time=0;', 'use_physrandgen=0;',
           'seed_pseudorandom=7;', 'thread_count=8;', 'statdist_traj_count=0;']
    for n in INPUTS:
        if n in frac:
            p = frac[n]; cfg.append(f'[{n}].istate = {p}[1],{round(1-p,4)}[0];')
        else:
            cfg.append(f'{n}.istate = {fix.get(n,0)};')
    for n in NODES:
        cfg.append(f'{n}.istate = 0;')
    for n in NODES:
        cfg.append(f'{n}.is_internal = 0;')
    open(d + '/m.cfg', 'w').write('\n'.join(cfg))
    arr, t, nodes, _ = cmaboss.MaBoSSSim(network=d + '/m.bnd', config=d + '/m.cfg').run().get_last_nodes_probtraj()
    shutil.rmtree(d, ignore_errors=True)
    last = {n: float(arr[0][i]) for i, n in enumerate(nodes)}
    return {n: round(last.get(n, 0.0), 3) for n in READOUTS}


OUT = {}
base = run({})
OUT['baseline'] = base

# ---- A. single-input isolation: each input ON alone ----
A = {}
for inp in INPUTS:
    r = run({inp: 1})
    A[inp] = {k: round(r[k] - base[k], 3) for k in READOUTS}
OUT['single_input_delta'] = A

# ---- B. dose-response (fractional activation) for key levers ----
levers = ['DR', 'Rapamycin', 'Exercise', 'leucine', 'Inflammaging', 'NMN_supp', 'Metformin', 'Semaglutide']
B = {}
for lv in levers:
    B[lv] = {}
    for p in [0.0, 0.25, 0.5, 0.75, 1.0]:
        r = run({}, frac={lv: p}) if p > 0 else base
        B[lv][p] = {k: r[k] for k in ['Sarcopenia', 'Anabolism', 'MSTN', 'Autophagy', 'NAD',
                                      'Insulin_sensitivity', 'Mito_biogenesis', 'Protein_synthesis']}
OUT['dose_response'] = B

# ---- C. pairwise interaction (synergy/antagonism) ----
IVs = ['DR', 'Rapamycin', 'Exercise', 'leucine', 'Metformin', 'NMN_supp', 'Semaglutide', 'Omega3_PUFA', 'Resveratrol']
key = ['Sarcopenia', 'Anabolism', 'Autophagy', 'NAD', 'Insulin_sensitivity', 'Mito_biogenesis']
single = {iv: {k: run({iv: 1})[k] - base[k] for k in key} for iv in IVs}
C = {}
for a, b in itertools.combinations(IVs, 2):
    comb = run({a: 1, b: 1})
    C[f'{a}+{b}'] = {k: {'combined_delta': round(comb[k] - base[k], 3),
                         'additive_pred': round(single[a][k] + single[b][k], 3),
                         'interaction': round((comb[k] - base[k]) - (single[a][k] + single[b][k]), 3)}
                     for k in key}
OUT['single_delta_key'] = {iv: {k: round(v, 3) for k, v in d.items()} for iv, d in single.items()}
OUT['pairwise_interaction'] = C

# ---- D. patient archetypes + interventions ----
ARCH = {
 'Healthy young': {},
 'Metabolically unhealthy (obese)': {'LEP': 1, 'FGL1': 1, 'Nutrients': 1},
 'Inflammaged': {'Inflammaging': 1, 'DNA_damage': 1},
 'Cirrhotic / hepatic (Schneider)': {'hyperammonemia': 1, 'FGL1': 1, 'Inflammaging': 1},
 'Sarcopenic-aged': {'Inflammaging': 1, 'hyperammonemia': 1, 'DNA_damage': 1},
}
IV_SET = {'(none)': {}, 'Exercise': {'Exercise': 1}, 'Exercise+Leucine': {'Exercise': 1, 'leucine': 1},
          'DR': {'DR': 1}, 'Rapamycin': {'Rapamycin': 1}, 'NMN': {'NMN_supp': 1},
          'Metformin': {'Metformin': 1}, 'Omega-3': {'Omega3_PUFA': 1}, 'GLP-1': {'Semaglutide': 1}}
D = {}
for an, astate in ARCH.items():
    D[an] = {}
    for ivn, iv in IV_SET.items():
        fix = dict(astate); fix.update(iv)
        r = run(fix)
        D[an][ivn] = {k: r[k] for k in ['Sarcopenia', 'Anabolism', 'MSTN', 'Insulin_resistance',
                                        'IL6', 'NAD', 'Autophagy', 'Mito_biogenesis']}
OUT['archetypes'] = D

# ---- E. boundary scan: dynamic range per readout; dead levers ----
rng = {}
for k in READOUTS:
    vals = [base[k]] + [base[k] + A[inp][k] for inp in INPUTS]
    rng[k] = round(max(vals) - min(vals), 3)
dead = {inp: round(sum(abs(A[inp][k]) for k in READOUTS), 3) for inp in INPUTS}
OUT['readout_dynamic_range'] = dict(sorted(rng.items(), key=lambda x: -x[1]))
OUT['input_total_influence'] = dict(sorted(dead.items(), key=lambda x: -x[1]))

json.dump(OUT, open('coupled_characterization.json', 'w'), indent=1)

# compact console summary
print('baseline anchored. readouts:', len(READOUTS), 'inputs:', len(INPUTS))
print('\nReadout dynamic range (max-min across single-input sweeps):')
for k, v in OUT['readout_dynamic_range'].items():
    print(f'  {k:22} {v:.3f}' + ('   <-- flat / model cannot move' if v < 0.05 else ''))
print('\nInput total influence (sum|delta| across readouts):')
for inp, v in OUT['input_total_influence'].items():
    print(f'  {inp:16} {v:.3f}' + ('   <-- dead lever' if v < 0.05 else ''))
print('\nStrongest single-input effects (|delta|>=0.15):')
for inp in INPUTS:
    big = [(k, A[inp][k]) for k in READOUTS if abs(A[inp][k]) >= 0.15]
    if big:
        print(f'  {inp:16} ' + ', '.join(f'{k}{v:+.2f}' for k, v in big))
print('\nNotable interactions (|interaction|>=0.05) on Sarcopenia/Anabolism/NAD:')
for pair, d in C.items():
    for k in ['Sarcopenia', 'Anabolism', 'NAD']:
        if abs(d[k]['interaction']) >= 0.05:
            tag = 'synergy' if d[k]['interaction'] * (1 if k == 'Anabolism' or k == 'NAD' else -1) > 0 else 'antagonism/redundant'
            print(f'  {pair:26} {k:12} interaction {d[k]["interaction"]:+.3f} ({tag})')
