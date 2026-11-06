"""Aging layer analysis: hallmark subscores across patient archetypes and interventions.
Writes aging_layer.json and two figures."""
import cmaboss_128n as cmaboss, json, tempfile, shutil
import coupled_model as CM
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np

INPUTS = CM.INPUTS; NODES = [x for x in CM.N if x not in INPUTS]
SUB = ['Aging_nutrient_sensing', 'Aging_mito', 'Aging_NAD', 'Aging_inflammation',
       'Aging_senescence_genomic', 'Aging_proteostasis', 'Aging_secretory']
SHORT = ['nutrient-sensing', 'mitochondrial', 'NAD decline', 'inflammaging',
         'senescence/genomic', 'proteostasis', 'secretory (GDF15)']
READ = ['Aging_index'] + SUB
BND = '\n'.join([f'Node {n} {{ logic={n}; rate_up=0.0; rate_down=0.0; }}' for n in INPUTS] +
                [f'Node {n} {{ rate_up = {" + ".join(CM.P1.term(w, e) for w, e in up)}; '
                 f'rate_down = {" + ".join(CM.P1.term(w, e) for w, e in dn)}; }}'
                 for n, (up, dn) in CM.N.items() if n not in INPUTS])


def run(fix, sc=800):
    d = tempfile.mkdtemp(); open(d + '/m.bnd', 'w').write(BND)
    cfg = ['time_tick=0.5;', 'max_time=40;', f'sample_count={sc};', 'discrete_time=0;', 'use_physrandgen=0;',
           'seed_pseudorandom=7;', 'thread_count=8;', 'statdist_traj_count=0;']
    for n in INPUTS: cfg.append(f'{n}.istate={fix.get(n,0)};')
    for n in NODES: cfg.append(f'{n}.istate=0;')
    for n in NODES: cfg.append(f'{n}.is_internal=0;')
    open(d + '/m.cfg', 'w').write('\n'.join(cfg))
    a, t, nodes, _ = cmaboss.MaBoSSSim(network=d + '/m.bnd', config=d + '/m.cfg').run().get_last_nodes_probtraj()
    shutil.rmtree(d, ignore_errors=True)
    L = {n: float(a[0][i]) for i, n in enumerate(nodes)}
    return {k: round(L.get(k, 0), 3) for k in READ}

ARCH = {'Healthy (fed)': {'Nutrients': 1}, 'Obese': {'LEP': 1, 'FGL1': 1, 'Nutrients': 1},
        'Inflammaged': {'Inflammaging': 1, 'DNA_damage': 1},
        'Cirrhotic/hepatic': {'hyperammonemia': 1, 'FGL1': 1, 'Inflammaging': 1},
        'Sarcopenic-aged': {'Inflammaging': 1, 'hyperammonemia': 1, 'DNA_damage': 1}}
IVS = {'(none)': {}, 'Exercise': {'Exercise': 1}, 'Ex+Leucine': {'Exercise': 1, 'leucine': 1},
       'DR': {'DR': 1}, 'Rapamycin': {'Rapamycin': 1}, 'NMN': {'NMN_supp': 1}, 'NR': {'NR_supp': 1},
       'Metformin': {'Metformin': 1}, 'Resveratrol': {'Resveratrol': 1}, 'Omega-3': {'Omega3_PUFA': 1},
       'CD38i': {'CD38i': 1}, 'GLP-1': {'Semaglutide': 1}}

# intervention effect on subscores (applied on healthy-fed background)
base = run({'Nutrients': 1})
IV = {}
for ivn, iv in IVS.items():
    if ivn == '(none)': continue
    fix = {'Nutrients': 1}; fix.update(iv)
    r = run(fix); IV[ivn] = {k: round(r[k] - base[k], 3) for k in READ}
# archetype profiles
ARCHP = {an: run(st) for an, st in ARCH.items()}
OUT = {'baseline_healthy_fed': base, 'intervention_delta': IV, 'archetype_profile': ARCHP, 'subscores': SUB}
json.dump(OUT, open('aging_layer.json', 'w'), indent=1)

# ---- figure 1: intervention × subscore delta heatmap ----
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9})
ivs = list(IV.keys())
M = np.array([[IV[iv][s] for s in SUB] for iv in ivs])
fig, ax = plt.subplots(figsize=(9, 6.6))
im = ax.imshow(M, cmap='RdBu_r', vmin=-0.4, vmax=0.4, aspect='auto')
ax.set_xticks(range(len(SUB))); ax.set_xticklabels(SHORT, rotation=40, ha='right', fontsize=8)
ax.set_yticks(range(len(ivs))); ax.set_yticklabels(ivs, fontsize=9)
for yi in range(len(ivs)):
    for xi in range(len(SUB)):
        v = M[yi, xi]
        if abs(v) >= 0.04:
            ax.text(xi, yi, f'{v:+.2f}', ha='center', va='center', fontsize=6.5,
                    color='white' if abs(v) > 0.25 else '#222')
ax.set_title('Which hallmark each intervention targets\n(Δ aging sub-score vs healthy-fed baseline; blue = lowers aging)', fontsize=10)
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label='Δ sub-score')
fig.tight_layout(); fig.savefig('fig_aging_interventions.png', dpi=170, bbox_inches='tight')
print('saved fig_aging_interventions.png')

# ---- figure 2: archetype × subscore profile heatmap ----
ans = list(ARCHP.keys())
M2 = np.array([[ARCHP[an][s] for s in SUB] for an in ans])
fig2, ax2 = plt.subplots(figsize=(9, 4.2))
im2 = ax2.imshow(M2, cmap='Oranges', vmin=0, vmax=0.9, aspect='auto')
ax2.set_xticks(range(len(SUB))); ax2.set_xticklabels(SHORT, rotation=40, ha='right', fontsize=8)
ax2.set_yticks(range(len(ans))); ax2.set_yticklabels(ans, fontsize=9)
for yi in range(len(ans)):
    for xi in range(len(SUB)):
        ax2.text(xi, yi, f'{M2[yi,xi]:.2f}', ha='center', va='center', fontsize=7,
                 color='white' if M2[yi, xi] > 0.55 else '#222')
ax2.set_title('Aging sub-score profile per patient archetype (higher = more aged on that hallmark)', fontsize=10)
fig2.colorbar(im2, ax=ax2, fraction=0.046, pad=0.02, label='sub-score')
fig2.tight_layout(); fig2.savefig('fig_aging_archetypes.png', dpi=170, bbox_inches='tight')
print('saved fig_aging_archetypes.png')

print('\nIndex: healthy-fed', base['Aging_index'], '| archetypes:',
      {a: ARCHP[a]['Aging_index'] for a in ARCHP})
print('\nIntervention Δ on Aging_index:', {iv: IV[iv]['Aging_index'] for iv in ivs})
