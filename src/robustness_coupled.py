"""Robustness of the coupled MUS-P1 qualitative conclusions to weight perturbation.
Each rep multiplies every edge weight by a random factor in [1-eps, 1+eps] and rechecks
whether the key intervention contrasts keep their sign. Reports the fraction preserved."""
import cmaboss_128n as cmaboss
import json, tempfile, shutil, random, copy
import coupled_model as CM

random.seed(1)
EPS = 0.30
REPS = 16
INPUTS = CM.INPUTS
RO = ['Sarcopenia', 'Catabolism', 'Anabolism', 'MSTN', 'Protein_synthesis', 'Autophagy', 'Insulin_resistance', 'IL6']

SC = {'baseline': {}, 'DR': {'DR': 1}, 'Leucine': {'leucine': 1}, 'Exercise': {'Exercise': 1},
      'Rapamycin': {'Rapamycin': 1}, 'Hyperammonemia': {'hyperammonemia': 1},
      'Aged': {'Inflammaging': 1, 'hyperammonemia': 1},
      'Aged+Ex+Leu': {'Inflammaging': 1, 'hyperammonemia': 1, 'Exercise': 1, 'leucine': 1},
      'Obese': {'LEP': 1, 'FGL1': 1}, 'Obese+ADIPOQ': {'LEP': 1, 'FGL1': 1, 'ADIPOQ': 1},
      'Obese+O3': {'LEP': 1, 'FGL1': 1, 'Omega3_PUFA': 1}}

# contrast: (name, scen_a, scen_b, readout, expected sign of a-b)  [+1 = a>b]
CONTRASTS = [
 ('Leucine raises Anabolism', 'Leucine', 'baseline', 'Anabolism', +1),
 ('DR lowers Anabolism (mTOR trade-off)', 'DR', 'baseline', 'Anabolism', -1),
 ('Rapamycin lowers Anabolism', 'Rapamycin', 'baseline', 'Anabolism', -1),
 ('Exercise lowers MSTN', 'Exercise', 'baseline', 'MSTN', -1),
 ('Hyperammonemia raises MSTN (liver bridge)', 'Hyperammonemia', 'baseline', 'MSTN', +1),
 ('Aged+Ex+Leu rescues Anabolism vs Aged', 'Aged+Ex+Leu', 'Aged', 'Anabolism', +1),
 ('Aged+Ex+Leu lowers Sarcopenia vs Aged', 'Aged+Ex+Leu', 'Aged', 'Sarcopenia', -1),
 ('Aged+Ex+Leu lowers MSTN vs Aged', 'Aged+Ex+Leu', 'Aged', 'MSTN', -1),
 ('Obesity raises insulin resistance', 'Obese', 'baseline', 'Insulin_resistance', +1),
 ('Adiponectin lowers insulin resistance (sign fix)', 'Obese+ADIPOQ', 'Obese', 'Insulin_resistance', -1),
 ('Omega-3 lowers IL-6', 'Obese+O3', 'Obese', 'IL6', -1),
]


def perturb(Nd, eps):
    out = {}
    for n, (up, dn) in Nd.items():
        f = lambda terms: [(round(w * random.uniform(1 - eps, 1 + eps), 4), e) for w, e in terms]
        out[n] = (f(up), f(dn))
    return out


def write_tmp(Nd, d):
    lines = [f'Node {n} {{ logic = {n}; rate_up = 0.0; rate_down = 0.0; }}' for n in INPUTS]
    for n, (up, dn) in Nd.items():
        if n in INPUTS:
            continue
        lines.append(f'Node {n} {{ rate_up = {" + ".join(CM.P1.term(w, e) for w, e in up)}; '
                     f'rate_down = {" + ".join(CM.P1.term(w, e) for w, e in dn)}; }}')
    open(d + '/m.bnd', 'w').write('\n'.join(lines) + '\n')
    cfg = ['time_tick=0.5;', 'max_time=40;', 'sample_count=300;', 'discrete_time=0;', 'use_physrandgen=0;',
           'seed_pseudorandom=7;', 'thread_count=8;', 'statdist_traj_count=0;']
    alln = INPUTS + [x for x in Nd if x not in INPUTS]
    for n in alln:
        cfg.append(f'{n}.istate=0;')
    for n in Nd:
        if n not in INPUTS:
            cfg.append(f'{n}.is_internal=0;')
    open(d + '/base.cfg', 'w').write('\n'.join(cfg))


def run(d, fix):
    cfg = open(d + '/base.cfg').read()
    for k, v in fix.items():
        cfg = cfg.replace(f'{k}.istate=0;', f'{k}.istate={v};')
    open(d + '/r.cfg', 'w').write(cfg)
    arr, times, nodes, _ = cmaboss.MaBoSSSim(network=d + '/m.bnd', config=d + '/r.cfg').run().get_last_nodes_probtraj()
    return {n: float(arr[0][i]) for i, n in enumerate(nodes)}


def eval_model(Nd):
    d = tempfile.mkdtemp()
    write_tmp(Nd, d)
    vals = {s: run(d, fix) for s, fix in SC.items()}
    shutil.rmtree(d, ignore_errors=True)
    res = {}
    for name, a, b, ro, sign in CONTRASTS:
        diff = vals[a].get(ro, 0) - vals[b].get(ro, 0)
        res[name] = (1 if diff > 0 else -1 if diff < 0 else 0, diff)
    return res


base = eval_model(CM.N)
counts = {c[0]: 0 for c in CONTRASTS}
mags = {c[0]: [] for c in CONTRASTS}
for r in range(REPS):
    res = eval_model(perturb(CM.N, EPS))
    for name, a, b, ro, sign in CONTRASTS:
        s, diff = res[name]
        if s == sign:
            counts[name] += 1
        mags[name].append(diff)

print(f'Robustness: {REPS} reps, +-{int(EPS*100)}% weight noise\n')
print(f"{'contrast':46} {'base':>7} {'preserved':>10} {'mean d':>8}")
summary = {}
for name, a, b, ro, sign in CONTRASTS:
    bs, bd = base[name]
    frac = counts[name] / REPS
    md = sum(mags[name]) / len(mags[name])
    ok = 'OK' if bs == sign else 'FLIP'
    print(f"{name:46} {bd:+7.3f} {frac*100:8.0f}% {md:+8.3f}  {ok}")
    summary[name] = {'baseline_diff': round(bd, 3), 'expected_sign': sign, 'frac_preserved': round(frac, 3),
                     'mean_diff': round(md, 3)}
json.dump({'eps': EPS, 'reps': REPS, 'contrasts': summary}, open('coupled_robustness.json', 'w'), indent=1)
