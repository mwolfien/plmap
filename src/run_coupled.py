"""Run lean-mass trade-off scenarios on the coupled MUS-P1 model (cmaboss_128n)."""
import cmaboss_128n as cmaboss
import json, tempfile, os, shutil
import coupled_model as CM

PRE, NC = CM.write()
BND, CFG = PRE + '.bnd', PRE + '.cfg'
RO = CM.READOUTS

SCEN = {
 'baseline':            {},
 'DR':                  {'DR': 1},
 'GLP1 (semaglutide)':  {'Semaglutide': 1},
 'Rapamycin':           {'Rapamycin': 1},
 'Exercise':            {'Exercise': 1},
 'Leucine':             {'leucine': 1},
 'Exercise+Leucine':    {'Exercise': 1, 'leucine': 1},
 'DR+Exercise':         {'DR': 1, 'Exercise': 1},
 'Rapamycin+Exercise':  {'Rapamycin': 1, 'Exercise': 1},
 'GLP1+Exercise':       {'Semaglutide': 1, 'Exercise': 1},
 'Inflammaging':        {'Inflammaging': 1},
 'Hyperammonemia(liver)': {'hyperammonemia': 1},
 'Hyperammonemia+Exercise': {'hyperammonemia': 1, 'Exercise': 1},
 'Aged (inflammaging+ammonia)': {'Inflammaging': 1, 'hyperammonemia': 1},
 'Aged + Exercise+Leucine': {'Inflammaging': 1, 'hyperammonemia': 1, 'Exercise': 1, 'leucine': 1},
 'Obese (leptin+FGL1)':  {'LEP': 1, 'FGL1': 1},
 'Obese + Adiponectin':  {'LEP': 1, 'FGL1': 1, 'ADIPOQ': 1},
 'Obese + Omega-3':      {'LEP': 1, 'FGL1': 1, 'Omega3_PUFA': 1},
 'Obese + Exercise':     {'LEP': 1, 'FGL1': 1, 'Exercise': 1},
}


def run(fix):
    d = tempfile.mkdtemp()
    b, c = d + '/m.bnd', d + '/m.cfg'
    shutil.copy(BND, b)
    cfg = open(CFG).read()
    for k, v in fix.items():
        cfg = cfg.replace(f'{k}.istate = 0;', f'{k}.istate = {v};')
    open(c, 'w').write(cfg)
    sim = cmaboss.MaBoSSSim(network=b, config=c)
    res = sim.run()
    arr, times, nodes, _ = res.get_last_nodes_probtraj()
    last = {n: float(arr[0][i]) for i, n in enumerate(nodes)}
    shutil.rmtree(d, ignore_errors=True)
    return {n: round(last.get(n, 0.0), 3) for n in RO}


out = {name: run(fix) for name, fix in SCEN.items()}
json.dump(out, open('coupled_scenarios.json', 'w'), indent=1)

print(f"{'scenario':30} " + ' '.join(f'{r[:9]:>9}' for r in ['Sarcopen', 'Catabol', 'Anabol', 'MSTN', 'ProtSyn', 'InsResist', 'IL6', 'AMPK']))
cols = ['Sarcopenia', 'Catabolism', 'Anabolism', 'MSTN', 'Protein_synthesis', 'Insulin_resistance', 'IL6', 'AMPK']
for name, r in out.items():
    print(f"{name:30} " + ' '.join(f'{r[c]:9.3f}' for c in cols))
