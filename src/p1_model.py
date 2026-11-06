"""P1 core MaBoSS model: hand-curated gates/weights over the curated P1 edges (+ reused NSCLC edges).
Each node: rate_up = basal_up + sum(weight * activator), rate_down = basal_down + sum(weight * inhibitor).
Stationary probability ~ up/(up+down) gives a graded activity in [0,1]. Comments name the supporting rid."""
INPUTS = ['Nutrients', 'DR', 'Inflammaging', 'DNA_damage', 'Rapamycin', 'Metformin', 'NMN_supp', 'NR_supp',
          'Resveratrol', 'Semaglutide', 'Canagliflozin', 'CD38i', 'NNMT']
N = {  # node: (up terms, down terms)  term = (weight, expr) ; expr is a node or '!node'
 'Insulin':   ([(1, 'Nutrients')], [(0.3, 1), (2, 'DR')]),                        # P1B-031
 'PPARA':     ([(0.1, 1), (1, 'DR')], [(1, 1)]),                                  # fasting proxy (assumption)
 'FGF21':     ([(0.1, 1), (1, 'DR'), (1, 'PPARA')], [(1, 1)]),                    # P1C-001, P1C-002
 'FGF21R':    ([(1, 'FGF21')], [(0.5, 1)]),                                        # P1C-003 KLB-FGFR1
 'IGF1':      ([(1, 'Nutrients')], [(0.3, 1), (1, 'DR'), (1, 'FGF21R')]),         # P1B-032 (mouse; human conflict), P1C-007
 'INSR':      ([(1, 'Insulin')], [(0.5, 1)]),                                      # P1B-001
 'IGF1R':     ([(1, 'IGF1')], [(0.5, 1)]),                                         # P1B-004
 'IRS1':      ([(1, 'INSR'), (1, 'IGF1R')], [(0.3, 1), (1, 'S6K1')]),             # P1B-002, -005, -015
 'PI3K':      ([(1, 'IRS1')], [(0.5, 1)]),                                         # P1B-003
 'mTORC2':    ([(1, 'PI3K')], [(0.3, 1), (2, 'Rapamycin')]),                       # P1B-021 (chronic)
 'Akt':       ([(0.5, 'PI3K'), (0.7, 'mTORC2')], [(0.5, 1)]),                      # M3-050, P1B-019
 'FOXO3':     ([(0.5, 1), (0.5, 'SIRT1')], [(0.2, 1), (2, 'Akt')]),               # P1A-018, P1B-006
 'Complex_I': ([(1, 1)], [(0.1, 1), (2, 'Metformin'), (2, 'Canagliflozin')]),     # P1B-025, P1C-015 (via complex I)
 'LKB1':      ([(1, 1), (1, 'SIRT1')], [(0.5, 1)]),                                # P1A-021
 'AMPK':      ([(0.1, 1), (2, '!Complex_I'), (1, 'DR'), (1, 'FGF21R'), (1, 'Resveratrol')], [(1, 1), (0.5, '!LKB1')]),  # M3-008, P1B-026/029, P1C-004, P1A-031
 'TSC2':      ([(0.3, 1), (1, 'AMPK')], [(0.3, 1), (1, 'Akt')]),                  # M3-006, P1B-008
 'RHEB':      ([(1, 1)], [(0.2, 1), (2, 'TSC2')]),                                 # P1B-009
 'mTOR':      ([(1.5, 'RHEB')], [(0.2, 1), (1, 'AMPK'), (3, 'Rapamycin'), (1, 'DR'), (1, '!Nutrients')]),  # P1B-010, M3-005, P1B-020, P1B-030
 'S6K1':      ([(1, 'mTOR')], [(0.5, 1)]),                                         # P1B-011
 'EIF4EBP1':  ([(1, 1)], [(0.2, 1), (2, 'mTOR')]),                                 # P1B-012 (active = unphosphorylated)
 'Protein_synthesis': ([(1, 'S6K1')], [(0.2, 1), (1, 'EIF4EBP1')]),               # P1B-013/014
 'TFEB':      ([(1, 1)], [(0.2, 1), (2, 'mTOR')]),                                 # P1B-017
 'ULK1':      ([(0.1, 1), (1, 'AMPK'), (1, '!mTOR')], [(0.3, 1), (1, 'mTOR')]),   # M5-053/054
 'Autophagy': ([(0.05, 1), (1, 'ULK1'), (1, 'TFEB'), (1, 'FOXO3'), (0.5, 'DR')], [(1.5, 1)]),  # M5-055, P1B-018, P1B-007, P1B-034
 'PGC1A':     ([(0.1, 1), (1, 'AMPK'), (1, 'SIRT1'), (0.5, 'FGF21R')], [(1, 1)]), # P1B-023, P1A-017, P1C-006
 'Mito_biogenesis': ([(1, 'PGC1A')], [(0.5, 1)]),                                  # P1B-024
 'NAM':       ([(1, 1)], [(0.2, 1), (1, 'NNMT')]),                                 # P1A-012
 'NAMPT':     ([(0.5, 1), (1, 'AMPK')], [(0.5, 1)]),                               # P1A-022
 'NMN':       ([(0.1, 1), (1, 'NAMPT'), (1, 'NR_supp'), (1.5, 'NMN_supp')], [(0.5, 1), (1, '!NAM'), (1.5, 'CD38')]),  # P1A-001/003/005/006, P1A-010
 'NAD':       ([(0.1, 1), (1, 'NMN'), (0.5, 'NR_supp')], [(0.3, 1), (1.5, 'CD38'), (1, 'PARP1')]),  # P1A-004/033, P1A-007/009
 'CD38':      ([(0.1, 1), (2, 'Inflammaging')], [(1, 1), (3, 'CD38i')]),           # P1A-013, P1A-029
 'PARP1':     ([(0.1, 1), (2, 'DNA_damage'), (0.3, 'SIRT6')], [(1, 1)]),           # M2-059 analogue, P1A-026
 'SIRT1':     ([(0.3, 'NAD'), (0.6, 'NAD & AMPK'), (0.6, 'NAD & DR'), (0.4, 'NAD & FGF21R'), (0.3, 'NAD & Resveratrol')], [(0.5, 1), (1, 'PARP1')]),  # P1A-014, -023, P1B-033, P1C-005, P1A-030, P1A-008
 'SIRT3':     ([(1, 'NAD')], [(0.5, 1)]),                                          # P1A-015
 'SIRT6':     ([(1, 'NAD')], [(0.5, 1)]),                                          # P1A-016
 'SOD2':      ([(0.2, 1), (1, 'SIRT3')], [(0.5, 1)]),                              # P1A-024
 'ROS':       ([(0.6, 1)], [(0.2, 1), (1, 'SOD2')]),          # P1A-025
 'p53':       ([(0.2, 1), (2, 'DNA_damage')], [(0.5, 1), (1, 'SIRT1')]),           # P1A-019
 'NFkB':      ([(0.2, 1), (2, 'Inflammaging')], [(0.5, 1), (1, 'SIRT1'), (1, 'SIRT6')]),  # P1A-020, P1A-027
 'GLP1R':     ([(0.1, 1), (2, 'Semaglutide')], [(1, 1)]),                          # P1C-011
 'Body_weight': ([(1, 'Nutrients')], [(0.3, 1), (2, 'GLP1R'), (1, 'DR')]),        # P1C-012
 'GDF15':     ([(0.2, 1), (2, 'Metformin')], [(1, 1)]),                            # P1B-027
 'Insulin_sensitivity': ([(0.3, 1), (1, 'DR'), (1, 'FGF21R'), (1, 'Metformin'), (1, 'NMN_supp')], [(0.3, 1), (1, 'S6K1'), (1, '!mTORC2')]),  # P1B-035, P1C-010, P1B-028, P1A-032, P1B-016, P1B-022
 'Infection_susceptibility': ([(0.1, 1), (1, 'DR')], [(1, 1)]),                    # P1C-018
 'Impaired_wound_healing': ([(0.1, 1), (0.3, 'DR'), (2, 'Rapamycin')], [(1, 1)]),  # P1C-020 (review-only, low weight), P1C-021
 'Cold_sensitivity': ([(0.1, 1), (1, 'DR'), (0.5, 'FGF21R')], [(1, 1)]),            # P1C-019, P1C-009
 'Vaccine_response': ([(0.5, 1), (1, 'Rapamycin')], [(1, 1)]),                      # P1C-022 (low-dose everolimus)
}

def term(w, e):
    if e == 1: return f'{w}'
    if ' & ' in e: return f'(({e}) ? {w} : 0.0)'
    if e.startswith('!'): return f'({e[1:]} ? 0.0 : {w})'
    return f'({e} ? {w} : 0.0)'

def write(prefix='P1_core'):
    lines = []
    for n in INPUTS:
        lines.append(f'Node {n} {{ logic = {n}; rate_up = 0.0; rate_down = 0.0; }}')
    for n, (up, dn) in N.items():
        lines.append(f'Node {n} {{ rate_up = {" + ".join(term(w, e) for w, e in up)}; rate_down = {" + ".join(term(w, e) for w, e in dn)}; }}')
    open(prefix + '.bnd', 'w').write('\n'.join(lines) + '\n')
    cfg = ['time_tick = 0.5;', 'max_time = 60;', 'sample_count = 4000;', 'discrete_time = 0;', 'use_physrandgen = 0;', 'seed_pseudorandom = 7;', 'thread_count = 4;', 'statdist_traj_count = 0;']
    for n in INPUTS: cfg.append(f'{n}.istate = 0;')
    for n in N: cfg.append(f'{n}.istate = 0;')
    for n in list(N): cfg.append(f'{n}.is_internal = 0;')
    open(prefix + '.cfg', 'w').write('\n'.join(cfg) + '\n')

if __name__ == '__main__':
    write(); print(len(INPUTS) + len(N), 'nodes')
