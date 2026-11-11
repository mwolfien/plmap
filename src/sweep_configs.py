#!/usr/bin/env python3
"""Grid-sweep feature subsets x index aggregations against a features.csv, in one run.

Answers "which setup extracts the most age / age-acceleration signal, and how much is even there?"
For each feature-subset config it runs the cohort through the coupled model, then scores several
index aggregations against (a) chronological age, (b) age-acceleration (clock residual), and the
(c) Aging_inflammation vs CRP EpiScore check (independent only for features.csv files built with
the default N:L-only inflammation input; see the inflammation_def column). It also prints the OLS ceiling -- the
maximum age signal linearly extractable from the input features -- as the reference to beat.

Usage:
  python3 sweep_configs.py --features features.csv                 # full cohort (slow)
  python3 sweep_configs.py --features features.csv --subsample 160 # fast age-stratified screen
  python3 sweep_configs.py --features features.csv --sc 300 --workers 2 --out sweep.json

Findings on GSE40279 (see validation/geo_validation_pipeline.md): the full feature set + the model's
built-in Aging_index node gives r~0 vs age because wrong-signed/flat features and an unstable
aggregation cancel the signal; the "core3" subset (inflammation, adiponectin, hepatic) with a mean of
sub-scores recovers r~0.24 (in-sample OLS ceiling ~0.40). For features.csv files built before inflammation
became N:L-only, infl_vs_CRP is NOT independent: CRP was part of the input (transduction, not validity).
"""
import csv, json, argparse
import numpy as np
import methylation_layer as ML
from multiprocessing import Pool

SUB = ML.SUB
COVERED = ['Aging_NAD', 'Aging_inflammation', 'Aging_secretory', 'Aging_nutrient_sensing', 'Aging_mito']
ALL = ['inflammation', 'immune_aging', 'nad_depletion', 'hepatic', 'adiposity_lep', 'adiponectin', 'damage_smoke']
# feature-subset configs (edit / extend freely)
CONFIGS = {
    'base_all': ALL,
    'drop_adipo': [f for f in ALL if f != 'adiposity_lep'],
    'drop_flat': [f for f in ALL if f not in ('immune_aging', 'damage_smoke')],
    'core3': ['inflammation', 'adiponectin', 'hepatic'],
    'core2': ['inflammation', 'adiponectin'],
}


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float('nan')


def pear(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float); m = np.isfinite(x) & np.isfinite(y)
    return float(np.corrcoef(x[m], y[m])[0, 1]) if m.sum() > 2 else float('nan')


def resid(y, x):
    x = np.asarray(x, float); y = np.asarray(y, float); m = np.isfinite(x) & np.isfinite(y)
    o = np.full_like(y, np.nan, float)
    if m.sum() > 2:
        b, a = np.polyfit(x[m], y[m], 1); o[m] = y[m] - (a + b * x[m])
    return o


def ranknorm(M):
    R = np.empty_like(M)
    for j in range(M.shape[1]):
        R[:, j] = np.argsort(np.argsort(M[:, j])) / (len(M) - 1)
    return R


def _worker(args):
    feat, sc = args                      # sc passed explicitly: globals don't reach spawn-based workers (macOS/Windows)
    return ML.run_patient(feat, sc=sc)


def ols_R(rows, names, age):
    X = np.array([[num(r[c]) for c in names] for r in rows]); m = np.all(np.isfinite(X), 1) & np.isfinite(age)
    if m.sum() < 5:
        return float('nan')
    Xm = np.column_stack([np.ones(m.sum()), X[m]]); beta, *_ = np.linalg.lstsq(Xm, age[m], rcond=None)
    return round(float(pear(Xm @ beta, age[m])), 3)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--features', required=True)
    ap.add_argument('--sc', type=int, default=250, help='MaBoSS sample_count per run (lower = faster screen)')
    ap.add_argument('--subsample', type=int, default=0, help='age-stratified subsample size (0 = full cohort)')
    ap.add_argument('--workers', type=int, default=2)
    ap.add_argument('--out', default='sweep_results.json')
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.features)))
    age = np.array([num(r['age']) for r in rows]); clock = np.array([num(r.get('clock_age', '')) for r in rows])
    vinf = np.array([num(r.get('val_inflammation_episcore', '')) for r in rows])
    if a.subsample and a.subsample < len(rows):
        idx = np.argsort(age)[::max(1, len(rows) // a.subsample)]
        rows = [rows[i] for i in idx]; age = age[idx]; clock = clock[idx]; vinf = vinf[idx]
    n = len(rows); accel = resid(clock, age)
    ceiling = ols_R(rows, [f for f in ALL if np.isfinite([num(r[f]) for r in rows]).any()], age)
    print(f"n={n}  OLS ceiling (all live features -> age) R = {ceiling}\n")
    print("NOTE: configs are chosen on THIS cohort, so the best r_age is optimistic (selection bias);"
          " confirm the pick on held-out samples or a second dataset.\n")
    print(f"{'config':11} {'aggregation':16} {'r_age':>7} {'r_accel':>8} {'infl_vs_CRP':>11}")

    result = {'n': n, 'ols_ceiling_R': ceiling, 'rows': []}
    ci = [SUB.index(s) for s in COVERED]
    for cname, flist in CONFIGS.items():
        feats = [{f: num(r[f]) for f in flist if np.isfinite(num(r[f]))} for r in rows]
        with Pool(a.workers) as p:
            outs = p.map(_worker, [(f, a.sc) for f in feats])
        M = np.array([[o[s] for s in SUB] for o in outs])
        infl = np.array([o['Aging_inflammation'] for o in outs]); ic = pear(vinf, infl)
        aggs = {'model_Aging_index': np.array([o['Aging_index'] for o in outs]),
                'mean_all7': M.mean(1), 'mean_covered': M[:, ci].mean(1),
                'ranknorm_all7': ranknorm(M).mean(1), 'ranknorm_covered': ranknorm(M[:, ci]).mean(1)}
        for an, v in aggs.items():
            ra, rac = pear(age, v), pear(accel, resid(v, age))
            print(f"{cname:11} {an:16} {ra:+7.3f} {rac:+8.3f} {ic:+11.3f}")
            result['rows'].append({'config': cname, 'aggregation': an, 'r_age': round(ra, 3),
                                   'r_accel': round(rac, 3), 'infl_vs_CRP': round(ic, 3)})
        print()
    result['rows'].sort(key=lambda d: -(d['r_age'] if d['r_age'] == d['r_age'] else -9))
    json.dump(result, open(a.out, 'w'), indent=1)
    print('TOP by r_age:')
    for d in result['rows'][:5]:
        print(f"  {d['r_age']:+.3f}  {d['config']} / {d['aggregation']}")
    print(f"\nsaved {a.out}")


if __name__ == '__main__':
    main()
