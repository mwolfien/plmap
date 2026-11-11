#!/usr/bin/env python3
"""Preflight check for a features.csv before the model run.

Catches the usual failure (empty EpiScore columns from a missing/mismatched weight table) at the
*input* stage, so you don't spend the model time to discover a flat signal afterwards. For every model
feature it reports coverage, variance and correlation with age; for every model input it reports what
drives it and its saturation; then it gives a VERDICT and a non-zero exit code if anything is wrong, so
it can gate a pipeline.

Usage:
  python3 check_features.py features.csv            # human report + exit code
  python3 check_features.py features.csv --min-features 5   # require >=5 live features to pass
  python3 check_features.py features.csv --quiet    # verdict + exit code only

Exit codes: 0 = pass; 1 = fail (an empty model feature, too few live features, or a missing required
column); 2 = usage / file error.
"""
import sys, csv, argparse

# feature -> {input: weight}; kept in sync with methylation_layer (imported when available).
_FALLBACK_MAP = {
    'inflammation':  {'Inflammaging': 1.0},
    'immune_aging':  {'Inflammaging': 0.5},
    'hepatic':       {'hyperammonemia': 0.9, 'FGL1': 0.8},
    'adiposity_lep': {'LEP': 1.0, 'Nutrients': 0.6},
    'adiponectin':   {'ADIPOQ': 1.0},
    'damage_smoke':  {'DNA_damage': 1.0},
    'nad_depletion': {'Inflammaging': 0.3},
}
try:                                                 # use the real map/blend if the package imports
    from methylation_layer import FEATURE_MAP as _MAP
except Exception:
    _MAP = _FALLBACK_MAP

MODEL_FEATURES = list(_FALLBACK_MAP)
# nad_depletion has no standard blood EpiScore (see episcore_weights/README.md): empty is tolerated
# (WARN, imputed downstream), not a hard failure. Everything else empty is a failure.
OPTIONAL = {'nad_depletion'}
REQUIRED = [f for f in MODEL_FEATURES if f not in OPTIONAL]
VAL_COLS = ['val_inflammation_episcore', 'val_gdf15_episcore']


def _num(x):
    return float(x) if str(x).strip() not in ('', 'NA', 'NaN', 'nan', 'None') else None


def _finite(xs):
    return [v for v in xs if v is not None and v == v]


def _pearson(xs, ys):
    pairs = [(a, b) for a, b in zip(xs, ys) if a is not None and b is not None and a == a and b == b]
    n = len(pairs)
    if n < 3:
        return None
    sx = sum(a for a, _ in pairs); sy = sum(b for _, b in pairs)
    mx, my = sx / n, sy / n
    cov = sum((a - mx) * (b - my) for a, b in pairs)
    vx = sum((a - mx) ** 2 for a, _ in pairs); vy = sum((b - my) ** 2 for _, b in pairs)
    return cov / (vx * vy) ** 0.5 if vx > 0 and vy > 0 else None


def _rank(vals):                                     # average ranks, ignoring None
    idx = [i for i, v in enumerate(vals) if v is not None and v == v]
    order = sorted(idx, key=lambda i: vals[i])
    r = [None] * len(vals); j = 0
    while j < len(order):
        k = j
        while k + 1 < len(order) and vals[order[k + 1]] == vals[order[j]]:
            k += 1
        avg = (j + k) / 2.0 + 1
        for t in range(j, k + 1):
            r[order[t]] = avg
        j = k + 1
    return r


def _spearman(xs, ys):
    return _pearson(_rank(xs), _rank(ys))


def _std(vals):
    f = _finite(vals); n = len(f)
    if n < 2:
        return None
    m = sum(f) / n
    return (sum((v - m) ** 2 for v in f) / (n - 1)) ** 0.5


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('features_csv')
    ap.add_argument('--min-features', type=int, default=len(REQUIRED),
                    help=f'minimum live model features to pass (default: {len(REQUIRED)} = all required)')
    ap.add_argument('--quiet', action='store_true', help='verdict + exit code only')
    a = ap.parse_args()

    try:
        rows = list(csv.DictReader(open(a.features_csv)))
    except OSError as e:
        print(f'ERROR: cannot read {a.features_csv}: {e}'); return 2
    if not rows:
        print('ERROR: features.csv has no rows'); return 2
    n = len(rows)
    cols = rows[0].keys()

    problems = []                                    # hard failures
    if 'age' not in cols:
        problems.append('required column "age" is missing')
    age = [_num(r.get('age')) for r in rows] if 'age' in cols else [None] * n

    # ---- model features ----
    live, empty, weak = [], [], []
    feat_rows = []
    for f in MODEL_FEATURES:
        if f not in cols:
            empty.append(f); feat_rows.append((f, 0, None, None, None, 'MISSING COLUMN')); continue
        vals = [_num(r.get(f)) for r in rows]
        nf = len(_finite(vals)); cover = nf / n
        if nf == 0:
            empty.append(f); feat_rows.append((f, 0, None, None, None, 'EMPTY')); continue
        live.append(f)
        sd = _std(vals); r = _pearson(vals, age); rho = _spearman(vals, age)
        verdict = 'ok'
        if cover < 0.5:
            verdict = f'WARN low coverage ({cover:.0%})'
        elif sd is not None and sd < 0.02:
            verdict = 'WARN near-constant'
        elif r is not None and abs(r) < 0.05:
            verdict = 'WARN no age signal'; weak.append(f)
        feat_rows.append((f, nf, sd, r, rho, verdict))

    # ---- model inputs (what the network actually receives), via the real blend ----
    input_cover = {}
    for f in live:
        vals = [_num(r.get(f)) for r in rows]
        present = sum(1 for v in vals if v is not None and v == v) / n
        for inp, w in _MAP.get(f, {}).items():
            input_cover.setdefault(inp, []).append((f, w, present))
    all_inputs = sorted({i for m in _MAP.values() for i in m})
    off_inputs = [i for i in all_inputs if i not in input_cover]

    # ---- report ----
    if not a.quiet:
        print(f'Feature preflight — {a.features_csv}  ({n} rows)\n')
        print(f'{"model feature":16}{"n":>7}{"cover":>8}{"std":>8}{"r_age":>8}{"rho_age":>9}  verdict')
        for f, nf, sd, r, rho, v in feat_rows:
            def fmt(x, p='{:.3f}'):
                return p.format(x) if isinstance(x, float) else ('-' if x is None else str(x))
            print(f'{f:16}{nf:>7}{(nf/n if nf else 0):>8.0%}{fmt(sd):>8}{fmt(r):>8}{fmt(rho):>9}  {v}')
        print(f'\nmodel inputs driven: {len(input_cover)}/{len(all_inputs)}'
              + (f'  |  OFF (no data): {", ".join(off_inputs)}' if off_inputs else ''))
        # validation + clock columns
        vstat = []
        for c in VAL_COLS + ['clock_age']:
            if c not in cols:
                vstat.append(f'{c}: missing')
            else:
                nf = len(_finite([_num(r.get(c)) for r in rows]))
                vstat.append(f'{c}: {"empty" if nf == 0 else f"{nf}/{n}"}')
        print('validation/external: ' + '  |  '.join(vstat))
        print()

    # ---- verdict ----
    empty_required = [f for f in empty if f not in OPTIONAL]
    empty_optional = [f for f in empty if f in OPTIONAL]
    if empty_required:
        problems.append(f'{len(empty_required)} required feature(s) empty: ' + ', '.join(empty_required))
    if len(live) < a.min_features:
        problems.append(f'only {len(live)} live features (< required {a.min_features})')

    warnings = []
    if empty_optional:
        warnings.append(f'optional feature(s) empty (imputed downstream): {", ".join(empty_optional)}')
    if weak:
        warnings.append(f'no age signal in {", ".join(weak)} — validate before trusting')

    if problems:
        print('VERDICT: FAIL')
        for p in problems:
            print(f'  - {p}')
        for w in warnings:
            print(f'  (warn) {w}')
        if empty_required:
            print('  -> empty columns usually mean the EpiScore weight tables are absent or their CpG IDs')
            print('     do not overlap the beta-matrix rownames. See episcore_weights/README.md, then')
            print('     re-run preprocess_geo.R and check the per-score "N/M CpGs" overlap lines.')
        return 1
    print('VERDICT: PASS')
    for w in warnings:
        print(f'  (warn) {w}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
