#!/usr/bin/env python3
"""Turn a published long-format EpiScore weights table into the per-score CSVs that
preprocess_geo.R expects (one file per score, columns: cpg,weight).

Works with any "long" weights table that has a CpG column, a weight/coefficient column, and
(optionally) a predictor/protein-name column — e.g.:
  - MethylDetectR  Predictors_Shiny_by_Groups.csv  (CpG_Site,Coefficient,Mean_Beta_Value,Predictor)
  - Gadd 2022 eLife Supplementary file 1E / the DanniGadd/EpiScores-for-protein-levels weights
  - a single-predictor table (just cpg,weight) for e.g. a smoking score
Columns are auto-detected; override with --cpg-col/--weight-col/--predictor-col if needed.

Typical workflow
  # 1. see which predictors the file contains (exact spellings) and their CpG counts:
  python3 make_episcore_weights.py Predictors_Shiny_by_Groups.csv --list

  # 2. (recommended) check overlap with YOUR array first, to avoid "0 CpGs overlap -> NA".
  #    In R once:  writeLines(rownames(beta), "beta_cpgs.txt")
  python3 make_episcore_weights.py Predictors_Shiny_by_Groups.csv --list --beta-cpgs beta_cpgs.txt

  # 3. write the per-score files (left = output file the pipeline wants, right = predictor name in THIS file):
  python3 make_episcore_weights.py Predictors_Shiny_by_Groups.csv --outdir episcore_weights \
      --map CRP.csv="CRP" GDF15.csv="GDF15" ADIPOQ.csv="Adiponectin" LEP.csv="Leptin" liver.csv="Liver ALT"

  # a single-predictor file (no predictor column) -> one output, no --map key needed:
  python3 make_episcore_weights.py smoking_weights.csv --outdir episcore_weights --single smoking.csv

Then: Rscript preprocess_geo.R --weights episcore_weights   (run it where that folder is)
"""
import sys, csv, re, os, argparse

CG = re.compile(r'^cg\d{6,}$', re.I)                       # Illumina 450K/EPIC probe id
INTERCEPT = {'(intercept)', 'intercept', 'const', '(const)'}


def detect_cols(header, rows):
    """Return (cpg_col, weight_col, predictor_col or None) by header name, then value sniffing."""
    low = {h: h.lower() for h in header}
    def find(pats):
        for p in pats:
            for h in header:
                if re.search(p, low[h]):
                    return h
        return None
    # weight: prefer explicit coefficient/weight/estimate over anything containing "beta"
    weight = find([r'^coef', r'coefficient', r'^weight$', r'weight', r'estimate', r'\bbeta\b(?!.*mean)'])
    cpg = find([r'^cpg$', r'cpg[_ ]?site', r'^probe', r'probe[_ ]?id', r'^ivl$', r'^cg$', r'marker'])
    pred = find([r'predictor', r'protein', r'^score$', r'trait', r'^name$', r'biomarker', r'episcore'])
    # value-sniff the CpG column if not named
    if cpg is None and rows:
        for h in header:
            vals = [r[h] for r in rows[:50] if r.get(h)]
            if vals and sum(bool(CG.match(str(v))) for v in vals) >= max(3, len(vals) // 2):
                cpg = h; break
    # weight fallback: the first numeric column that isn't the cpg column
    if weight is None and rows:
        for h in header:
            if h == cpg:
                continue
            ok = 0
            for r in rows[:50]:
                try:
                    float(r[h]); ok += 1
                except (TypeError, ValueError):
                    pass
            if ok >= 25:
                weight = h; break
    return cpg, weight, pred


def load(path):
    for enc in ('utf-8-sig', 'latin-1'):                 # some cohort files are latin-1
        try:
            with open(path, newline='', encoding=enc) as f:
                rows = list(csv.DictReader(f))
            break
        except UnicodeDecodeError:
            continue
    else:
        sys.exit(f'ERROR: could not decode {path} as utf-8 or latin-1')
    if not rows:
        sys.exit(f'ERROR: {path} has no data rows')
    return list(rows[0].keys()), rows


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('weights_csv')
    ap.add_argument('--outdir', default='episcore_weights')
    ap.add_argument('--map', nargs='*', default=[], metavar='FILE=PREDICTOR',
                    help='output file = predictor name in the weights file (repeatable)')
    ap.add_argument('--single', metavar='FILE.csv',
                    help='weights file has no predictor column: write all rows to this one output')
    ap.add_argument('--list', action='store_true', help='list predictors + CpG counts and exit')
    ap.add_argument('--beta-cpgs', help='file with your array CpG IDs (one per line) to report overlap')
    ap.add_argument('--cpg-col'); ap.add_argument('--weight-col'); ap.add_argument('--predictor-col')
    a = ap.parse_args()

    header, rows = load(a.weights_csv)
    cpg = a.cpg_col or None; wgt = a.weight_col or None; pred = a.predictor_col or None
    d_cpg, d_wgt, d_pred = detect_cols(header, rows)
    cpg = cpg or d_cpg; wgt = wgt or d_wgt; pred = pred or d_pred
    if not cpg or not wgt:
        sys.exit(f'ERROR: could not detect columns. Header = {header}\n'
                 f'  pass --cpg-col / --weight-col (and --predictor-col) explicitly.')
    print(f'columns: cpg="{cpg}"  weight="{wgt}"  predictor="{pred or "(none)"}"')

    beta = None
    if a.beta_cpgs:
        with open(a.beta_cpgs) as f:
            beta = {ln.strip() for ln in f if ln.strip()}
        print(f'array CpGs loaded: {len(beta)}')

    def overlap(cpgs):
        return '' if beta is None else f'  overlap={sum(c in beta for c in cpgs)}/{len(cpgs)}'

    # group rows by predictor (or a single bucket)
    groups = {}
    for r in rows:
        key = (r.get(pred) or '').strip() if pred else '__all__'
        c = (r.get(cpg) or '').strip()
        if not c or c.lower() in INTERCEPT:
            groups.setdefault(key, {'intercept': None, 'rows': []})
            if c.lower() in INTERCEPT:
                groups[key]['intercept'] = num(r.get(wgt))
            continue
        g = groups.setdefault(key, {'intercept': None, 'rows': []})
        w = num(r.get(wgt))
        if w is not None:
            g['rows'].append((c, w))

    if a.list or (not a.map and not a.single):
        print('\npredictors in this file (name : n CpGs):')
        for k in sorted(groups, key=lambda k: -len(groups[k]['rows'])):
            print(f'  {k!r:40} {len(groups[k]["rows"]):5d}{overlap([c for c, _ in groups[k]["rows"]])}')
        if not a.list:
            print('\n(nothing written — pass --map FILE=PREDICTOR ... or --single FILE.csv)')
        return 0

    os.makedirs(a.outdir, exist_ok=True)
    # case-insensitive predictor lookup
    by_lower = {k.lower(): k for k in groups}
    jobs = []
    if a.single:
        jobs.append((a.single, '__all__'))
    for item in a.map:
        if '=' not in item:
            sys.exit(f'ERROR: --map entry "{item}" must be FILE=PREDICTOR')
        fname, predname = item.split('=', 1)
        key = predname if predname in groups else by_lower.get(predname.lower())
        if key is None:
            print(f'  SKIP {fname}: predictor {predname!r} not found (use --list to see names)')
            continue
        jobs.append((fname, key))

    for fname, key in jobs:
        g = groups.get(key, {'intercept': None, 'rows': []})
        out = os.path.join(a.outdir, fname)
        with open(out, 'w', newline='') as f:
            w = csv.writer(f); w.writerow(['cpg', 'weight'])
            if g['intercept'] is not None:
                w.writerow(['(Intercept)', g['intercept']])
            for c, val in g['rows']:
                w.writerow([c, val])
        print(f'  wrote {out}: {len(g["rows"])} CpGs'
              + ('' if g['intercept'] is None else ' (+intercept)')
              + overlap([c for c, _ in g['rows']]))
    print('\nNext: Rscript preprocess_geo.R --weights ' + a.outdir)
    return 0


if __name__ == '__main__':
    sys.exit(main())
