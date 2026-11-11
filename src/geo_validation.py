"""GEO validation for the coupled model's aging layer.

TARGET DATASET: GSE40279 (Hannum 2013, PMID 23177740) - 656 whole-blood 450K samples, ages 19-101.
Chosen because: whole blood (matches the blood-readout design), wide age span, public, clock-grade.

TWO MODES
  (A) REAL DATA:   python3 geo_validation.py --features features.csv
        Consumes the per-sample feature matrix produced by preprocess_gse40279.R (R/Bioconductor:
        EpiDISH cell composition + protein EpiScores + an epigenetic clock -> cohort-percentile
        features). Runs every sample through the coupled model, then reports:
          - Aging_index + sub-scores vs chronological age (Pearson/Spearman)
          - Aging_index vs the clock (external, non-circular - the clock is NOT a model input)
          - NON-CIRCULAR molecular checks: Aging_inflammation vs CRP/IL-6 EpiScore,
            Aging_secretory vs GDF15 EpiScore
          - age-acceleration: residual(clock~age) vs residual(index~age)
  (B) SIMULATED:   python3 geo_validation.py        (no --features)
        A literature-calibrated synthetic cohort that exercises the full machinery end-to-end and
        shows the EXPECTED signal (aging index rises with age) - the pipeline proof used before the
        real matrix is attached. Clearly a simulation, not real data.

The real-data branch is the "~15-line" wiring promised in SETUP.md, implemented properly: features.csv
rows -> methylation_layer.run_patient -> correlations. Heavy preprocessing (A's R step) runs where the
data + RAM live; this script is light Python and runs on a laptop.
"""
import os, sys, csv, json, argparse
import numpy as np
import methylation_layer as ML

MODEL_FEATURES = ['inflammation', 'immune_aging', 'nad_depletion', 'hepatic',
                  'adiposity_lep', 'adiponectin', 'damage_smoke']
METRICS = ['Aging_index'] + ML.SUB + ['Sarcopenia']


# ----------------------------------------------------------------------------- stats helpers
def pearson(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    return float(np.corrcoef(x[m], y[m])[0, 1]) if m.sum() > 2 else float('nan')


def spearman(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() <= 2:
        return float('nan')
    rx = np.argsort(np.argsort(x[m])); ry = np.argsort(np.argsort(y[m]))
    return float(np.corrcoef(rx, ry)[0, 1])


def residual(y, x):
    """residual of y regressed on x (for age-acceleration)."""
    m = np.isfinite(x) & np.isfinite(y); out = np.full_like(y, np.nan, dtype=float)
    if m.sum() > 2:
        b, a = np.polyfit(x[m], y[m], 1); out[m] = y[m] - (a + b * x[m])
    return out


def run_cohort(cohort, ages):
    """cohort: list of {'readout': {...}}; returns {metric: {pearson,spearman,slope_per_decade}}."""
    ag = np.asarray(ages, float); corr = {}
    for m in METRICS:
        y = np.array([c['readout'][m] for c in cohort], float)
        sl = float(np.polyfit(ag[np.isfinite(y)], y[np.isfinite(y)], 1)[0] * 10)
        corr[m] = {'pearson_r': round(pearson(ag, y), 3), 'spearman_r': round(spearman(ag, y), 3),
                   'slope_per_decade': round(sl, 3)}
    return corr


# ============================================================================= REAL DATA
def real_mode(path, use_features=None):
    rows = list(csv.DictReader(open(path)))
    n = len(rows)
    feature_list = [f for f in MODEL_FEATURES if use_features is None or f in use_features]
    print(f'REAL DATA: {path}  ({n} samples)')
    if use_features is not None:
        print(f'  feature subset: {", ".join(feature_list)}')

    # parse columns; impute each model feature's NA/blank to the cohort median (fallback 0.5)
    def col(name):
        return np.array([float(r[name]) if r.get(name, '') not in ('', 'NA', 'NaN', None) else np.nan
                         for r in rows], float)

    ages = col('age')
    clock = col('clock_age') if 'clock_age' in rows[0] else np.full(n, np.nan)
    feats, excluded, imputed = {}, [], []
    for f in feature_list:
        v = col(f) if f in rows[0] else np.full(n, np.nan)
        if not np.isfinite(v).any():                 # whole column empty -> feature UNUSED
            excluded.append(f)                       # (its inputs stay OFF, never imputed to a constant)
            continue
        nmiss = int(np.isnan(v).sum())
        if nmiss:                                     # partial gaps -> per-sample cohort-median impute
            med = float(np.nanmedian(v))
            v = np.where(np.isnan(v), med, v)
            imputed.append(f'{f} ({nmiss}/{n} -> {med:.2f})')
        feats[f] = np.clip(v, 0.0, 1.0)
    if excluded:
        print('  EXCLUDED (no data in column; input left OFF, NOT imputed): ' + ', '.join(excluded))
    if imputed:
        print('  partially imputed (cohort median): ' + ', '.join(imputed))
    if len(feats) <= max(1, len(feature_list) // 2):
        print(f'  WARNING: only {len(feats)}/{len(feature_list)} features have data '
              f'({", ".join(feats)}). Headline correlations reflect a degraded input set; '
              f'supply the EpiScore weight tables (see episcore_weights/README.md).')

    # run every sample through the coupled model (only features with data are passed;
    # missing ones are simply absent so their inputs are not driven)
    cohort = []
    for i in range(n):
        feat = {f: float(feats[f][i]) for f in feats}
        r = ML.run_patient(feat, sc=500)
        cohort.append({'sample_id': rows[i].get('sample_id', str(i)), 'age': float(ages[i]),
                       'clock_age': float(clock[i]) if np.isfinite(clock[i]) else None,
                       'features': feat, 'readout': r})
        if (i + 1) % 50 == 0:
            print(f'    ...{i+1}/{n}')

    corr_age = run_cohort(cohort, ages)
    print(f"\n{'readout':26} {'r vs age':>10} {'rho vs age':>11} {'/decade':>9}")
    for m in METRICS:
        c = corr_age[m]
        print(f"{m:26} {c['pearson_r']:10.3f} {c['spearman_r']:11.3f} {c['slope_per_decade']:9.3f}")

    # robust aggregate: mean of the 7 hallmark sub-scores (more stable than the model's Aging_index node,
    # which can collapse under some feature sets; see validation/geo_validation_pipeline.md Findings)
    SUBM = np.array([[c['readout'][s] for s in ML.SUB] for c in cohort], float)
    idx_mean = SUBM.mean(1)
    corr_age['Aging_index_mean'] = {'pearson_r': round(pearson(ages, idx_mean), 3),
                                    'spearman_r': round(spearman(ages, idx_mean), 3)}

    # OLS ceiling: maximum age signal linearly extractable from the INPUT features (reference line)
    def _ols_R(names):
        if not names:
            return float('nan')
        X = np.column_stack([feats[f] for f in names] + [np.ones(n)])
        m = np.all(np.isfinite(X), 1) & np.isfinite(ages)
        if m.sum() < 5:
            return float('nan')
        beta, *_ = np.linalg.lstsq(X[m], ages[m], rcond=None)
        return round(float(pearson(X[m] @ beta, ages[m])), 3)
    ols = _ols_R(list(feats))
    print(f"\nAging_index (model node)      r vs age = {corr_age['Aging_index']['pearson_r']:+.3f}")
    print(f"Aging_index_mean (sub-scores) r vs age = {corr_age['Aging_index_mean']['pearson_r']:+.3f}"
          f"   [OLS ceiling from features = {ols}]")

    out = {'dataset': 'GSE40279 (real)', 'n': n, 'corr_vs_age': corr_age, 'ols_ceiling_R': ols,
           'features_used': list(feats), 'excluded_features': excluded,
           'imputed_features': imputed, 'cohort': cohort}

    idx = idx_mean   # use the robust mean index for clock / age-acceleration comparisons

    # clock comparison (external, non-circular)
    if np.isfinite(clock).sum() > 2:
        out['Aging_index_vs_clock'] = {'pearson_r': round(pearson(clock, idx), 3)}
        acc_clock = residual(clock, ages); acc_idx = residual(idx, ages)
        out['age_acceleration'] = {'pearson_r_resid_index_vs_resid_clock':
                                   round(pearson(acc_clock, acc_idx), 3),
                                   'note': 'do people the clock calls old-for-age also score high mechanistically?'}
        print(f"\nAging_index vs clock r = {out['Aging_index_vs_clock']['pearson_r']}"
              f" | age-accel (resid-resid) r = {out['age_acceleration']['pearson_r_resid_index_vs_resid_clock']}")

    # NON-CIRCULAR molecular checks
    checks = {}
    if 'val_inflammation_episcore' in rows[0]:
        v = col('val_inflammation_episcore')
        y = np.array([c['readout']['Aging_inflammation'] for c in cohort], float)
        checks['Aging_inflammation_vs_CRP_IL6_episcore'] = round(pearson(v, y), 3)
    if 'val_gdf15_episcore' in rows[0]:
        v = col('val_gdf15_episcore')
        y = np.array([c['readout']['Aging_secretory'] for c in cohort], float)
        checks['Aging_secretory_vs_GDF15_episcore'] = round(pearson(v, y), 3)
    if checks:
        out['noncircular_checks'] = checks
        print('\nnon-circular molecular checks (sub-score vs independent EpiScore):')
        for k, v in checks.items():
            print(f'  {k:46} r = {v}')

    # a test-subset run (features_test.csv) gets its own outputs so it never overwrites the full run
    tag = '_test' if os.path.basename(path).startswith('features_test') else ''
    fjson, fpng = f'geo_validation_real{tag}.json', f'fig_geo_validation_real{tag}.png'
    json.dump(out, open(fjson, 'w'), indent=1, default=str)
    _figure(cohort, ages, title_suffix='GSE40279 (real data' + (', TEST subset)' if tag else ')'), fname=fpng)
    print(f'\nsaved {fjson} + {fpng}')
    return out


# ============================================================================= SIMULATED
def simulated_mode():
    rng = np.random.default_rng(42)
    N = 60
    ages = rng.uniform(22, 88, N)

    def frac_at_age(age, slope_lo, slope_hi, noise=0.12, lo=0.03, hi=0.97):
        t = (age - 20) / 70.0
        return float(np.clip(slope_lo + (slope_hi - slope_lo) * t + rng.normal(0, noise), lo, hi))

    cohort = []
    for age in ages:
        smoker = rng.random() < 0.25
        feat = {'inflammation':  frac_at_age(age, 0.18, 0.78),
                'immune_aging':  frac_at_age(age, 0.12, 0.85),
                'nad_depletion': frac_at_age(age, 0.20, 0.80),
                'hepatic':       frac_at_age(age, 0.15, 0.45, noise=0.15),
                'adiposity_lep': frac_at_age(age, 0.30, 0.50, noise=0.18),
                'adiponectin':   frac_at_age(age, 0.65, 0.40, noise=0.15),
                'damage_smoke':  float(np.clip((0.75 if smoker else 0.15) + rng.normal(0, 0.1), 0.03, 0.97))}
        cohort.append({'age': float(age), 'smoker': bool(smoker), 'features': feat,
                       'readout': ML.run_patient(feat, sc=500)})

    corr = run_cohort(cohort, ages)
    ag = np.array(ages)
    print(f'SIMULATED literature-calibrated cohort, N={N}, ages {ag.min():.0f}-{ag.max():.0f}')
    print(f"{'readout':26} {'r vs age':>10} {'/decade':>9}")
    for m in METRICS:
        print(f"{m:26} {corr[m]['pearson_r']:10.3f} {corr[m]['slope_per_decade']:9.3f}")
    json.dump({'dataset_target': 'GSE40279 (simulation pending real matrix)', 'n': N,
               'corr_vs_age': corr, 'cohort': cohort}, open('geo_validation.json', 'w'), indent=1)
    _figure(cohort, ages, title_suffix='literature-calibrated simulated cohort (pipeline proof; GSE40279 pending)',
            fname='fig_geo_validation.png')
    print('\nsaved geo_validation.json + fig_geo_validation.png')
    return corr


# ============================================================================= figure
def _figure(cohort, ages, title_suffix, fname):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    ag = np.asarray(ages, float)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, m, col in zip(axes, ['Aging_index', 'Aging_NAD', 'Aging_inflammation'],
                          ['#16233B', '#2E8B7F', '#C0392B']):
        y = np.array([c['readout'][m] for c in cohort], float)
        ax.scatter(ag, y, c=col, alpha=0.6, s=28, edgecolor='white', linewidth=0.5)
        mfin = np.isfinite(ag) & np.isfinite(y)
        b, a = np.polyfit(ag[mfin], y[mfin], 1); xs = np.array([ag.min(), ag.max()])
        ax.plot(xs, a + b * xs, color=col, lw=2)
        ax.set_title(f'{m.replace("Aging_", "Aging: ")}  (r={pearson(ag, y):.2f})', fontsize=11)
        ax.set_xlabel('chronological age (y)'); ax.set_ylabel('model readout (0-1)')
        ax.grid(alpha=0.25, lw=0.5); ax.set_ylim(0, 1)
    fig.suptitle(f'Coupled-model aging readouts vs age - {title_suffix}', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(fname, dpi=160, bbox_inches='tight')


# ============================================================================= main
if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--features', help='features.csv from preprocess_gse40279.R (real-data mode)')
    ap.add_argument('--use-features', help='comma-separated subset of model features to use, e.g. '
                    '"inflammation,adiponectin,hepatic" (the GSE40279-recommended "core3"); default: all present')
    a = ap.parse_args()
    if a.features:
        uf = [s.strip() for s in a.use_features.split(',')] if a.use_features else None
        real_mode(a.features, use_features=uf)
    else:
        simulated_mode()
