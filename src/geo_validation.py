"""GEO validation scaffold for the coupled model's aging layer.

TARGET DATASET: GSE40279 (Hannum 2013, PMID 23177740) — 656 whole-blood 450K samples, ages 19-101.
Chosen because: whole blood (matches the blood-readout design), wide age span, public, clock-grade.

Two parts:
  (1) REAL-DATA PIPELINE (documented in geo_validation_pipeline.md) — the exact steps to derive the
      methylation features the model layer needs (EpiDISH cell composition, protein EpiScores, an
      epigenetic clock) from the beta matrix and correlate the model aging index with age/clock.
      Requires the beta matrix + R methylation stack -> run where the data lives.
  (2) THIS SCRIPT — a literature-calibrated SIMULATED cohort that exercises the full machinery
      end-to-end and shows the EXPECTED validation signal (aging index rises with age), so the
      pipeline is proven before the real matrix is attached. Clearly a simulation, not real data.

Age->feature slopes use documented blood-aging directions (percentile at age, 20->90):
  inflammation / immune_aging : rise with age (inflammaging; naive-T decline)   [Franceschi; Varadhan 2014]
  nad_depletion               : rises (CD38 up, NAD down with age)              [Camacho-Pereira 2016]
  hepatic, adiposity_lep      : mild midlife rise                               [general]
  adiponectin(protective)     : modeled via ADIPOQ input, mild age decline
  damage_smoke                : smoking status, age-independent (random)
"""
import numpy as np, json
import methylation_layer as ML

rng = np.random.default_rng(42)
N = 60
ages = rng.uniform(22, 88, N)


def frac_at_age(age, slope_lo, slope_hi, noise=0.12, lo=0.03, hi=0.97):
    t = (age - 20) / 70.0
    return float(np.clip(slope_lo + (slope_hi - slope_lo) * t + rng.normal(0, noise), lo, hi))


cohort = []
for age in ages:
    smoker = rng.random() < 0.25
    feat = {
        'inflammation':  frac_at_age(age, 0.18, 0.78),
        'immune_aging':  frac_at_age(age, 0.12, 0.85),
        'nad_depletion': frac_at_age(age, 0.20, 0.80),
        'hepatic':       frac_at_age(age, 0.15, 0.45, noise=0.15),
        'adiposity_lep': frac_at_age(age, 0.30, 0.50, noise=0.18),
        'adiponectin':   frac_at_age(age, 0.65, 0.40, noise=0.15),   # declines with age
        'damage_smoke':  float(np.clip((0.75 if smoker else 0.15) + rng.normal(0, 0.1), 0.03, 0.97)),
    }
    r = ML.run_patient(feat, sc=500)
    cohort.append({'age': float(age), 'smoker': bool(smoker), 'features': feat, 'readout': r})

json.dump({'dataset_target': 'GSE40279 (simulation pending real matrix)', 'n': N, 'cohort': cohort},
          open('geo_validation.json', 'w'), indent=1)

ag = np.array([c['age'] for c in cohort])
metrics = ['Aging_index'] + ML.SUB + ['Sarcopenia']
print(f'Literature-calibrated simulated cohort, N={N}, ages {ag.min():.0f}-{ag.max():.0f}')
print(f"{'readout':26} {'Pearson r vs age':>16} {'slope/decade':>13}")
corr = {}
for m in metrics:
    y = np.array([c['readout'][m] for c in cohort])
    r = float(np.corrcoef(ag, y)[0, 1]); sl = float(np.polyfit(ag, y, 1)[0] * 10)
    corr[m] = {'pearson_r': round(r, 3), 'slope_per_decade': round(sl, 3)}
    print(f'{m:26} {r:16.3f} {sl:13.3f}')
json.dump(corr, open('geo_validation_corr.json', 'w'), indent=1)

# figure: aging index vs age + a couple of subscores
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
for ax, m, col in zip(axes, ['Aging_index', 'Aging_NAD', 'Aging_inflammation'],
                      ['#16233B', '#2E8B7F', '#C0392B']):
    y = np.array([c['readout'][m] for c in cohort])
    ax.scatter(ag, y, c=col, alpha=0.6, s=28, edgecolor='white', linewidth=0.5)
    b, a = np.polyfit(ag, y, 1); xs = np.array([ag.min(), ag.max()])
    ax.plot(xs, a + b * xs, color=col, lw=2)
    ax.set_title(f'{m.replace("Aging_","Aging: ")}  (r={np.corrcoef(ag,y)[0,1]:.2f})', fontsize=11)
    ax.set_xlabel('chronological age (y)'); ax.set_ylabel('model readout (0-1)')
    ax.grid(alpha=0.25, lw=0.5); ax.set_ylim(0, 1)
fig.suptitle('Coupled-model aging readouts vs age — literature-calibrated simulated cohort (pipeline proof; GSE40279 pending)', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig('fig_geo_validation.png', dpi=160, bbox_inches='tight')
print('\nsaved fig_geo_validation.png')
