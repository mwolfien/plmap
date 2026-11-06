"""Characterization figures: (1) input x readout effect heatmap, (2) readout dynamic-range bars
showing where the model has signal vs is compressed (boundary)."""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

O = json.load(open('coupled_characterization.json'))
A = O['single_input_delta']
RNG = O['readout_dynamic_range']
INF = O['input_total_influence']

readouts = list(RNG.keys())                      # sorted by dynamic range desc
inputs = list(INF.keys())                        # sorted by influence desc
M = np.array([[A[i][r] for r in readouts] for i in inputs])

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9})
fig = plt.figure(figsize=(13, 9))
gs = fig.add_gridspec(1, 2, width_ratios=[3.1, 1], wspace=0.32)

# ---- heatmap ----
ax = fig.add_subplot(gs[0])
vmax = 0.6
im = ax.imshow(M, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
ax.set_xticks(range(len(readouts))); ax.set_xticklabels(readouts, rotation=55, ha='right', fontsize=8)
ax.set_yticks(range(len(inputs))); ax.set_yticklabels(inputs, fontsize=8)
for yi in range(len(inputs)):
    for xi in range(len(readouts)):
        v = M[yi, xi]
        if abs(v) >= 0.15:
            ax.text(xi, yi, f'{v:+.2f}', ha='center', va='center', fontsize=6,
                    color='white' if abs(v) > 0.38 else '#222')
ax.set_title('What each lever does in isolation\n(Δ activity vs baseline; single input ON)', fontsize=10, pad=8)
ax.set_xlabel('readout'); ax.set_ylabel('input (by total influence)')
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01); cb.set_label('Δ vs baseline', fontsize=8)

# ---- dynamic-range bars ----
ax2 = fig.add_subplot(gs[1])
vals = [RNG[r] for r in readouts]
colors = ['#c0392b' if v < 0.05 else '#e08a2e' if v < 0.12 else '#2e6da4' for v in vals]
y = range(len(readouts))
ax2.barh(list(y), vals, color=colors, height=0.72)
ax2.set_yticks(list(y)); ax2.set_yticklabels(readouts, fontsize=8)
ax2.invert_yaxis()
ax2.axvline(0.05, ls=':', color='#c0392b', lw=1)
ax2.set_title('Dynamic range per readout\n(max−min across all single-input sweeps)', fontsize=10, pad=8)
ax2.set_xlabel('range of activity (0–1)')
for i, v in enumerate(vals):
    ax2.text(v + 0.008, i, f'{v:.2f}', va='center', fontsize=7)
ax2.text(0.052, len(readouts) - 0.5, 'compressed → model\ncannot yet move these', color='#c0392b', fontsize=7, va='top')
for s in ['top', 'right']:
    ax2.spines[s].set_visible(False)

fig.suptitle('Coupled MUS–P1 model — systematic behaviour map (isolated input battery)', fontsize=12, y=0.99)
fig.savefig('fig_behaviour_map.png', dpi=170, bbox_inches='tight')
print('saved fig_behaviour_map.png', M.shape)

# ---- dose-response small multiples ----
B = O['dose_response']
levers = list(B.keys())
metrics = ['Anabolism', 'Sarcopenia', 'MSTN', 'NAD', 'Insulin_sensitivity', 'Autophagy']
fig2, axes = plt.subplots(2, 4, figsize=(14, 6.4), sharex=True)
xs = [0.0, 0.25, 0.5, 0.75, 1.0]
pal = ['#2e6da4', '#c0392b', '#2e8b57', '#8e44ad', '#e08a2e', '#1a9e9e']
for ax, lv in zip(axes.flat, levers):
    for mi, m in enumerate(metrics):
        ys = [B[lv][str(p) if str(p) in B[lv] else p][m] for p in xs]
        ax.plot(xs, ys, marker='o', ms=4, lw=1.6, color=pal[mi], label=m)
    ax.set_title(lv, fontsize=9); ax.set_ylim(0, 1); ax.grid(alpha=0.25, lw=0.5)
    ax.set_xlabel('intervention intensity (fraction ON)', fontsize=7)
for ax in axes.flat[len(levers):]:
    ax.axis('off')
axes.flat[0].legend(fontsize=7, loc='upper left', framealpha=0.9)
fig2.suptitle('Dose-response (fractional activation) of key levers', fontsize=12)
fig2.tight_layout(rect=[0, 0, 1, 0.96])
fig2.savefig('fig_dose_response.png', dpi=150, bbox_inches='tight')
print('saved fig_dose_response.png')
