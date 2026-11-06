"""MUS v2 workbook: primary-evidence re-curation + coupled-model results."""
import json
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

V = json.load(open('MUS_module_v3.json'))
SC = json.load(open('coupled_scenarios.json'))
RB = json.load(open('coupled_robustness.json'))
refs = V['references']

HDR = PatternFill('solid', fgColor='1F3864')
SUB = PatternFill('solid', fgColor='D6DCE5')
WHITE = Font(color='FFFFFF', bold=True, size=11)
BOLD = Font(bold=True)
WRAP = Alignment(wrap_text=True, vertical='top')
TOP = Alignment(vertical='top')
thin = Side(style='thin', color='BFBFBF')
BORD = Border(left=thin, right=thin, top=thin, bottom=thin)
MODEL_FILL = {'include': 'E2EFDA', 'use P1 edge/path': 'E2EFDA', 'map only (review-level)': 'FFF2CC',
              'exclude (unsupported)': 'FBE5D6', 'exclude (sign conflict)': 'F8CBAD', 'drop': 'F4B183',
              'membership (merge into parent node)': 'DDEBF7'}

wb = Workbook()


def style_header(ws, row, ncol, fill=HDR, font=WHITE):
    for c in range(1, ncol + 1):
        cell = ws.cell(row=row, column=c)
        cell.fill = fill; cell.font = font; cell.alignment = WRAP; cell.border = BORD


# ---------- 1. README ----------
ws = wb.active; ws.title = 'README'
ws.column_dimensions['A'].width = 120
rows = [
 ('MUS module — Sarcopenia / muscle layer, v3 (all review-only edges resolved)', 14, '1F3864'),
 ('Prevention & Longevity Map · muscle module coupled to the P1 nutrient-sensing/NAD+ core', 11, '808080'),
 ('', 10, None),
 ('What changed in v3', 12, '1F3864'),
 ('All 90 remaining review-only edges are resolved: 47 canonical biochemistry/signalling included as', 10, None),
 ('established model structure, 6 phenotype edges as structural glue, 17 regulatory edges re-curated with', 10, None),
 ('verified primary PubMed papers, 20 kept honestly map-only (no retrievable primary), 1 adipokine sign fixed', 10, None),
 ('(ADIPOQ improves insulin sensitivity). The executable model now rests on 149 primary/canonical edges.', 10, None),
 ('', 10, None),
 ('What changed from v1', 12, '1F3864'),
 ('v1 flagged that 158/191 muscle edges cited only reviews. In v2, 92 edges were re-curated against PRIMARY', 10, None),
 ('literature by three curation agents, 9 edges inherit verified P1 evidence, and all new PMIDs were re-fetched', 10, None),
 ('from PubMed and cross-checked by title. Primary-backed edges rose from 22 to 92 of 191.', 10, None),
 ('', 10, None),
 ('The two broken source PMIDs are resolved:', 11, '1F3864'),
 ('  • 31747227 (Glucose→Glycogen) — not in PubMed. Replaced with Shulman 1990 (PMID 2403659), Tsao 1996', 10, None),
 ('    (8522056) and Bouskila 2007 (18003720); cirrhosis context Selberg 1993 (8486761).', 10, None),
 ('  • 32119304 (Glycogenolysis⊣Glycogen) — not in PubMed. Replaced with Nogales-Gadea 2012 (22730558),', 10, None),
 ('    Real-Martinez 2019 (30914683) for muscle; Rothman 1991 (1948033), Petersen 1999/1998 (10070020/9502760)', 10, None),
 ('    for the liver submap. A full typo screen (56 numeric variants) found no relevant record.', 10, None),
 ('', 10, None),
 ('Tabs', 12, '1F3864'),
 ('  Edges — all 191 reactions with verdict, re-curated PMIDs, E/C/R score and model routing', 10, None),
 ('  Bridge edges — the muscle↔P1 links, now primary-backed (the ones that matter for coupling)', 10, None),
 ('  Coupled model — lean-mass trade-off scenarios on the MUS+P1 MaBoSS model (78 nodes)', 10, None),
 ('  Robustness — sign-stability of the conclusions under ±30% weight noise', 10, None),
 ('  References — every PMID cited in v2, with title/journal/year/DOI', 10, None),
 ('', 10, None),
 ('Evidence is from PubMed; DOIs are in the References tab. Confidence = min(5, max(1, E+C+R)):', 10, '808080'),
 ('E 3 direct / 2 perturbation / 1 correlative-review · C +1 aging context · R +1 ≥2 independent primary PMIDs.', 10, '808080'),
]
for i, (t, sz, col) in enumerate(rows, 1):
    c = ws.cell(row=i, column=1, value=t)
    c.font = Font(bold=(sz >= 12), size=sz, color=col or '000000'); c.alignment = WRAP


# ---------- 2. Edges ----------
ws = wb.create_sheet('Edges')
cols = [('rid', 9), ('source', 18), ('effect', 9), ('target', 20), ('status_vs_P1', 16), ('verdict', 16),
        ('model', 20), ('model_edge (if re-routed)', 24), ('E', 4), ('C', 4), ('R', 4), ('conf', 5),
        ('evidence_type', 13), ('tissue', 22), ('PMIDs (v2)', 22), ('notes', 60)]
ws.append([c[0] for c in cols])
for i, (_, w) in enumerate(cols, 1):
    ws.column_dimensions[get_column_letter(i)].width = w
style_header(ws, 1, len(cols))
for e in sorted(V['reactions'], key=lambda x: x['rid']):
    ws.append([e['rid'], e['source'], e['effect'], e['target'], e.get('status_vs_P1', ''), e['verdict'],
               e['model'], e.get('model_edge', ''), e.get('E'), e.get('C'), e.get('R'), e.get('confidence'),
               e.get('evidence_type', ''), e.get('tissue', ''), ', '.join(e['pmids']), e.get('notes', '')])
    r = ws.max_row
    fill = MODEL_FILL.get(e['model'])
    for c in range(1, len(cols) + 1):
        cell = ws.cell(row=r, column=c); cell.alignment = WRAP if c in (14, 15, 16) else TOP; cell.border = BORD
        if fill and c == 7:
            cell.fill = PatternFill('solid', fgColor=fill)
ws.freeze_panes = 'A2'; ws.auto_filter.ref = f'A1:{get_column_letter(len(cols))}{ws.max_row}'


# ---------- 3. Bridge edges ----------
ws = wb.create_sheet('Bridge edges')
ws.append(['The muscle↔P1 links — these are the edges the coupling depends on, now backed by primary papers'])
ws.cell(1, 1).font = Font(bold=True, size=12, color='1F3864')
ws.append([])
head = ['rid', 'edge', 'verdict', 'conf', 'PMID (key)', 'first author / year', 'what it shows (primary)']
ws.append(head); style_header(ws, 3, len(head))
for w, col in zip([9, 34, 16, 5, 13, 22, 70], 'ABCDEFG'):
    ws.column_dimensions[col].width = w
BRIDGE = ['MUS-023', 'MUS-081', 'MUS-020', 'MUS-022', 'MUS-024', 'MUS-065', 'MUS-096', 'MUS-008', 'MUS-012',
          'MUS-030', 'MUS-033', 'MUS-034', 'MUS-043', 'MUS-047', 'MUS-139', 'MUS-113', 'MUS-140', 'MUS-089',
          'MUS-100', 'MUS-063', 'MUS-071', 'MUS-161', 'MUS-179']
byid = {e['rid']: e for e in V['reactions']}
for rid in BRIDGE:
    e = byid[rid]
    p = e['pmids'][0] if e['pmids'] else ''
    m = refs.get(p, {})
    fa = f"{m.get('fa','')} {m.get('year','')}"
    ws.append([rid, f"{e['source']} {e['effect']} {e['target']}", e['verdict'], e.get('confidence'), p, fa,
               e.get('notes', '')[:300]])
    rr = ws.max_row
    for c in range(1, len(head) + 1):
        ws.cell(rr, c).alignment = WRAP if c == 7 else TOP; ws.cell(rr, c).border = BORD
ws.freeze_panes = 'A4'


# ---------- 4. Coupled model ----------
ws = wb.create_sheet('Coupled model')
ws.append(['Coupled MUS + P1 MaBoSS model — lean-mass trade-off scenarios'])
ws.cell(1, 1).font = Font(bold=True, size=12, color='1F3864')
ws.append(['78 nodes (P1 nutrient/NAD+ core reused + muscle layer), cmaboss_128n. Values = last-state activity 0–1.'])
ws.cell(2, 1).font = Font(italic=True, color='808080')
ws.append(['Read relatively, not absolutely: the model is qualitative. Focus on how a readout MOVES between rows.'])
ws.cell(3, 1).font = Font(italic=True, color='808080')
ws.append([])
rolabels = ['Sarcopenia', 'Catabolism', 'Anabolism', 'Proteolysis', 'Protein_synthesis', 'MSTN', 'mTOR',
            'AMPK', 'FOXO3', 'Autophagy', 'PGC1A', 'Insulin_sensitivity', 'NAD', 'FGF21']
ws.append(['scenario'] + rolabels)
style_header(ws, 5, len(rolabels) + 1)
ws.column_dimensions['A'].width = 28
for i in range(len(rolabels)):
    ws.column_dimensions[get_column_letter(2 + i)].width = 11
base = SC['baseline']
for name, r in SC.items():
    row = [name] + [r.get(c, '') for c in rolabels]
    ws.append(row)
    rr = ws.max_row
    for ci, c in enumerate(rolabels):
        cell = ws.cell(rr, 2 + ci); cell.alignment = TOP; cell.border = BORD
        v = r.get(c); b = base.get(c)
        if isinstance(v, (int, float)) and isinstance(b, (int, float)) and name != 'baseline':
            d = v - b
            if abs(d) >= 0.03:
                cell.fill = PatternFill('solid', fgColor='C6EFCE' if (
                    (d < 0) == (c in ('Sarcopenia', 'Catabolism', 'Proteolysis', 'MSTN'))) else 'FFC7CE')
    ws.cell(rr, 1).border = BORD; ws.cell(rr, 1).font = BOLD if 'Aged' in name or name == 'baseline' else Font()
ws.append([])
ws.append(['Green = change in the lean-mass-favourable direction; red = unfavourable (vs baseline).'])
ws.cell(ws.max_row, 1).font = Font(italic=True, color='808080')
ws.freeze_panes = 'B6'


# ---------- 5. Robustness ----------
ws = wb.create_sheet('Robustness')
ws.append([f"Sign-stability under ±{int(RB['eps']*100)}% random weight noise, {RB['reps']} reps"])
ws.cell(1, 1).font = Font(bold=True, size=12, color='1F3864')
ws.append([])
head = ['conclusion', 'baseline Δ', 'expected sign', '% reps preserved', 'mean Δ']
ws.append(head); style_header(ws, 3, len(head))
for w, col in zip([52, 12, 13, 16, 10], 'ABCDE'):
    ws.column_dimensions[col].width = w
for name, d in RB['contrasts'].items():
    ws.append([name, d['baseline_diff'], '+' if d['expected_sign'] > 0 else '−',
               f"{d['frac_preserved']*100:.0f}%", d['mean_diff']])
    rr = ws.max_row
    for c in range(1, 6):
        ws.cell(rr, c).alignment = WRAP if c == 1 else TOP; ws.cell(rr, c).border = BORD
    if d['frac_preserved'] >= 0.9:
        ws.cell(rr, 4).fill = PatternFill('solid', fgColor='C6EFCE')


# ---------- 6. References ----------
ws = wb.create_sheet('References')
head = ['PMID', 'first author', 'year', 'journal', 'title', 'DOI']
ws.append(head); style_header(ws, 1, len(head))
for w, col in zip([12, 18, 7, 24, 70, 32], 'ABCDEF'):
    ws.column_dimensions[col].width = w
for p in sorted(refs, key=lambda x: int(x) if x.isdigit() else 0):
    m = refs[p]
    ws.append([p, m.get('fa', ''), m.get('year', ''), m.get('journal', ''), m.get('title', ''),
               m.get('doi', '')])
    rr = ws.max_row
    for c in range(1, 7):
        ws.cell(rr, c).alignment = WRAP if c == 5 else TOP; ws.cell(rr, c).border = BORD
ws.freeze_panes = 'A2'; ws.auto_filter.ref = f'A1:F{ws.max_row}'

out = 'MUS_module_v3_resolved.xlsx'
wb.save(out)
print('saved', out, '| sheets', wb.sheetnames)
