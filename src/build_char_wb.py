"""Model characterization workbook: behaviour map, single-input matrix, dose-response,
interactions, patient archetypes, boundaries & the methylation-personalization extension map."""
import json
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
from openpyxl.drawing.image import Image as XLImage

O = json.load(open('coupled_characterization.json'))
HDR = PatternFill('solid', fgColor='1F3864'); WHITE = Font(color='FFFFFF', bold=True)
WRAP = Alignment(wrap_text=True, vertical='top'); TOP = Alignment(vertical='top')
thin = Side(style='thin', color='CCCCCC'); BORD = Border(left=thin, right=thin, top=thin, bottom=thin)
CEN = Alignment(horizontal='center', vertical='center')


def hrow(ws, row, n):
    for c in range(1, n + 1):
        x = ws.cell(row=row, column=c); x.fill = HDR; x.font = WHITE; x.alignment = WRAP; x.border = BORD


def heat(v, scale=0.6):
    if v is None or v == '':
        return None
    t = max(-1, min(1, v / scale))
    if t >= 0:
        r, g, b = 255, int(255 - t * 120), int(255 - t * 150)
    else:
        t = -t; r, g, b = int(255 - t * 150), int(255 - t * 90), 255
    return PatternFill('solid', fgColor=f'{r:02X}{g:02X}{b:02X}')


wb = Workbook()

# 1. Overview
ws = wb.active; ws.title = 'Overview'; ws.column_dimensions['A'].width = 118
rows = [
 ('Coupled MUS–P1 model — systematic behaviour characterization', 14, '1F3864'),
 ('Isolated-setup battery: what the model can show, and where its current boundaries are', 11, '808080'),
 ('', 10, None),
 ('Why this exists', 12, '1F3864'),
 ('Before using the model to interpret patient methylation profiles, we need to know what signals it can', 10, None),
 ('actually produce and where it is flat. This runs the 22 inputs through isolated and combined setups and', 10, None),
 ('maps every effect on 21 readouts, so reviewers and clinicians see real signals and honest limits.', 10, None),
 ('', 10, None),
 ('Headline findings', 12, '1F3864'),
 ('1. Strong, mechanistically correct signals at the mid-layer: inflammaging/leptin drive IL-6, TNF, NF-kB,', 10, None),
 ('   insulin resistance and myostatin; DR/metformin/resveratrol/canagliflozin drive AMPK; exercise cuts', 10, None),
 ('   myostatin and raises PGC1A/mito; leucine drives protein synthesis; NMN/NR raise NAD; the sign-fixed', 10, None),
 ('   adiponectin and omega-3 arms behave correctly (lower insulin resistance / IL-6).', 10, None),
 ('2. BOUNDARY: the top-level phenotypes Sarcopenia (range 0.11) and Catabolism (0.06) are compressed —', 10, 'C0392B'),
 ('   over-buffered by high basal rates, so the aggregation layer, not the mechanism layer, is the weak point.', 10, 'C0392B'),
 ('3. No dead levers: every input has downstream influence (weakest: GLP-1, glucocorticoids).', 10, None),
 ('4. Interventions are mostly SUB-ADDITIVE (diminishing returns when stacked on one axis, e.g. exercise+leucine', 10, None),
 ('   on anabolism); a few favourable synergies on sarcopenia (DR+leucine, leucine+NMN, NMN+GLP-1).', 10, None),
 ('5. Archetypes: exercise+leucine is the best anti-sarcopenia lever in every patient state; omega-3 is the', 10, None),
 ('   targeted IL-6 reducer; NMN restores depleted NAD in inflammaged/cirrhotic states; rapamycin costs lean mass.', 10, None),
 ('', 10, None),
 ('Tabs: Behaviour map · Single-input matrix · Dose-response · Interactions · Archetypes · Boundaries & extension', 10, '808080'),
 ('All values are last-state activity (0–1) from the 90-node cmaboss model; read them relatively. Evidence base: PubMed.', 10, '808080'),
]
for i, (t, sz, col) in enumerate(rows, 1):
    c = ws.cell(row=i, column=1, value=t); c.font = Font(bold=sz >= 12, size=sz, color=col or '000000'); c.alignment = WRAP

# 2. Behaviour map (image + boundary table)
ws = wb.create_sheet('Behaviour map')
ws['A1'] = 'Behaviour map — each lever in isolation, and dynamic range per readout'
ws['A1'].font = Font(bold=True, size=12, color='1F3864')
try:
    img = XLImage('fig_behaviour_map.png'); img.width = 980; img.height = 760
    ws.add_image(img, 'A3')
except Exception as e:
    ws['A3'] = f'[figure: fig_behaviour_map.png] {e}'

# 3. Single-input matrix
ws = wb.create_sheet('Single-input matrix')
A = O['single_input_delta']; base = O['baseline']
readouts = list(O['readout_dynamic_range'].keys())
inputs = list(O['input_total_influence'].keys())
ws.append(['input \\ readout (Δ vs baseline)'] + readouts)
hrow(ws, 1, len(readouts) + 1); ws.column_dimensions['A'].width = 18
for i in range(len(readouts)):
    ws.column_dimensions[get_column_letter(2 + i)].width = 7.5
for inp in inputs:
    ws.append([inp] + [A[inp][r] for r in readouts])
    rr = ws.max_row
    ws.cell(rr, 1).font = Font(bold=True)
    for ci, r in enumerate(readouts):
        cell = ws.cell(rr, 2 + ci); cell.alignment = CEN; cell.border = BORD
        f = heat(A[inp][r])
        if f and abs(A[inp][r]) >= 0.03:
            cell.fill = f
ws.append([]); ws.append(['baseline'] + [base[r] for r in readouts])
ws.freeze_panes = 'B2'

# 4. Dose-response
ws = wb.create_sheet('Dose-response')
ws['A1'] = 'Dose-response — fractional activation (0→1) of each lever'; ws['A1'].font = Font(bold=True, size=12, color='1F3864')
try:
    img2 = XLImage('fig_dose_response.png'); img2.width = 1000; img2.height = 460
    ws.add_image(img2, 'A3')
except Exception as e:
    ws['A3'] = f'[figure] {e}'

# 5. Interactions
ws = wb.create_sheet('Interactions')
ws.append(['Pairwise interactions (combined − additive prediction); negative = sub-additive/redundant'])
ws.cell(1, 1).font = Font(bold=True, size=12, color='1F3864'); ws.append([])
keys = ['Sarcopenia', 'Anabolism', 'Autophagy', 'NAD', 'Insulin_sensitivity', 'Mito_biogenesis']
head = ['pair'] + [f'{k} (int.)' for k in keys]
ws.append(head); hrow(ws, 3, len(head)); ws.column_dimensions['A'].width = 26
for i in range(len(keys)):
    ws.column_dimensions[get_column_letter(2 + i)].width = 16
for pair, d in O['pairwise_interaction'].items():
    ws.append([pair] + [d[k]['interaction'] for k in keys])
    rr = ws.max_row
    for ci, k in enumerate(keys):
        cell = ws.cell(rr, 2 + ci); cell.alignment = CEN; cell.border = BORD
        v = d[k]['interaction']
        if abs(v) >= 0.05:
            cell.fill = PatternFill('solid', fgColor='FBE5D6' if v < 0 else 'E2EFDA')
ws.freeze_panes = 'B4'

# 6. Archetypes
ws = wb.create_sheet('Archetypes')
ws.append(['Patient archetypes × interventions (activity 0–1)']); ws.cell(1, 1).font = Font(bold=True, size=12, color='1F3864')
ws.append([])
cols = ['Sarcopenia', 'Anabolism', 'MSTN', 'Insulin_resistance', 'IL6', 'NAD', 'Autophagy', 'Mito_biogenesis']
r = 3
for an, ivs in O['archetypes'].items():
    ws.cell(r, 1, an).font = Font(bold=True, size=11, color='1F3864'); r += 1
    ws.cell(r, 1, 'intervention');
    for ci, c in enumerate(cols):
        ws.cell(r, 2 + ci, c)
    hrow(ws, r, len(cols) + 1); r += 1
    base_iv = ivs['(none)']
    for ivn, vals in ivs.items():
        ws.cell(r, 1, ivn).alignment = TOP
        for ci, c in enumerate(cols):
            cell = ws.cell(r, 2 + ci, vals[c]); cell.alignment = CEN; cell.border = BORD
            d = vals[c] - base_iv[c]
            if ivn != '(none)' and abs(d) >= 0.03:
                good = (d < 0) == (c in ('Sarcopenia', 'MSTN', 'Insulin_resistance', 'IL6'))
                cell.fill = PatternFill('solid', fgColor='E2EFDA' if good else 'FBE5D6')
        r += 1
    r += 1
ws.column_dimensions['A'].width = 18
for i in range(len(cols)):
    ws.column_dimensions[get_column_letter(2 + i)].width = 11

# 7. Boundaries & extension
ws = wb.create_sheet('Boundaries & extension')
ws.column_dimensions['A'].width = 120
B = [
 ('What the model can do, and what it would need to do more', 13, '1F3864'),
 ('', 10, None),
 ('Can do now (clear signals)', 12, '1F3864'),
 ('• Inflammation / insulin-resistance axis: inflammaging, leptin, FGL1, hyperammonemia → IL-6, TNF, NF-kB, insulin resistance.', 10, None),
 ('• Nutrient-sensing / NAD axis: DR, rapamycin, metformin, canagliflozin, resveratrol, NMN/NR → AMPK, mTOR, NAD, SIRT1, PGC1A.', 10, None),
 ('• Muscle-mass axis: exercise, leucine, follistatin, decorin, irisin, myostatin → MSTN, protein synthesis, proteolysis, anabolism.', 10, None),
 ('• Liver→muscle bridge (Schneider): hyperammonemia/FGL1 → myostatin / insulin resistance.', 10, None),
 ('• Intervention comparison and ranking per patient archetype; dose-response and pairwise interaction structure.', 10, None),
 ('', 10, None),
 ('Current boundaries (honest limits)', 12, 'C0392B'),
 ('• Sarcopenia and Catabolism are compressed (range 0.11 / 0.06): the phenotype-aggregation layer is over-buffered', 10, None),
 ('  and needs recalibration before the top-level outcome is quantitatively trustworthy. Mechanism-layer readouts are fine.', 10, None),
 ('• Boolean/continuous steady state only: no real dose (fractional activation is a proxy), no time-to-effect, no dynamics.', 10, None),
 ('• Weights are hand-set, not fit to data: use the model for direction and ranking, not effect sizes.', 10, None),
 ('• No methylation layer yet: clocks, EpiScores, cell composition and TFBS-CpGs are not wired as inputs/readouts.', 10, None),
 ('• DR lean-mass effect is small and feedback-gated (S6K1⊣IRS1); it is a hypothesis, not a calibrated prediction.', 10, None),
 ('', 10, None),
 ('Methylation personalization — how a single patient would parameterize the model', 12, '1F3864'),
 ('Model input / node            ←  blood-methylation feature that would set it', 11, '1F3864'),
 ('Inflammaging / NF-kB / IL-6   ←  inflammatory EpiScores (CRP, IL-6), immune cell-type composition (CD8/naive)', 10, None),
 ('Insulin resistance / sensitivity ←  insulin-resistance EpiScore, HbA1c/HOMA proxy methylation scores', 10, None),
 ('hyperammonemia / FGL1 / LECT2 ←  liver-function / hepatokine EpiScores (the Schneider axis)', 10, None),
 ('NAD / SIRT1 / mitochondria    ←  age-related NAD-pathway and mitochondrial EpiScores', 10, None),
 ('Exercise / DR / omega-3 / drugs ←  lifestyle & exposure scores + the patient’s actual regimen (as inputs)', 10, None),
 ('Clocks (Horvath, TFMethyl)    ←  READOUTS: compare model-predicted aging state to measured epigenetic age', 10, None),
 ('', 10, None),
 ('Minimum extension to reach the original goal (methylation → per-patient simulation)', 12, '1F3864'),
 ('1. Add a methylation input layer: map EpiScores/cell-composition/TFBS-CpGs to the node initial states above.', 10, None),
 ('2. Recalibrate the phenotype-aggregation layer so Sarcopenia/Catabolism have usable dynamic range.', 10, None),
 ('3. Add clocks as explicit readouts and fit a thin ML layer from node states → clock / outcome.', 10, None),
 ('4. Calibrate key weights (DR vs rapamycin mTOR tone; feedback strength) against human data.', 10, None),
]
for i, (t, sz, col) in enumerate(B, 1):
    c = ws.cell(row=i, column=1, value=t); c.font = Font(bold=sz >= 12, size=sz, color=col or '000000'); c.alignment = WRAP

out = 'MUS_P1_model_characterization.xlsx'
wb.save(out); print('saved', out, wb.sheetnames)

# ===== appended: aging-layer tabs =====
from openpyxl import load_workbook
AL=json.load(open('aging_layer.json'))
wb2=load_workbook(out)
SUBN={'Aging_nutrient_sensing':'Nutrient-sensing','Aging_mito':'Mitochondrial','Aging_NAD':'NAD decline',
 'Aging_inflammation':'Inflammaging','Aging_senescence_genomic':'Senescence/genomic',
 'Aging_proteostasis':'Proteostasis','Aging_secretory':'Secretory (GDF15)'}
SUB=AL['subscores']
# Aging layer overview
ws=wb2.create_sheet('Aging layer')
L=[('Aging layer — hallmarks of aging as mechanistic sub-scores (no clocks → no circularity)',13,'1F3864'),
 ('',10,None),
 ('Seven hallmark sub-scores + an integrated index are read out from pathway activity and secretory markers',10,None),
 ('already in the map, NOT from epigenetic clocks. Clocks stay an external comparison; the ML layer predicts',10,None),
 ('clock/outcome from node states, so agreement with measured epigenetic age is a result, not an assumption.',10,None),
 ('',10,None),
 ('Sub-score   →  model nodes (direction: higher sub-score = more aged)   |  key literature (PubMed)',11,'1F3864'),
 ('Nutrient-sensing   ↑mTOR, ↑IGF1, ↑insulin resistance, ↓AMPK            |  López-Otín 2013 (23746838)',10,None),
 ('Mitochondrial      ↑Mito_dysfunction, ↑ROS, ↓PGC1A/Mito_biogenesis     |  López-Otín 2013/2023',10,None),
 ('NAD decline        ↑CD38, ↓NAD                                          |  Camacho-Pereira 2016 (27304511)',10,None),
 ('Inflammaging       ↑NF-κB, ↑IL-6, ↑TNF                                  |  Varadhan 2014 (23689826)',10,None),
 ('Senescence/genomic ↑p53, ↑PARP1, ↑DNA damage                           |  López-Otín 2013 (23746838)',10,None),
 ('Proteostasis       ↑proteolysis, ↓autophagy, ↓protein synthesis        |  López-Otín 2023 (36599349)',10,None),
 ('Secretory (GDF15)  ↑GDF15, ↑myostatin                                   |  Wiklund 2010 (20854422)',10,None),
 ('',10,None),
 ('DOIs: 10.1016/j.cell.2013.05.039 · 10.1016/j.cell.2022.11.001 · 10.1016/j.cmet.2016.05.006 ·',9,'808080'),
 ('10.1093/gerona/glt023 · 10.1111/j.1474-9726.2010.00629.x',9,'808080')]
ws.column_dimensions['A'].width=118
for i,(t,sz,c) in enumerate(L,1):
    cc=ws.cell(i,1,t); cc.font=Font(bold=sz>=12,size=sz,color=c or '000000'); cc.alignment=WRAP
try:
    im=XLImage('fig_aging_interventions.png'); im.width=720; im.height=560; ws.add_image(im,'A20')
    im2=XLImage('fig_aging_archetypes.png'); im2.width=720; im2.height=340; ws.add_image(im2,'A52')
except Exception as e: ws['A20']=str(e)
# Archetype profile table
ws=wb2.create_sheet('Aging by archetype')
ws.append(['Aging sub-score profile per patient archetype (higher = more aged on that hallmark)'])
ws.cell(1,1).font=Font(bold=True,size=12,color='1F3864'); ws.append([])
head=['archetype','Aging index']+[SUBN[s] for s in SUB]; ws.append(head); hrow(ws,3,len(head))
ws.column_dimensions['A'].width=20
for s in range(len(head)-1): ws.column_dimensions[get_column_letter(2+s)].width=14
for an,prof in AL['archetype_profile'].items():
    ws.append([an,prof['Aging_index']]+[prof[s] for s in SUB]); rr=ws.max_row
    for ci in range(1,len(head)):
        cell=ws.cell(rr,1+ci); cell.alignment=CEN; cell.border=BORD
        v=cell.value
        if isinstance(v,(int,float)): cell.fill=heat((v-0.4)*1.5)
ws.freeze_panes='B4'
# Intervention delta table
ws=wb2.create_sheet('Aging by intervention')
ws.append(['Intervention effect on aging sub-scores (Δ vs healthy-fed; negative/blue = lowers aging)'])
ws.cell(1,1).font=Font(bold=True,size=12,color='1F3864'); ws.append([])
head=['intervention','Δ Aging index']+[SUBN[s] for s in SUB]; ws.append(head); hrow(ws,3,len(head))
ws.column_dimensions['A'].width=16
for s in range(len(head)-1): ws.column_dimensions[get_column_letter(2+s)].width=14
for iv,d in AL['intervention_delta'].items():
    ws.append([iv,d['Aging_index']]+[d[s] for s in SUB]); rr=ws.max_row
    for ci in range(1,len(head)):
        cell=ws.cell(rr,1+ci); cell.alignment=CEN; cell.border=BORD
        v=cell.value
        if isinstance(v,(int,float)) and abs(v)>=0.03: cell.fill=heat(v)
ws.freeze_panes='B4'
wb2.save(out); print('aging tabs added ->',out, wb2.sheetnames)
