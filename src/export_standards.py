"""Export the coupled model to community standards for publication / FAIR reuse:
  - SBML-qual (SBML L3 qualitative models) — signed influence graph + canonical Boolean function
  - SIF + node/edge TSV — universal import (Cytoscape -> SBGN-AF via plugins)
  - SED-ML descriptor (template) + COMBINE/OMEX archive bundling everything with a manifest
Honest scope: SBML-qual carries the STRUCTURE (nodes, signed regulation) and a canonical sign-logic;
the EXACT executable dynamics (weighted MaBoSS rates) live in MUS_P1_coupled.bnd/.cfg, which are
bundled alongside. That two-file split (standard structural model + exact engine files) is normal and
is what the OMEX archive is for.
"""
import libsbml, os, zipfile, datetime
import coupled_model as CM

N, INPUTS = CM.N, set(CM.INPUTS)
READOUTS = set(CM.READOUTS)


def parse_regulators(node):
    """Return (activators, inhibitors) sets and (basal_up, basal_down) flags for a node's gates."""
    up, dn = N[node]
    act, inh = set(), set()
    basal_up = basal_dn = False
    def handle(terms, positive):
        nonlocal basal_up, basal_dn
        for w, e in terms:
            if e == 1:
                if positive: basal_up = True
                else: basal_dn = True
                continue
            parts = e.split(' & ') if ' & ' in str(e) else [e]
            for p in parts:
                neg = str(p).startswith('!')
                name = str(p)[1:] if neg else str(p)
                # positive gate: plain -> activator, !x -> inhibitor ; negative gate: reversed
                is_act = (positive and not neg) or (not positive and neg)
                (act if is_act else inh).add(name)
    handle(up, True); handle(dn, False)
    act.discard(1); inh.discard(1)
    return act, inh, basal_up, basal_dn


KNOWN = set(N) | INPUTS
# ---- collect signed edges ----
edges = []   # (regulator, sign, target)
for n in N:
    if n in INPUTS:
        continue
    act, inh, bu, bd = parse_regulators(n)
    for a in sorted(act):
        if a in KNOWN: edges.append((a, '+', n))
    for i in sorted(inh):
        if i in KNOWN: edges.append((i, '-', n))

all_nodes = list(CM.INPUTS) + list(N.keys())
print(f'{len(all_nodes)} nodes, {len(edges)} signed edges ({len(INPUTS)} inputs)')

# ---- SIF + TSV ----
with open('PLM_network.sif', 'w') as f:
    for s, sign, t in edges:
        f.write(f'{s}\t{"activates" if sign=="+" else "inhibits"}\t{t}\n')
with open('PLM_edges.tsv', 'w') as f:
    f.write('source\tsign\ttarget\n')
    for s, sign, t in edges:
        f.write(f'{s}\t{sign}\t{t}\n')
with open('PLM_nodes.tsv', 'w') as f:
    f.write('node\trole\n')
    for n in all_nodes:
        role = 'input' if n in INPUTS else ('readout' if n in READOUTS else 'internal')
        f.write(f'{n}\t{role}\n')

# ---- SBML-qual ----
NS = libsbml.SBMLNamespaces(3, 1, 'qual', 1)
doc = libsbml.SBMLDocument(NS)
doc.setPackageRequired('qual', True)
m = doc.createModel(); m.setId('Prevention_Longevity_Map_coupled')
m.setName('Prevention & Longevity Map — coupled P1+MUS Boolean model (influence graph)')
c = m.createCompartment(); c.setId('cell'); c.setConstant(True); c.setSpatialDimensions(3); c.setSize(1)
qm = m.getPlugin('qual')


def sid(x):
    import re
    s = re.sub(r'[^0-9a-zA-Z_]', '_', str(x))
    return s if s[0].isalpha() or s[0] == '_' else 's_' + s


idmap = {n: sid(n) for n in all_nodes}
for n in all_nodes:
    qs = qm.createQualitativeSpecies()
    qs.setId(idmap[n]); qs.setName(str(n)); qs.setCompartment('cell')
    qs.setConstant(n in INPUTS)                 # inputs are fixed (no transition updates them)
    qs.setMaxLevel(1); qs.setInitialLevel(0)

for n in all_nodes:
    if n in INPUTS:
        continue
    act, inh, bu, bd = parse_regulators(n)
    act = [a for a in act if a in KNOWN]; inh = [i for i in inh if i in KNOWN]
    tr = qm.createTransition(); tr.setId('tr_' + idmap[n])
    for r in sorted(set(act + inh)):
        inp = tr.createInput()
        inp.setId(f'in_{idmap[n]}_{idmap[r]}')
        inp.setQualitativeSpecies(idmap[r])
        inp.setTransitionEffect(libsbml.INPUT_TRANSITION_EFFECT_NONE)
        inp.setSign(libsbml.INPUT_SIGN_POSITIVE if r in act else libsbml.INPUT_SIGN_NEGATIVE)
    out = tr.createOutput(); out.setId(f'out_{idmap[n]}')
    out.setQualitativeSpecies(idmap[n]); out.setTransitionEffect(libsbml.OUTPUT_TRANSITION_EFFECT_ASSIGNMENT_LEVEL)
    # canonical sign-logic:  ON if (any activator) and not (any inhibitor); basal-driven if no activators
    def anyexpr(nodes):
        return '(' + ' || '.join(f'{idmap[x]} >= 1' for x in sorted(nodes)) + ')' if nodes else ''
    if act:
        formula = anyexpr(act) + (f' && !({anyexpr(inh)})' if inh else '')
        default = 0
    elif bu:                                    # constant drive, repressible
        formula = f'!({anyexpr(inh)})' if inh else '1'
        default = 0 if inh else 1
    else:
        formula = '0'; default = 0
    dt = tr.createDefaultTerm(); dt.setResultLevel(default)
    if formula not in ('0', '1') or (formula == '1' and default == 0):
        ft = tr.createFunctionTerm(); ft.setResultLevel(1 if formula != '0' else 0)
        ast = libsbml.parseL3Formula(formula if formula not in ('1',) else 'true')
        if ast is not None:
            ft.setMath(ast)

note = ('<body xmlns="http://www.w3.org/1999/xhtml"><p>Qualitative (SBML-qual) representation of the '
        'Prevention &amp; Longevity Map coupled P1+MUS Boolean model: nodes and SIGNED regulatory '
        'influences, with a canonical sign-consistent Boolean function (a node is active when an '
        'activator is active and no inhibitor is active; basal-driven nodes are active unless '
        'repressed). The EXACT executable dynamics use MaBoSS weighted transition rates and are in the '
        'accompanying MUS_P1_coupled.bnd / .cfg files. Evidence for each edge is PubMed-curated '
        '(see MUS_module_v3.json). Generated with libSBML ' + libsbml.getLibSBMLDottedVersion() + '.</p></body>')
m.setNotes(note)
libsbml.writeSBMLToFile(doc, 'PLM_model.sbml')
errs = doc.getNumErrors(libsbml.LIBSBML_SEV_ERROR)
doc.checkConsistency()
cerr = sum(doc.getError(i).getSeverity() == libsbml.LIBSBML_SEV_ERROR for i in range(doc.getNumErrors()))
print(f'SBML-qual written: PLM_model.sbml  (write errors={errs}, consistency errors={cerr})')
for i in range(doc.getNumErrors()):
    e = doc.getError(i)
    if e.getSeverity() == libsbml.LIBSBML_SEV_ERROR:
        print('  ERROR:', e.getMessage()[:140])
