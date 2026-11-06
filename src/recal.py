"""Recalibration harness: patch the aggregation stack, report anchor scenarios,
dynamic range and sign contrasts. Edit STACK to iterate."""
import cmaboss_128n as cmaboss, copy, tempfile, shutil, sys
import coupled_model as CM
INPUTS=CM.INPUTS; NODES=[x for x in CM.N if x not in INPUTS]
ANCH={'Healthy-fasted':{}, 'Healthy-fed':{'Nutrients':1},
      'Aged':{'Inflammaging':1,'hyperammonemia':1,'DNA_damage':1},
      'Cachectic':{'Inflammaging':1,'hyperammonemia':1,'glucocorticoids':1},
      'Ex+Leu (fed)':{'Nutrients':1,'Exercise':1,'leucine':1},
      'Obese':{'LEP':1,'FGL1':1,'Nutrients':1}}
SHOW=['Sarcopenia','Catabolism','Proteolysis','Anabolism','Protein_synthesis','FOXO3','UPS','MSTN']
def runN(N,fix,sc=800):
    d=tempfile.mkdtemp()
    open(d+'/m.bnd','w').write('\n'.join([f'Node {n} {{ logic={n}; rate_up=0.0; rate_down=0.0; }}' for n in INPUTS]+
      [f'Node {n} {{ rate_up = {" + ".join(CM.P1.term(w,e) for w,e in up)}; rate_down = {" + ".join(CM.P1.term(w,e) for w,e in dn)}; }}' for n,(up,dn) in N.items() if n not in INPUTS]))
    cfg=['time_tick=0.5;','max_time=40;',f'sample_count={sc};','discrete_time=0;','use_physrandgen=0;','seed_pseudorandom=7;','thread_count=8;','statdist_traj_count=0;']
    for n in INPUTS: cfg.append(f'{n}.istate={fix.get(n,0)};')
    for n in NODES: cfg.append(f'{n}.istate=0;')
    for n in NODES: cfg.append(f'{n}.is_internal=0;')
    open(d+'/m.cfg','w').write('\n'.join(cfg))
    a,t,nodes,_=cmaboss.MaBoSSSim(network=d+'/m.bnd',config=d+'/m.cfg').run().get_last_nodes_probtraj()
    shutil.rmtree(d,ignore_errors=True)
    return {n:round(float(a[0][i]),3) for i,n in enumerate(nodes)}
def report(N,label):
    print(f'\n===== {label} =====')
    print('anchor'.ljust(16)+''.join(s[:10].rjust(11) for s in SHOW))
    rows={}
    for an,fix in ANCH.items():
        r=runN(N,fix); rows[an]=r
        print(an.ljust(16)+''.join(f'{r[s]:11.3f}' for s in SHOW))
    # dynamic range of Sarcopenia/Catabolism across anchors
    for k in ['Sarcopenia','Catabolism','Proteolysis']:
        vs=[rows[a][k] for a in ANCH]; print(f'  range {k}: {max(vs)-min(vs):.3f}')
    # sign checks
    hf=rows['Healthy-fed']['Sarcopenia']; exl=rows['Ex+Leu (fed)']['Sarcopenia']; ag=rows['Aged']['Sarcopenia']; ca=rows['Cachectic']['Sarcopenia']
    print(f'  Healthy-fed {hf} | Ex+Leu {exl} (want <hf) | Aged {ag} (want >hf) | Cachectic {ca} (want highest)')
    return rows

report(CM.N,'CURRENT')

def recalibrated():
    N=copy.deepcopy(CM.N)
    # bottom-up: lower resting of each stack node, balance formulations, matched gains
    up,dn=N['FOXO3']; N['FOXO3']=([(0.12,1)]+[(w,e) for w,e in up if e!=1],dn)       # drop high up-basal 0.5->0.12
    N['TRIM63']=([(0.08,1),(1,'NFkB'),(1,'FOXO3')],[(1.0,1)])                          # rest low
    N['UPS']=([(0.08,1),(1,'NFkB'),(1,'TRIM63'),(1,'smad2_3'),(1.5,'glucocorticoids')],[(1.1,1),(1,'IL15')])
    N['Proteolysis']=([(0.05,1),(1.8,'UPS'),(0.6,'Autophagy'),(1.3,'FOXO3')],[(1.4,1),(1.0,'FNDC5')])
    N['Anabolism']=([(0.15,1),(2.0,'Protein_synthesis')],[(0.35,1),(1.0,'Proteolysis')])
    N['Catabolism']=([(0.08,1),(2.0,'Proteolysis'),(1.0,'Mito_dysfunction')],[(0.4,1),(1.4,'Anabolism')])
    N['Sarcopenia']=([(0.07,1),(2.6,'Catabolism')],[(0.07,1),(2.6,'Anabolism')])
    return N

if 'recal' in sys.argv:
    report(recalibrated(),'RECALIBRATED v2')
