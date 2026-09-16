# Robustness of period-2 value functions: ll tech, (noise, ell) combos; appends to lead_v2.json
import numpy as np, json, time, sys
from lead import Par, profile, u
D=json.load(open('lead_v2.json')); deltas=np.array(D['deltas']); out=D['vals']; t0=time.time()
combos=[(0.3,60.0),(0.1,40.0),(0.1,100.0)]
for noise,ell in combos:
    for r in [(1,1),(0,0),(1,0),(0,1)]:
        key=f"ll|{noise}|{ell:g}|{r[0]}{r[1]}"
        if key in out: continue
        rows=[]
        for dl in deltas:
            p=Par(delta=float(dl),noise=noise,ell=ell,tech='ll'); P=profile(p,r); s=P['s']
            rows.append([float(dl),P['V'][0],P['V'][1],u(s[0],s,p,0,p.m*r[0]),u(s[1],s,p,1,p.m*r[1]),P['Q'],P['piU'],float(s[0]),float(s[1])])
        out[key]=rows; json.dump(dict(deltas=deltas.tolist(),vals=out),open('lead_v2.json','w'))
        print(key,round(time.time()-t0)); sys.stdout.flush()
