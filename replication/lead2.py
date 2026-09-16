# Period-2 value functions on a delta grid, for both regimes and both techs; saved for the dynamic model
import numpy as np, json, time, sys
from lead import Par, profile
deltas=np.round(np.arange(-0.6,1.41,0.1),2)
out={}
t0=time.time()
for tech in ['ll','avg']:
    for noise in [0.1]:
        for r in [(1,1),(0,0),(1,0),(0,1)]:
            key=f"{tech}|{noise}|{r[0]}{r[1]}"; rows=[]
            for dl in deltas:
                p=Par(delta=float(dl),noise=noise,ell=60.0,tech=tech); P=profile(p,r)
                # lab period-2 payoffs
                s=P['s']; qq=P['Q']
                from lead import pis, u
                labU = u(s[0],s,p,0,p.m*r[0]); labC = u(s[1],s,p,1,p.m*r[1])
                rows.append([float(dl),P['V'][0],P['V'][1],labU,labC,qq,P['piU'],float(s[0]),float(s[1])])
            out[key]=rows; print(key,round(time.time()-t0)); sys.stdout.flush()
json.dump(dict(deltas=deltas.tolist(),vals=out),open('lead_v2.json','w'))
