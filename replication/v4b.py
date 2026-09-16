from core4 import *
import json, time, sys
t=time.time()
for tech in ['wl','ll']:
    f=lambda x: (lambda c: c['DR']+c['DA'])(cells(Par(tech=tech,ell=x)))
    print(tech,"ell_rd=",round(brentq(f,5,150,xtol=0.05),1)); sys.stdout.flush()
# decisiveness under corrected ll
for r in [1.0,2.0,3.0]:
    for ell in [30,60]:
        c=cells(Par(tech='ll',r=r,ell=ell)); print(f"ll r={r} ell={ell}: {c['kind']} DR={c['DR']:+.2f} DA={c['DA']:+.2f} sRR={c['sRR'][0]:.2f} sAR=({c['sAR'][0]:.2f},{c['sAR'][1]:.2f})")
sys.stdout.flush()
# belief map
ratios=np.geomspace(0.02,2.0,14); q0s=np.linspace(0.02,0.7,8); res={}
kmap={'PD':0,'Deadlock':0,'Chicken':1,'SH':2,'Rdom':3}
for tech in ['avg','wl','ll']:
    M=np.zeros((len(q0s),len(ratios)),int)
    for i,q0 in enumerate(q0s):
        for j,ra in enumerate(ratios):
            M[i,j]=kmap[cells(Par(tech=tech,W=100*ra,ell=100.0,q0=q0))['kind']]
    res[tech]=M.tolist(); print(tech,"map done",round(time.time()-t)); sys.stdout.flush()
json.dump(dict(ratios=ratios.tolist(),q0s=q0s.tolist(),res=res),open('wmap4.json','w'))
