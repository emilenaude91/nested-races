import numpy as np, json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from lead5 import run
d0s=np.arange(0.1,0.91,0.1); rec={}
for theta in [0.0,0.2]:
    for reg in [(1,1),(1,0)]:
        rows=[]
        for d0 in d0s:
            P,ne=run(0.1,60.0,1.0,theta,float(d0)); r=P[reg]
            rows.append([d0,r['s2'][0],r['piU2'],r['Q2'],r['d1'],r['s1'][0]])
        rec[f"{theta}|{reg[0]}{reg[1]}"]=rows
json.dump(rec,open('leadfig.json','w'))
fig,axes=plt.subplots(1,3,figsize=(9.2,3.1))
sty={'0.0|11':('k-','both regulate, no erosion'),'0.0|10':('k--','leader alone, no erosion'),'0.2|11':('k-','both regulate, erosion 0.2'),'0.2|10':('k--','leader alone, erosion 0.2')}
for k,(ls,lab) in sty.items():
    A=np.array(rec[k]); lw=1.6 if k.startswith('0.0') else 0.9
    axes[0].plot(A[:,0],A[:,1],ls,lw=lw,label=lab); axes[1].plot(A[:,0],A[:,2],ls,lw=lw); axes[2].plot(A[:,0],A[:,3],ls,lw=lw)
axes[0].set_title("leader's safety, period 2",fontsize=10); axes[1].set_title("leader's win probability",fontsize=10); axes[2].set_title("global catastrophe risk",fontsize=10)
for ax in axes: ax.set_xlabel(r"initial lead $\delta_0$")
axes[0].legend(fontsize=7,frameon=False,loc='upper left'); axes[1].set_ylim(0,1.05); axes[2].set_ylim(0.05,0.21)
plt.tight_layout(); plt.savefig('paper/fig_lead.pdf'); print("ok")
