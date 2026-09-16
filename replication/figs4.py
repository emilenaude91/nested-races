import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Rectangle, Patch
from core4 import *
# Fig 2: D_R, D_A vs ell for avg and wl
fig,axes=plt.subplots(1,2,figsize=(9,3.3))
for ax,tech,ttl in zip(axes,['avg','wl'],['average risk','weakest-link risk']):
    ells=np.linspace(2,200,60); DR=[];DA=[]
    for l in ells:
        c=cells(Par(tech=tech,ell=l)); DR.append(c['DR']); DA.append(c['DA'])
    ax.plot(ells,DR,'k-',label=r"$D_R=u_{RR}-u_{AR}$"); ax.plot(ells,DA,'k--',label=r"$D_A=u_{RA}-u_{AA}$")
    ax.axhline(0,color='grey',lw=0.7); ax.set_xlabel(r"shared catastrophe loss $\ell$"); ax.set_title(ttl,fontsize=10); ax.set_ylim(-2,8); ax.legend(fontsize=8,loc='upper left')
axes[0].set_ylabel("gain from Regulate"); plt.tight_layout(); plt.savefig('paper/fig_taxonomy.pdf')
# Fig 3: belief map
d=json.load(open('wmap4.json')); ratios=np.array(d['ratios']); q0s=np.array(d['q0s'])
cmap=ListedColormap(['#F0997B','#FAC775','#AFA9EC','#9FE1CB']); norm=BoundaryNorm([-.5,.5,1.5,2.5,3.5],4)
fig,axes=plt.subplots(1,3,figsize=(9.2,3.5),sharey=True)
titles={'avg':'average risk','wl':'weakest-link risk','ll':'leader-link risk'}
xe=np.geomspace(ratios[0]/1.2,ratios[-1]*1.2,len(ratios)+1); ye=np.linspace(q0s[0]-0.045,q0s[-1]+0.045,len(q0s)+1)
for ax,tech in zip(axes,['avg','wl','ll']):
    M=np.array(d['res'][tech]); ax.pcolormesh(xe,ye,M,cmap=cmap,norm=norm,shading='flat'); ax.set_xscale('log'); ax.set_title(titles[tech],fontsize=10); ax.set_xlabel(r'prize / catastrophe loss, $W/\ell$')
    ax.add_patch(Rectangle((0.5,0.0),1.9,0.2,fill=False,ls='--',lw=1,ec='k')); ax.text(0.52,0.02,'"durable\nhegemony"',fontsize=7,va='bottom')
    ax.add_patch(Rectangle((0.02,0.2),0.13,0.5,fill=False,ls='--',lw=1,ec='k')); ax.text(0.024,0.72,'"temporary lead,\nhigh risk"',fontsize=7,va='bottom')
    ax.add_patch(Rectangle((0.5,0.3),1.9,0.4,fill=False,ls=':',lw=1.2,ec='k')); ax.text(0.52,0.72,'"pace the\nfrontier"',fontsize=7,va='bottom'); ax.set_ylim(0,0.85)
axes[0].set_ylabel(r'baseline catastrophe probability, $q_0$')
fig.legend(handles=[Patch(color='#F0997B',label="Prisoner's Dilemma"),Patch(color='#FAC775',label='Chicken'),Patch(color='#AFA9EC',label='Stag Hunt'),Patch(color='#9FE1CB',label='Regulate dominant')],loc='lower center',ncol=4,fontsize=8,frameon=False,bbox_to_anchor=(0.5,-0.02))
plt.tight_layout(rect=(0,0.06,1,1)); plt.savefig('paper/fig_beliefs.pdf')
# Fig lead
R=json.load(open('lead_v4.json'))['0.1|60|0.0']
d0=[r['d0'] for r in R]; 
fig,axes=plt.subplots(1,3,figsize=(9.2,3.0))
axes[0].plot(d0,[r['RR']['sU2'] for r in R],'k-',label='both states regulate'); axes[0].plot(d0,[r['RA']['sU1'] for r in R],'k--',label='leader regulates alone')
axes[1].plot(d0,[r['RR']['piU2'] for r in R],'k-'); axes[1].plot(d0,[r['RA']['piU2'] for r in R],'k--')
axes[2].plot(d0,[r['RR']['Q2'] for r in R],'k-'); axes[2].plot(d0,[r['RA']['Q2'] for r in R],'k--')
axes[0].set_title("leader's safety",fontsize=10); axes[1].set_title("leader's win probability, period 2",fontsize=10); axes[2].set_title("global catastrophe risk, period 2",fontsize=10)
for ax in axes: ax.set_xlabel(r"initial lead $\delta_0$")
axes[0].legend(fontsize=8,frameon=False,loc='upper left'); axes[0].set_ylim(0,1); axes[1].set_ylim(0,1.05); axes[2].set_ylim(0.05,0.21)
plt.tight_layout(); plt.savefig('paper/fig_lead.pdf'); print('figs ok')
