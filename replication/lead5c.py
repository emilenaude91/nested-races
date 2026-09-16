"""Period-2 regime re-choice: states re-choose r in period 2 given the realised lead; laboratories anticipate it."""
import numpy as np, json, sys
from lead4c import P, vals2, DG, Qll, eq1, cont, G1, u1_vec
REGS=[(1,1),(1,0),(0,1),(0,0)]
def static_ne(p,delta):
    V={reg:vals2(p,reg,delta) for reg in REGS}
    alt={ (1,1):((0,1),(1,0)), (1,0):((0,0),(1,1)), (0,1):((1,1),(0,0)), (0,0):((1,0),(0,1)) }
    ne=[reg for reg in REGS if V[reg][0]>=V[alt[reg][0]][0]-1e-9 and V[reg][1]>=V[alt[reg][1]][1]-1e-9]
    # selection: Pareto-best pure NE for states (sum of V), else (0,0)
    if ne: sel=max(ne,key=lambda r: V[r][0]+V[r][1])
    else: sel=(0,0)
    return sel,V[sel],ne
def build_rechoice_table(p):
    T=[];regs=[]
    for dl in DG:
        sel,v,ne=static_ne(p,float(dl)); T.append(v); regs.append(sel)
    return np.array(T),regs
def game1_rechoice(p,T,prev):
    Tdict={reg:T for reg in REGS}   # same continuation for all period-1 regimes (period-2 regime is re-chosen)
    Pm={reg:eq1(p,reg,Tdict,prev.get(reg,(0.3,0.3))) for reg in REGS}
    alt={ (1,1):((0,1),(1,0)), (1,0):((0,0),(1,1)), (0,1):((1,1),(0,0)), (0,0):((1,0),(0,1)) }
    ne=[reg for reg in REGS if Pm[reg]['VU']>=Pm[alt[reg][0]]['VU']-1e-9 and Pm[reg]['VC']>=Pm[alt[reg][1]]['VC']-1e-9]
    return Pm,ne
lab={(1,1):'RR',(1,0):'RA',(0,1):'AR',(0,0):'AA'}
for noise,ell in [(0.1,60.0)]:
    p=P(noise=noise,ell=ell); T,regs=build_rechoice_table(p)
    print("period-2 regime chosen by lead:", " ".join(f"{d:+.1f}:{lab[r]}" for d,r in zip(DG,regs) if abs(d*10-round(d*10))<1e-6 and -0.4<=d<=1.6 and (round(d*10)%2==0)))
    for theta in [0.0,0.2,0.4]:
        prev={}
        for d0 in [0.2,0.3,0.4,0.6,0.9]:
            p=P(noise=noise,ell=ell,theta=theta,delta=d0); Pm,ne=game1_rechoice(p,T,prev); prev={k:v['s1'] for k,v in Pm.items()}
            rr=Pm[(1,1)]; d1=rr['d1']; sel2=regs[int(np.argmin(np.abs(DG-np.clip(d1,DG[0],DG[-1]))))]
            print(f"theta={theta} d0={d0}: period-1 NE={[lab[n] for n in ne]} | under RR1: s1=({rr['s1'][0]:.2f},{rr['s1'][1]:.2f}) d1={d1:+.2f} -> period-2 regime {lab[sel2]}, piU2={rr['piU2']:.2f} Q2={rr['Q2']:.3f}")
            sys.stdout.flush()
