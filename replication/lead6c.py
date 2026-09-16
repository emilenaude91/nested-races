"""Check whether the period-2 regime game ever lacks a pure NE on the lead grid; implement mixed fallback and recompute the re-choice results."""
import numpy as np, sys
from lead4c import P, vals2, DG, eq1
REGS=[(1,1),(1,0),(0,1),(0,0)]
alt={ (1,1):((0,1),(1,0)), (1,0):((0,0),(1,1)), (0,1):((1,1),(0,0)), (0,0):((1,0),(0,1)) }
def period2(p,delta):
    V={reg:vals2(p,reg,delta) for reg in REGS}
    ne=[reg for reg in REGS if V[reg][0]>=V[alt[reg][0]][0]-1e-9 and V[reg][1]>=V[alt[reg][1]][1]-1e-9]
    if ne:
        sel=max(ne,key=lambda r: V[r][0]+V[r][1]); return V[sel],('pure',sel),ne
    # mixed NE of 2x2: U regulates w.p. pU making C indifferent; C regulates w.p. pC making U indifferent
    # C's payoff: V[(rU,rC)][1]; indifference for C: pU*V[(1,1)][1]+(1-pU)*V[(0,1)][1] = pU*V[(1,0)][1]+(1-pU)*V[(0,0)][1]
    a=V[(1,1)][1]-V[(1,0)][1]; b=V[(0,0)][1]-V[(0,1)][1]; pU=b/(a+b)
    c=V[(1,1)][0]-V[(0,1)][0]; d=V[(0,0)][0]-V[(1,0)][0]; pC=d/(c+d)
    w={(1,1):pU*pC,(1,0):pU*(1-pC),(0,1):(1-pU)*pC,(0,0):(1-pU)*(1-pC)}
    EV=sum(w[r]*V[r] for r in REGS)
    return EV,('mixed',(round(pU,3),round(pC,3))),ne
p=P(noise=0.1,ell=60.0)
n_mixed=0; tab=[]
for dl in DG:
    v,kind,ne=period2(p,float(dl)); tab.append(v)
    if kind[0]=='mixed': n_mixed+=1; print(f"delta={dl:+.1f}: NO pure NE; mixed pU,pC={kind[1]}")
print("grid points with no pure period-2 equilibrium:",n_mixed,"of",len(DG))
# also check on a finer grid
fine=np.round(np.arange(-0.8,1.61,0.02),2); nm=0; where=[]
for dl in fine:
    _,kind,_=period2(p,float(dl))
    if kind[0]=='mixed': nm+=1; where.append(float(dl))
print("fine grid (step 0.02): no-pure-NE points:",nm, where[:20])
