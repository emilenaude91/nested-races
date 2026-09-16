"""Which laboratory prefers which regulatory regime? Static period-2 payoffs from the lead model (leader-link, winner-causes accounting)."""
import numpy as np
from lead4c import P, vals2
print("delta | regime | labU | labC | piU | Q   (leader U, follower C)")
for noise in [0.1,0.3]:
    for delta in [0.0,0.3,0.6,0.9]:
        p=P(noise=noise,ell=60.0,delta=delta); rows={}
        for reg in [(1,1),(1,0),(0,1),(0,0)]:
            v=vals2(p,reg,delta); rows[reg]=v
        lab={(1,1):'RR',(1,0):'RA',(0,1):'AR',(0,0):'AA'}
        best_U=max(rows,key=lambda r: rows[r][2]); best_C=max(rows,key=lambda r: rows[r][3])
        print(f"noise={noise} delta={delta}: " + " | ".join(f"{lab[r]}: U={rows[r][2]:+.2f} C={rows[r][3]:+.2f} piU={rows[r][5]:.2f} Q={rows[r][4]:.3f}" for r in rows) + f"  => leader prefers {lab[best_U]}, follower prefers {lab[best_C]}; U: RR-AA={rows[(1,1)][2]-rows[(0,0)][2]:+.2f}, C: RR-AA={rows[(1,1)][3]-rows[(0,0)][3]:+.2f}")
# average-risk symmetric baseline: lab payoff RR vs AA as function of lambda
from core4 import Par, lab_eq, Q, pis, qj
print("\nsymmetric average risk: lab payoff RR vs AA")
for lam in [0.05,0.2,0.5]:
    for ell in [30,100]:
        p=Par(lam=lam,ell=ell); out={}
        for reg in [(1,1),(0,0)]:
            s,_=lab_eq(p,p.m*reg[0],p.m*reg[1]); qq=Q(s[0],s[1],p); pU=pis(s[0],s[1],p)
            u=(1-qq)*pU*p.WL - p.kappa*s[0]**2/2 - p.lam*qq*p.ell - p.rho*p.m*reg[0]*(1-s[0]); out[reg]=u
        print(f"lam={lam} ell={ell}: u_lab(RR)={out[(1,1)]:+.2f} u_lab(AA)={out[(0,0)]:+.2f} -> labs prefer {'RR' if out[(1,1)]>out[(0,0)] else 'AA'}")
