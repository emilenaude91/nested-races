import numpy as np
from fast import Par, lab_eq, Q, pi_all
def payoff(p, rrow, rcol):
    e=np.array([p.m*rrow, p.m*rcol]); s=lab_eq(p,e); q=Q(s,p); P=pi_all(s,p)
    return (1-q)*P[0]*p.W - p.Lam*q*p.ell - p.d*rrow, q, s
R=np.linspace(0,1,11)
def sustainable(p):
    U=np.zeros((11,11)); Qm=np.zeros((11,11))
    for i,rr in enumerate(R):
        for j,rc in enumerate(R):
            U[i,j],Qm[i,j],_=payoff(p,rr,rc)
    sym=[]
    for j,rc in enumerate(R):
        if U[j,j]>=U[:,j].max()-1e-9: sym.append(rc)
    return sym, U, Qm
for g in [1.0,2.0]:
    for ell in [15,25,30,40]:
        for m in [0.5,1.0]:
            p=Par(gamma=g,ell=ell,m=m); sym,U,Qm=sustainable(p)
            rmax=max(sym) if sym else None
            qr = Qm[int(round(rmax*10)),int(round(rmax*10))] if rmax is not None else None
            print(f"gamma={g} ell={ell} m={m}: symmetric NE r-levels={ [round(x,1) for x in sym] } r_max={rmax} Q(r_max)={qr if qr is None else round(qr,3)} Q(0)={Qm[0,0]:.3f} Q(1)={Qm[10,10]:.3f}")
