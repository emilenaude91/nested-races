from core4 import *
import sys
R=np.linspace(0,1,11)
def ladder(p):
    U=np.zeros((11,11))
    for i,rr in enumerate(R):
        for j,rc in enumerate(R):
            s,_=lab_eq(p,p.m*rr,p.m*rc); qq=Q(s[0],s[1],p); U[i,j]=prizeU(s[0],s[1],p)*p.W-p.Lam*qq*p.ell-p.d*rr
    sym=[round(float(R[j]),1) for j in range(11) if U[:,j].max()-U[j,j]<=1e-6]
    gains=[round(float(U[:,j].max()-U[j,j]),3) for j in range(11)]
    return sym,gains
for m in [1.0,0.5]:
    for ell in [20,40,60,100]:
        sym,g=ladder(Par(tech='wl',m=m,ell=ell)); print(f"wl m={m} ell={ell}: symmetric NE r={sym}; deviation gains by rung={g}"); sys.stdout.flush()
