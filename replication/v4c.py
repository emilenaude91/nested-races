from core4 import *
import sys, time
R=np.linspace(0,1,21)
def payoff(p,rrow,rcol):
    s,_=lab_eq(p,p.m*rrow,p.m*rcol); qq=Q(s[0],s[1],p)
    return prizeU(s[0],s[1],p)*p.W - p.Lam*qq*p.ell - p.d*rrow
def ladder(p):
    U=np.zeros((21,21))
    for i,rr in enumerate(R):
        for j,rc in enumerate(R): U[i,j]=payoff(p,rr,rc)
    sym=[]; gains=[]
    for j in range(21):
        g=U[:,j].max()-U[j,j]; gains.append(g)
        if g<=1e-6: sym.append(round(float(R[j]),2))
    return sym, gains
t=time.time()
for tech in ['avg','wl']:
    for m in [0.5,0.75,1.0]:
        for ell in ([25,30,35,40,45,50] if tech=='avg' else [10,20,30,60,100]):
            p=Par(tech=tech,m=m,ell=ell); sym,g=ladder(p)
            # also smallest deviation gain among non-equilibrium interior rungs
            interior=[(round(float(R[j]),2),round(g[j],3)) for j in range(1,20) if g[j]>1e-6]
            mn=min(interior,key=lambda x:x[1]) if interior else None
            print(f"{tech} m={m} ell={ell}: symmetric NE r = {sym}; smallest deviation gain at a non-NE interior rung {mn}")
            sys.stdout.flush()
print("t",round(time.time()-t))
