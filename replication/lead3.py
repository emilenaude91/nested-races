# Two-period lead-as-stock model. States commit to a regime r=(rU,rC) for both periods.
# Period 1: labs choose s1; lead evolves d1 = d0 + tau*(xU1 - xC1) - theta; period-1 risk Q1 (leader-link on current lead).
# Period 2: static game at d1 (precomputed). Backward induction; period-1 lab equilibrium by grid BR.
import numpy as np, json, sys
from lead import Par, pis, q, Q
D=json.load(open('lead_v2.json')); deltas=np.array(D['deltas'])
def cont(tech,r,dl):
    rows=np.array(D['vals'][f"{tech}|0.1|{r[0]}{r[1]}"]); dl=np.clip(dl,deltas[0],deltas[-1])
    return [np.interp(dl,rows[:,0],rows[:,k]) for k in (1,2,3,4,5,6,7,8)]  # VU,VC,labU,labC,Q2,piU2,sU2,sC2
def lead_next(s1,p):
    x=1-np.asarray(s1)+p.eps; return p.delta + p.tau*(x[0]-x[1]) - p.theta
def u1(sj,s1,p,j,e,r):
    s=np.array(s1,float); s[j]=sj; q1=Q(s,p); d1=lead_next(s,p); c=cont(p.tech,r,d1)
    return -p.kappa*sj**2/2 - p.lam*q1*p.ell - p.rho*e*(1-sj) + (1-q1)*c[2+j]
GR=np.linspace(0,1,101)
def br1(s1,p,j,e,r):
    vals=[u1(g,s1,p,j,e,r) for g in GR]; return GR[int(np.argmax(vals))]
def eq1(p,r):
    e=(p.m*r[0],p.m*r[1]); s=np.array([0.3,0.3])
    for it in range(200):
        sn=np.array([br1(s,p,j,e[j],r) for j in range(2)])
        if np.max(np.abs(sn-s))<1e-9: break
        s=0.5*s+0.5*sn
    s=np.round(s*100)/100
    q1=Q(s,p); d1=lead_next(s,p); c=cont(p.tech,r,d1)
    VU=-p.Lam*q1*p.ell - p.d*r[0] + (1-q1)*c[0]; VC=-p.Lam*q1*p.ell - p.d*r[1] + (1-q1)*c[1]
    return dict(s1=s,Q1=q1,d1=d1,VU=VU,VC=VC,Q2=c[4],piU2=c[5],s2=(c[6],c[7]))
def game1(p):
    P={r:eq1(p,r) for r in [(1,1),(1,0),(0,1),(0,0)]}; ne=[]
    for rU in (0,1):
        for rC in (0,1):
            if P[(rU,rC)]['VU']>=P[(1-rU,rC)]['VU']-1e-9 and P[(rU,rC)]['VC']>=P[(rU,1-rC)]['VC']-1e-9: ne.append((rU,rC))
    return P,ne
if __name__=="__main__":
    for tech in ['ll','avg']:
        print(f"\n===== tech={tech}, noise=0.1, ell=60, tau=1 =====")
        for theta in [0.0,0.2,0.4]:
            for d0 in [0.3,0.6]:
                p=Par(delta=d0,noise=0.1,ell=60.0,tech=tech); p.tau=1.0; p.theta=theta
                P,ne=game1(p)
                f=lambda r: f"s1=({P[r]['s1'][0]:.2f},{P[r]['s1'][1]:.2f}) d1={P[r]['d1']:+.2f} s2=({P[r]['s2'][0]:.2f},{P[r]['s2'][1]:.2f}) piU2={P[r]['piU2']:.2f} Q1={P[r]['Q1']:.3f} Q2={P[r]['Q2']:.3f}"
                print(f"theta={theta} d0={d0} NE={ne}\n   RA(U paces): {f((1,0))}\n   AA        : {f((0,0))}\n   RR        : {f((1,1))}")
                sys.stdout.flush()
