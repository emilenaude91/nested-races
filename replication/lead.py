import numpy as np, json, sys, time
from scipy.optimize import brentq

class Par:
    def __init__(self, **kw):
        d = dict(W=10.0, sigma_share=0.5, kappa=6.0, q0=0.2, alpha=0.8, eps=0.05, tech='ll',
                 lam=0.05, Lam=1.0, rho=4.0, m=1.0, ell=100.0, d=0.2, delta=0.0, noise=0.3)
        d.update(kw); self.__dict__.update(d)
    @property
    def WL(self): return self.sigma_share*self.W

def piU(s,p):
    x=1-np.asarray(s,float)+p.eps; z=(x[0]-x[1]+p.delta)/p.noise
    return 1/(1+np.exp(-z))
def pis(s,p):
    pu=piU(s,p); return np.array([pu,1-pu])
def q(sj,p): return p.q0*(1-p.alpha*sj)
def Q(s,p):
    s=np.asarray(s,float); pi=pis(s,p)
    if p.tech=='ll': return float(pi[0]*q(s[0],p)+pi[1]*q(s[1],p))
    if p.tech=='avg': return p.q0*(1-p.alpha*s.mean())
    if p.tech=='wl': return p.q0*(1-p.alpha*s.mean()**2)
def u(sj,s,p,j,e):
    s=np.array(s,float); s[j]=sj; qq=Q(s,p); pi=pis(s,p)[j]
    return (1-qq)*pi*p.WL - p.kappa*sj**2/2 - p.lam*qq*p.ell - p.rho*e*(1-sj)
GR=np.linspace(0,1,201)
def br(s,p,j,e):
    vals=np.array([u(g,s,p,j,e) for g in GR]); k=int(np.argmax(vals)); a,b=GR[max(k-1,0)],GR[min(k+1,200)]
    for _ in range(35):
        c=b-(b-a)*0.618; dd=a+(b-a)*0.618
        if u(c,s,p,j,e)>u(dd,s,p,j,e): b=dd
        else: a=c
    return (a+b)/2
def lab_eq(p,e):
    s=np.array([0.3,0.3])
    for it in range(400):
        sn=np.array([br(s,p,j,e[j]) for j in range(2)])
        if np.max(np.abs(sn-s))<1e-6: return sn
        s=0.5*s+0.5*sn
    return s
def profile(p,r):
    e=np.array([p.m*r[0],p.m*r[1]]); s=lab_eq(p,e); qq=Q(s,p); pi=pis(s,p)
    V=[(1-qq)*pi[i]*p.W - p.Lam*qq*p.ell - p.d*r[i] for i in range(2)]
    return dict(V=V,s=s,Q=qq,piU=float(pi[0]))
def game(p):
    P={r:profile(p,r) for r in [(1,1),(1,0),(0,1),(0,0)]}  # (rU,rC)
    # U best response to C's action
    ne=[]
    for rU in (0,1):
        for rC in (0,1):
            okU = P[(rU,rC)]['V'][0] >= P[(1-rU,rC)]['V'][0]-1e-9
            okC = P[(rU,rC)]['V'][1] >= P[(rU,1-rC)]['V'][1]-1e-9
            if okU and okC: ne.append((rU,rC))
    # U's gains from regulating
    DU_R = P[(1,1)]['V'][0]-P[(0,1)]['V'][0]; DU_A = P[(1,0)]['V'][0]-P[(0,0)]['V'][0]
    DC_R = P[(1,1)]['V'][1]-P[(1,0)]['V'][1]; DC_A = P[(0,1)]['V'][1]-P[(0,0)]['V'][1]
    return P, ne, (DU_R,DU_A,DC_R,DC_A)
if __name__=="__main__":
    t0=time.time()
    for noise in [0.3,0.1]:
        print(f"\n=== leader-link risk, decisiveness noise={noise}, ell=60 ===")
        for delta in [0.0,0.1,0.2,0.3,0.5,0.8]:
            p=Par(delta=delta,noise=noise,ell=60.0); P,ne,D=game(p)
            lab=lambda r: f"s=({P[r]['s'][0]:.2f},{P[r]['s'][1]:.2f}) Q={P[r]['Q']:.3f} piU={P[r]['piU']:.2f}"
            print(f"delta={delta:.1f} NE={ne} DU(R,A)=({D[0]:+.2f},{D[1]:+.2f}) DC(R,A)=({D[2]:+.2f},{D[3]:+.2f}) | AA {lab((0,0))} | RA(U reg) {lab((1,0))} | RR {lab((1,1))}")
            sys.stdout.flush()
    print("time",round(time.time()-t0))
