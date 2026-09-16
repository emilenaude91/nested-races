import numpy as np
from scipy.optimize import minimize_scalar, brentq

class Par:
    def __init__(self, **kw):
        d = dict(W=10.0, sigma=0.5, kappa=6.0, q0=0.2, alpha=0.8, eps=0.05, gamma=1.0, r=1.0,
                 lam=0.05, Lam=1.0, rho=4.0, m=1.0, ell=100.0, d=0.2, prize_survives=False, eta=0.0, tau=0.05)
        d.update(kw); self.__dict__.update(d)
    @property
    def WL(self): return self.sigma*self.W

def Q(s,p): return p.q0*(1-p.alpha*np.mean(s)**p.gamma)
def pis(s,p):
    x=(1-np.asarray(s)+p.eps)**p.r; return x/x.sum()
def u(sj,s,p,j,e):
    s=np.array(s,float); s[j]=sj; q=Q(s,p); pi=pis(s,p)[j]
    prize = pi*p.WL if p.prize_survives else (1-q)*pi*p.WL
    x=1-s+p.eps; behind = x[1-j]-x[j]   # positive when j is slower
    soft = p.tau*np.log1p(np.exp(behind/p.tau))  # softplus
    return prize - p.kappa*sj**2/2 - p.lam*q*p.ell - p.rho*e*(1-sj) - p.eta*soft
def br(s,p,j,e):
    grid=np.linspace(0,1,401); vals=np.array([u(g,s,p,j,e) for g in grid]); k=int(np.argmax(vals))
    lo,hi=grid[max(k-1,0)],grid[min(k+1,400)]
    res=minimize_scalar(lambda x:-u(x,s,p,j,e),bounds=(lo,hi),method='bounded',options=dict(xatol=1e-10))
    # unimodality check: count local maxima on grid
    lm = np.sum((vals[1:-1]>vals[:-2])&(vals[1:-1]>=vals[2:]))
    return res.x, lm
def lab_eq(p,e,s0=(0.3,0.3)):
    s=np.array(s0,float); mx=0
    for it in range(3000):
        out=[br(s,p,j,e[j]) for j in range(2)]; s_new=np.array([o[0] for o in out]); mx=max(mx,max(o[1] for o in out))
        if np.max(np.abs(s_new-s))<1e-9: return s_new,mx
        s=0.5*s+0.5*s_new
    return s,mx
def cells(p):
    out={}; mxs=0
    for name,r in [('RR',(1,1)),('AA',(0,0)),('AR',(0,1))]:
        e=np.array([p.m*r[0],p.m*r[1]]); s,mx=lab_eq(p,e); mxs=max(mxs,mx); q=Q(s,p); pi=pis(s,p)
        V=[(1-q)*pi[i]*p.W - p.Lam*q*p.ell - p.d*r[i] for i in range(2)]
        out[name]=(V,s,q)
    uRR=out['RR'][0][0]; uAA=out['AA'][0][0]; uAR,uRA=out['AR'][0]
    DR=uRR-uAR; DA=uRA-uAA; Delta=2*out['AR'][2]-out['RR'][2]-out['AA'][2]
    psi=(out['AR'][1][0]-out['AA'][1][0])+(out['AR'][1][1]-out['RR'][1][0])
    kind = "Rdom" if (DR>0 and DA>0) else ("PD" if (DR<0 and DA<0 and uRR>uAA) else ("Deadlock" if (DR<0 and DA<0) else ("StagHunt" if DR>0>DA else "Chicken")))
    return dict(DR=DR,DA=DA,Delta=Delta,psi=psi,kind=kind,sAA=out['AA'][1],sRR=out['RR'][1],sAR=out['AR'][1],maxlocal=mxs)
def thr(base, lo=1, hi=2000):
    fA=lambda x: cells(Par(**{**base,'ell':x}))['DA']; fR=lambda x: cells(Par(**{**base,'ell':x}))['DR']
    try: a=brentq(fA,lo,hi,xtol=0.05)
    except Exception: a=None
    try: b=brentq(fR,lo,hi,xtol=0.05)
    except Exception: b=None
    return a,b
