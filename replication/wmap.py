import numpy as np
from scipy.optimize import brentq
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

class Par:
    def __init__(self, **kw):
        d = dict(W=10.0, sigma=0.5, kappa=6.0, q0=0.2, alpha=0.8, eps=0.05, tech='avg',
                 lam=0.05, Lam=1.0, rho=4.0, m=1.0, ell=100.0, d=0.2)
        d.update(kw); self.__dict__.update(d)
    @property
    def WL(self): return self.sigma*self.W

def pi_all(s,p):
    x=1-np.asarray(s)+p.eps; return x/x.sum()
def Q(s,p):
    s=np.asarray(s,float)
    if p.tech=='avg': return p.q0*(1-p.alpha*s.mean())
    if p.tech=='wl':  return p.q0*(1-p.alpha*s.mean()**2)
    if p.tech=='ll':  return float(np.sum(pi_all(s,p)*p.q0*(1-p.alpha*s)))
def u(sj,s,p,j,e):
    s=np.array(s,float); s[j]=sj; q=Q(s,p); pi=pi_all(s,p)[j]
    return (1-q)*pi*p.WL - p.kappa*sj**2/2 - p.lam*q*p.ell - p.rho*e*(1-sj)
GR=np.linspace(0,1,201)
def br(s,p,j,e):
    vals=np.array([u(g,s,p,j,e) for g in GR]); k=int(np.argmax(vals))
    lo=GR[max(k-1,0)]; hi=GR[min(k+1,200)]
    # golden-section refine
    a,b=lo,hi
    for _ in range(40):
        c=b-(b-a)*0.618; d=a+(b-a)*0.618
        if u(c,s,p,j,e)>u(d,s,p,j,e): b=d
        else: a=c
    return (a+b)/2
def lab_eq(p,e):
    s=np.array([0.3,0.3])
    for it in range(400):
        sn=np.array([br(s,p,j,e[j]) for j in range(2)])
        if np.max(np.abs(sn-s))<1e-6: return sn
        s=0.5*s+0.5*sn
    return s
def cells(p):
    out={}
    for name,r in [('RR',(1,1)),('AA',(0,0)),('AR',(0,1))]:
        e=np.array([p.m*r[0],p.m*r[1]]); s=lab_eq(p,e); q=Q(s,p); pi=pi_all(s,p)
        V=[(1-q)*pi[i]*p.W - p.Lam*q*p.ell - p.d*r[i] for i in range(2)]; out[name]=(V,s,q)
    uRR=out['RR'][0][0]; uAA=out['AA'][0][0]; uAR,uRA=out['AR'][0]
    DR=uRR-uAR; DA=uRA-uAA
    if DR>0 and DA>0: k=3
    elif DR<0 and DA<0: k=0
    elif DR>0>DA: k=2
    else: k=1
    return k, DR, DA, out
if __name__=="__main__":
    import sys, json, time
    t0=time.time()
    ratios=np.geomspace(0.02,2.0,16); q0s=np.linspace(0.02,0.7,9)
    res={}
    for tech in ['avg','wl','ll']:
        M=np.zeros((len(q0s),len(ratios)),int)
        for i,q0 in enumerate(q0s):
            for j,r in enumerate(ratios):
                M[i,j]=cells(Par(tech=tech,W=100*r,ell=100.0,q0=q0))[0]
        res[tech]=M.tolist(); print(tech, "done", round(time.time()-t0)); sys.stdout.flush()
    json.dump(dict(ratios=ratios.tolist(),q0s=q0s.tolist(),res=res),open('wmap.json','w'))
