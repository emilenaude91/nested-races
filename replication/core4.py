"""Unified solver (v4): techs 'avg' (linear mean), 'ir' (increasing returns, mean^gamma), 'wl' (soft-min weakest link),
'll' (leader-link with winner-causes-catastrophe accounting). Vectorised grid best responses + golden refine."""
import numpy as np
from scipy.optimize import brentq
class Par:
    def __init__(self, **kw):
        d=dict(W=10.0,sigma=0.5,kappa=6.0,q0=0.2,alpha=0.8,eps=0.05,gamma=2.0,tau=0.05,r=1.0,tech='avg',
               lam=0.05,Lam=1.0,rho=4.0,m=1.0,ell=100.0,d=0.2)
        d.update(kw); self.__dict__.update(d)
    @property
    def WL(self): return self.sigma*self.W
def pis(sU,sC,p):
    xU=(1-sU+p.eps)**p.r; xC=(1-sC+p.eps)**p.r; return xU/(xU+xC)
def qj(s,p): return p.q0*(1-p.alpha*s)
def Q(sU,sC,p):
    if p.tech=='avg': return p.q0*(1-p.alpha*(sU+sC)/2)
    if p.tech=='ir':  return p.q0*(1-p.alpha*((sU+sC)/2)**p.gamma)
    if p.tech=='wl':
        S=-p.tau*np.log((np.exp(-sU/p.tau)+np.exp(-sC/p.tau))/2); return p.q0*(1-p.alpha*S)
    if p.tech=='ll':
        pU=pis(sU,sC,p); return pU*qj(sU,p)+(1-pU)*qj(sC,p)
def prizeU(sU,sC,p):
    pU=pis(sU,sC,p)
    if p.tech=='ll': return pU*(1-qj(sU,p))
    return (1-Q(sU,sC,p))*pU
def u_vec(sgrid,s_other,p,j,e):
    """payoff of lab j over a grid of own safety, rival fixed."""
    if j==0: sU,sC=sgrid,s_other
    else: sU,sC=s_other,sgrid
    if p.tech=='ll':
        pj=pis(sU,sC,p) if j==0 else 1-pis(sU,sC,p); prize=pj*(1-qj(sgrid,p))
    else:
        pj=pis(sU,sC,p) if j==0 else 1-pis(sU,sC,p); prize=(1-Q(sU,sC,p))*pj
    return prize*p.WL - p.kappa*sgrid**2/2 - p.lam*Q(sU,sC,p)*p.ell - p.rho*e*(1-sgrid)
G=np.linspace(0,1,401)
def br(s_other,p,j,e):
    v=u_vec(G,s_other,p,j,e); k=int(np.argmax(v)); a,b=G[max(k-1,0)],G[min(k+1,400)]
    for _ in range(30):
        c=b-(b-a)*0.618; dd=a+(b-a)*0.618
        if u_vec(np.array([c]),s_other,p,j,e)[0]>u_vec(np.array([dd]),s_other,p,j,e)[0]: b=dd
        else: a=c
    nloc=int(np.sum((v[1:-1]>v[:-2])&(v[1:-1]>=v[2:])))
    return (a+b)/2, nloc
def lab_eq(p,eU,eC,s0=(0.3,0.3)):
    s=np.array(s0,float); nl=0
    for it in range(600):
        bU,n1=br(s[1],p,0,eU); bC,n2=br(s[0],p,1,eC); nl=max(nl,n1,n2)
        sn=np.array([bU,bC])
        if np.max(np.abs(sn-s))<1e-8: return sn,nl
        s=0.5*s+0.5*sn
    return s,nl
def cells(p):
    out={}; nl=0
    for name,(rU,rC) in [('RR',(1,1)),('AA',(0,0)),('AR',(0,1))]:
        s,n=lab_eq(p,p.m*rU,p.m*rC); nl=max(nl,n); qq=Q(s[0],s[1],p)
        VU=prizeU(s[0],s[1],p)*p.W - p.Lam*qq*p.ell - p.d*rU
        # column state (C): prize with roles swapped
        pC=1-pis(s[0],s[1],p); prC=(pC*(1-qj(s[1],p))) if p.tech=='ll' else (1-qq)*pC
        VC=prC*p.W - p.Lam*qq*p.ell - p.d*rC
        out[name]=dict(s=s,Q=qq,VU=VU,VC=VC)
    uRR=out['RR']['VU']; uAA=out['AA']['VU']; uAR=out['AR']['VU']; uRA=out['AR']['VC']
    DR=uRR-uAR; DA=uRA-uAA; Delta=2*out['AR']['Q']-out['RR']['Q']-out['AA']['Q']
    if DR>0 and DA>0: kind='Rdom'
    elif DR<0 and DA<0: kind='PD' if uRR>uAA else 'Deadlock'
    elif DR>0>DA: kind='SH'
    else: kind='Chicken'
    return dict(uRR=uRR,uAA=uAA,uAR=uAR,uRA=uRA,DR=DR,DA=DA,Delta=Delta,kind=kind,
                sAA=out['AA']['s'],sRR=out['RR']['s'],sAR=out['AR']['s'],QAA=out['AA']['Q'],QRR=out['RR']['Q'],QAR=out['AR']['Q'],
                idchk=(DR-DA)-(p.W/2+p.Lam*p.ell)*Delta,nloc=nl)
def thr(base,key='ell',lo=2,hi=400):
    f=lambda name: (lambda x: cells(Par(**{**base,key:x}))[name])
    res={}
    for name in ['DA','DR']:
        try: res[name]=brentq(f(name),lo,hi,xtol=0.05)
        except Exception: res[name]=None
    return res
