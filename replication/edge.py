"""Asymmetric resources: U's lab has a speed multiplier a (compute/talent) -> x_U = a(1-s_U+eps). Grid best-response solver."""
import numpy as np
P=dict(W=10,sigma=0.5,kappa=6,q0=0.2,alpha=0.8,eps=0.05,lam=0.05,Lam=1,rho=4,d=0.2,tau=0.05)
G=np.linspace(0,1,401)
def speeds(sU,sC,a): return a*(1-sU+P['eps']), (1-sC+P['eps'])
def pi(sU,sC,a): xU,xC=speeds(sU,sC,a); return xU/(xU+xC)
def q(s): return P['q0']*(1-P['alpha']*s)
def Q(sU,sC,a,tech):
    if tech=='avg': return P['q0']*(1-P['alpha']*(sU+sC)/2)
    if tech=='wl':
        t=P['tau']; S=-t*np.log(0.5*(np.exp(-sU/t)+np.exp(-sC/t))); return P['q0']*(1-P['alpha']*S)
    pu=pi(sU,sC,a); return pu*q(sU)+(1-pu)*q(sC)
def u(j,so,sr,e,a,tech,ell):
    sU,sC=(so,sr) if j==0 else (sr,so); pu=pi(sU,sC,a); pj=pu if j==0 else 1-pu; qq=Q(sU,sC,a,tech)
    prize=pj*(1-q(so)) if tech=='ll' else (1-qq)*pj
    return prize*P['sigma']*P['W']-P['kappa']*so**2/2-P['lam']*qq*ell-P['rho']*e*(1-so)
def br(j,sr,e,a,tech,ell): vals=[u(j,g,sr,e,a,tech,ell) for g in G]; return G[int(np.argmax(vals))]
def eq(eU,eC,a,tech,ell):
    s=[0.3,0.3]
    for _ in range(200):
        n=[br(0,s[1],eU,a,tech,ell),br(1,s[0],eC,a,tech,ell)]
        if abs(n[0]-s[0])<1e-9 and abs(n[1]-s[1])<1e-9: return n
        s=[0.5*s[0]+0.5*n[0],0.5*s[1]+0.5*n[1]]
    return s
def cells(a,tech,ell,m=1.0):
    out={}
    for name,(rU,rC) in dict(RR=(1,1),RA=(1,0),AR=(0,1),AA=(0,0)).items():
        s=eq(m*rU,m*rC,a,tech,ell); pu=pi(s[0],s[1],a); qq=Q(s[0],s[1],a,tech)
        prU=pu*(1-q(s[0])) if tech=='ll' else (1-qq)*pu; prC=(1-pu)*(1-q(s[1])) if tech=='ll' else (1-qq)*(1-pu)
        out[name]=dict(VU=prU*P['W']-P['Lam']*qq*ell-P['d']*rU, VC=prC*P['W']-P['Lam']*qq*ell-P['d']*rC, s=s, Q=qq, pu=pu)
    DU=(out['RR']['VU']-out['AR']['VU'], out['RA']['VU']-out['AA']['VU']); DC=(out['RR']['VC']-out['RA']['VC'], out['AR']['VC']-out['AA']['VC'])
    def disp(DR,DA): return 'restrains' if DR>=0 and DA>=0 else 'races' if DR<0 and DA<0 else 'follows' if DR>=0 else 'compensates'
    ne=[]
    if DU[0]>=0 and DC[0]>=0: ne.append('RR')
    if DU[1]<=0 and DC[1]<=0: ne.append('AA')
    if DU[1]>=0 and DC[0]<=0: ne.append('RA')
    if DU[0]<=0 and DC[1]>=0: ne.append('AR')
    return out,DU,DC,disp(*DU),disp(*DC),ne
for tech in ['avg','wl','ll']:
    for ell in [30,60,120]:
        for a in [1.0,1.5,2.5]:
            out,DU,DC,dU,dC,ne=cells(a,tech,ell)
            print(f"{tech} ell={ell:3d} a={a}: piU(AA)={out['AA']['pu']:.2f} sAA=({out['AA']['s'][0]:.2f},{out['AA']['s'][1]:.2f}) sRR=({out['RR']['s'][0]:.2f},{out['RR']['s'][1]:.2f}) DU=({DU[0]:+.2f},{DU[1]:+.2f}) DC=({DC[0]:+.2f},{DC[1]:+.2f}) U:{dU:<9} C:{dC:<9} NE={ne}",flush=True)
