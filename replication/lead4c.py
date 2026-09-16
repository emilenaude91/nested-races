"""Head-start model v4: difference-form contest with lead, leader-link risk with winner-causes-catastrophe accounting.
Two periods, absorbing catastrophe, no discounting, regime commitment, inertia-based equilibrium selection."""
import numpy as np, json, sys, time
class P:
    def __init__(self,**kw):
        d=dict(W=10.0,sigma=0.5,kappa=6.0,q0=0.2,alpha=0.8,eps=0.05,lam=0.05,Lam=1.0,rho=4.0,m=1.0,ell=60.0,d=0.2,
               delta=0.3,noise=0.1,tau=1.0,theta=0.0); d.update(kw); self.__dict__.update(d)
    @property
    def WL(self): return self.sigma*self.W
def piU(sU,sC,p,delta): return 1/(1+np.exp(-((1-sU)-(1-sC)+delta)/p.noise))
def q(s,p): return p.q0*(1-p.alpha*s)
def Qll(sU,sC,p,delta):
    pu=piU(sU,sC,p,delta); return pu*q(sU,p)+(1-pu)*q(sC,p)
G=np.linspace(0,1,201)
def u2_vec(sg,so,p,j,e,delta):   # period-2 lab payoff over own grid
    sU,sC=(sg,so) if j==0 else (so,sg); pu=piU(sU,sC,p,delta); pj=pu if j==0 else 1-pu
    return pj*(1-q(sg,p))*p.WL - p.kappa*sg**2/2 - p.lam*Qll(sU,sC,p,delta)*p.ell - p.rho*e*(1-sg)
def eq2(p,eU,eC,delta,s0=(0.3,0.3)):
    s=np.array(s0,float)
    for it in range(500):
        bU=G[np.argmax(u2_vec(G,s[1],p,0,eU,delta))]; bC=G[np.argmax(u2_vec(G,s[0],p,1,eC,delta))]
        sn=np.array([bU,bC])
        if np.max(np.abs(sn-s))<1e-9: return sn
        s=0.5*s+0.5*sn
    return s
def vals2(p,reg,delta):
    s=eq2(p,p.m*reg[0],p.m*reg[1],delta); pu=piU(s[0],s[1],p,delta); qq=Qll(s[0],s[1],p,delta)
    VU=pu*(1-q(s[0],p))*p.W-p.Lam*qq*p.ell-p.d*reg[0]; VC=(1-pu)*(1-q(s[1],p))*p.W-p.Lam*qq*p.ell-p.d*reg[1]
    LU=u2_vec(np.array([s[0]]),s[1],p,0,p.m*reg[0],delta)[0]; LC=u2_vec(np.array([s[1]]),s[0],p,1,p.m*reg[1],delta)[0]
    return np.array([VU,VC,LU,LC,qq,pu,s[0],s[1]])
DG=np.round(np.arange(-0.8,1.61,0.1),2)
def build_tables(p):
    T={}
    for reg in [(1,1),(1,0),(0,1),(0,0)]: T[reg]=np.array([vals2(p,reg,float(dl)) for dl in DG])
    return T
def cont(T,reg,dl):
    dl=np.clip(dl,DG[0],DG[-1]); return np.array([np.interp(dl,DG,T[reg][:,k]) for k in range(8)])
def u1_vec(sg,so,p,j,e,reg,T):
    sU,sC=(sg,so) if j==0 else (so,sg); q1=Qll(sU,sC,p,p.delta); d1=p.delta+p.tau*((1-sU)-(1-sC))-p.theta
    d1c=np.clip(d1,DG[0],DG[-1]); Lj=np.interp(d1c,DG,T[reg][:,2+j])
    return -p.kappa*sg**2/2 - p.lam*q1*p.ell - p.rho*e*(1-sg) + (1-q1)*Lj
G1=np.linspace(0,1,201)
def eq1(p,reg,T,s0):
    s=np.array(s0,float)
    for it in range(500):
        bU=G1[np.argmax(u1_vec(G1,s[1],p,0,p.m*reg[0],reg,T))]; bC=G1[np.argmax(u1_vec(G1,s[0],p,1,p.m*reg[1],reg,T))]
        sn=np.array([bU,bC])
        if np.max(np.abs(sn-s))<1e-9: break
        s=0.5*s+0.5*sn
    q1=Qll(s[0],s[1],p,p.delta); d1=p.delta+p.tau*((1-s[0])-(1-s[1]))-p.theta; c=cont(T,reg,d1)
    VU=-p.Lam*q1*p.ell-p.d*reg[0]+(1-q1)*c[0]; VC=-p.Lam*q1*p.ell-p.d*reg[1]+(1-q1)*c[1]
    return dict(s1=s,Q1=q1,d1=d1,VU=VU,VC=VC,Q2=c[4],piU2=c[5],s2=(c[6],c[7]))
def game1(p,T,prev):
    Pm={}; 
    for reg in [(1,1),(1,0),(0,1),(0,0)]: Pm[reg]=eq1(p,reg,T,prev.get(reg,(0.3,0.3)))
    ne=[]
    for rU in (0,1):
        for rC in (0,1):
            if Pm[(rU,rC)]['VU']>=Pm[(1-rU,rC)]['VU']-1e-9 and Pm[(rU,rC)]['VC']>=Pm[(rU,1-rC)]['VC']-1e-9: ne.append((rU,rC))
    return Pm,ne
if __name__=="__main__":
    t0=time.time(); results={}
    for noise,ell in [(0.1,60.0),(0.3,60.0),(0.1,100.0)]:
        p=P(noise=noise,ell=ell); T=build_tables(p); print("tables",noise,ell,round(time.time()-t0)); sys.stdout.flush()
        for theta in [0.0,0.2,0.4]:
            prev={}; row=[]
            for d0 in np.arange(0.1,0.91,0.1):
                p=P(noise=noise,ell=ell,theta=theta,delta=float(d0)); Pm,ne=game1(p,T,prev); prev={k:v['s1'] for k,v in Pm.items()}
                lab={(1,1):'RR',(1,0):'RA',(0,1):'AR',(0,0):'AA'}; tag='/'.join(lab[n] for n in ne) or 'none'
                ra=Pm[(1,0)]; rr=Pm[(1,1)]
                row.append(dict(d0=round(float(d0),1),ne=tag,RA=dict(sU1=round(float(ra['s1'][0]),2),d1=round(float(ra['d1']),2),piU2=round(float(ra['piU2']),2),Q2=round(float(ra['Q2']),3)),
                                RR=dict(s1=[round(float(x),2) for x in rr['s1']],d1=round(float(rr['d1']),2),sU2=round(float(rr['s2'][0]),2),piU2=round(float(rr['piU2']),2),Q2=round(float(rr['Q2']),3))))
                print(f"noise={noise} ell={ell:g} theta={theta} d0={d0:.1f} NE={tag:6s} | RA sU1={ra['s1'][0]:.2f} d1={ra['d1']:+.2f} piU2={ra['piU2']:.2f} Q2={ra['Q2']:.3f} | RR s1=({rr['s1'][0]:.2f},{rr['s1'][1]:.2f}) d1={rr['d1']:+.2f} sU2={rr['s2'][0]:.2f} piU2={rr['piU2']:.2f} Q2={rr['Q2']:.3f}")
                sys.stdout.flush()
            results[f"{noise}|{ell:g}|{theta}"]=row
    json.dump(results,open('lead_v4.json','w')); print("t",round(time.time()-t0))
