import numpy as np, json, sys
import lead3
from lead import Par, Q
D=json.load(open('lead_v2.json')); deltas=np.array(D['deltas'])
def cont_key(tech,noise,ell,r):
    k1=f"{tech}|{noise}|{ell:g}|{r[0]}{r[1]}"; k2=f"{tech}|{noise}|{r[0]}{r[1]}"
    return D['vals'][k1] if k1 in D['vals'] else D['vals'][k2]
def make_cont(noise,ell):
    def cont(tech,r,dl):
        rows=np.array(cont_key(tech,noise,ell,r)); dl=np.clip(dl,deltas[0],deltas[-1])
        return [np.interp(dl,rows[:,0],rows[:,k]) for k in (1,2,3,4,5,6,7,8)]
    return cont
def run(noise,ell,tau,theta,d0,tech='ll'):
    lead3.cont=make_cont(noise,ell)
    p=Par(delta=d0,noise=noise,ell=ell,tech=tech); p.tau=tau; p.theta=theta
    P,ne=lead3.game1(p); return P,ne
print("=== robustness: NE and RR/RA outcomes ===")
for (noise,ell,tau) in [(0.1,60.0,1.0),(0.3,60.0,1.0),(0.1,40.0,1.0),(0.1,100.0,1.0),(0.1,60.0,0.5)]:
    for theta in [0.0,0.2]:
        for d0 in [0.3,0.6]:
            P,ne=run(noise,ell,tau,theta,d0)
            ra=P[(1,0)]; rr=P[(1,1)]
            print(f"noise={noise} ell={ell:g} tau={tau} theta={theta} d0={d0}: NE={ne} | RA: sU1={ra['s1'][0]:.2f} d1={ra['d1']:+.2f} piU2={ra['piU2']:.2f} Q2={ra['Q2']:.3f} | RR: s1=({rr['s1'][0]:.2f},{rr['s1'][1]:.2f}) d1={rr['d1']:+.2f} sU2={rr['s2'][0]:.2f} piU2={rr['piU2']:.2f} Q2={rr['Q2']:.3f}")
            sys.stdout.flush()
print("\n=== who-paces map (ll, noise=0.1, ell=60, tau=1) ===")
grid={}
for theta in [0.0,0.1,0.2,0.3,0.4]:
    row=[]
    for d0 in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]:
        P,ne=run(0.1,60.0,1.0,theta,d0)
        # classify: lead lost under RA?
        lab={(1,1):'RR',(1,0):'RA',(0,1):'AR',(0,0):'AA'}
        tag='/'.join(lab[n] for n in ne) if ne else 'none'
        if 'RA' in tag and P[(1,0)]['piU2']<0.3: tag=tag.replace('RA','RA*')
        row.append(tag)
    grid[theta]=row; print(f"theta={theta}: "+" ".join(f"{t:>6}" for t in row)); sys.stdout.flush()
json.dump(grid,open('whopaces.json','w'))
