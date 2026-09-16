from robust import *
import sys
print("--- Tullock exponent r (gamma=1) ---")
for r in [0.5,1.0,1.5,2.0]:
    ks=[]; 
    for ell in [10,30,60,150]:
        c=cells(Par(r=r,ell=ell)); ks.append((ell,c['kind'],round(c['Delta'],4),round(c['psi'],4),c['maxlocal']))
    print("r=",r,ks, "thr(DA,DR)=",thr(dict(r=r)))
print("--- prize survives catastrophe (gamma=1) ---")
ks=[]
for ell in [10,30,60,150]:
    c=cells(Par(prize_survives=True,ell=ell)); ks.append((ell,c['kind'],round(c['Delta'],4),round(c['psi'],4)))
print(ks,"thr=",thr(dict(prize_survives=True)))
print("--- gamma=2 unimodality check ---")
for ell in [20,60,120,200]:
    c=cells(Par(gamma=2.0,ell=ell)); print(ell,c['kind'],"max local maxima on grid:",c['maxlocal'])
