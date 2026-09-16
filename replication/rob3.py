from fast import *
print("eps sensitivity (thresholds ell_A, ell_R):")
for eps in [0.01,0.05,0.2]: print(" eps",eps,[round(x,1) for x in thresholds(dict(eps=eps))])
print("d sensitivity:")
for d in [0.0,0.2,1.0]:
    try: print(" d",d,[round(x,1) for x in thresholds(dict(d=d))])
    except Exception as ex: print(" d",d,"err",ex)
print("rho (speed asymmetry) -> Chicken band width:")
for rho in [2,4,8,16]:
    a,b=thresholds(dict(rho=rho)); print(" rho",rho,round(a,1),round(b,1),"width",round(b-a,2))
print("gamma=0.5 band:",[round(x,1) for x in thresholds(dict(gamma=0.5))])
# r=2 D_R threshold search
from robust import cells, Par as P2
from scipy.optimize import brentq
f=lambda x: cells(P2(r=2.0,ell=x))['DR']
print("r=2 DR at 60,100,150:",round(f(60),3),round(f(100),3),round(f(150),3)); print(" r=2 ell_R:",round(brentq(f,60,150,xtol=0.1),1))
