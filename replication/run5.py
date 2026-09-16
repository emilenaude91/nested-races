from fast import *
pp=Par(kappa=20, ell=60, rho=16); n=2
e0 = -( pp.WL*(pp.q0*pp.alpha/n**2 - (1-pp.q0)*(n-1)/(n**2*(1+pp.eps))) + pp.lam*pp.q0*pp.alpha*pp.ell/n )/pp.rho
sc = min(1.0, pp.q0*pp.alpha*(pp.W/2 + pp.Lam*pp.ell)/pp.kappa); q = pp.q0*(1-pp.alpha*sc)
ebar = -( pp.WL*(pp.q0*pp.alpha/n**2 - (1-q)*(n-1)/(n**2*(1-sc+pp.eps))) - pp.kappa*sc + pp.lam*pp.q0*pp.alpha*pp.ell/n )/pp.rho
print("rho=16: e0=%.3f sc=%.3f ebar=%.3f s*(m=ebar)=%.3f" % (max(e0,0), sc, ebar, lab_eq(pp, np.array([ebar,ebar]))[0]))
# lambda-bar baseline
p=Par(); print("lambda-bar(0) n=2:", p.WL*((1-p.q0)/(1+p.eps) - p.q0*p.alpha)/(2*p.q0*p.alpha*p.ell))
# check u_RR vs u_AA at m=0.2, ell=30 (deadlock) and PD boundaries
for ell in [30,40,50]:
    r=interstate(Par(ell=ell, m=0.2)); print(ell, round(r['uRR'],3), round(r['uAA'],3), r['kind'])
