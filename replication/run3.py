from fast import *
import numpy as np
p = Par()
# --- verify cross partial formula (gamma=1, n=2)
def u(sj, sk, e=0.0):
    s = np.array([sj, sk]); v = 1 - s + p.eps; V = v.sum(); q = Q(s, p)
    return (1-q)*(v[0]/V)*p.WL - p.kappa*sj**2/2 - p.lam*q*p.ell - p.rho*e*(1-sj)
h=1e-5
for (sj, sk) in [(0.2,0.5),(0.5,0.2),(0.3,0.3)]:
    num = (u(sj+h,sk+h)-u(sj+h,sk-h)-u(sj-h,sk+h)+u(sj-h,sk-h))/(4*h*h)
    s=np.array([sj,sk]); v=1-s+p.eps; V=v.sum(); q=Q(s,p); B=v[1]
    ana = p.WL*(v[0]-B)*(p.q0*p.alpha/(2*V**2) + (1-q)/V**3)
    print("cross partial num/ana", round(num,6), round(ana,6))
# --- verify symmetric FOC formula for n labs
def F(s, n, e=0.0):
    q = p.q0*(1-p.alpha*s)
    return p.WL*(p.q0*p.alpha/n**2 - (1-q)*(n-1)/(n**2*(1-s+p.eps))) - p.kappa*s + p.lam*p.q0*p.alpha*p.ell/n + p.rho*e
for n in [2,3,5]:
    s=np.full(n,0.3); print("FOC num/ana n=",n, round(du(0.3,s,p,0,0.0),6), round(F(0.3,n),6))
# --- n-lab comparative static at symmetric eq (e=0 and e=m), various lam
def sym_eq(n, e, lam):
    pp = Par(lam=lam)
    f = lambda s: (lambda q: pp.WL*(pp.q0*pp.alpha/n**2 - (1-q)*(n-1)/(n**2*(1-s+pp.eps))) - pp.kappa*s + pp.lam*pp.q0*pp.alpha*pp.ell/n + pp.rho*e)(pp.q0*(1-pp.alpha*s))
    if f(0) <= 0: return 0.0
    if f(1) >= 0: return 1.0
    return brentq(f, 0, 1)
print("\n n-lab symmetric equilibrium safety (e=m=1), lam=0.05 / 0.3 / 1.0; and e=0 lam=1.0")
for n in [2,3,4,6,8,12]:
    print(n, round(sym_eq(n,1.0,0.05),4), round(sym_eq(n,1.0,0.3),4), round(sym_eq(n,1.0,1.0),4), "| e=0:", round(sym_eq(n,0.0,1.0),4), round(sym_eq(n,0.0,3.0),4))
# --- implementation thresholds under (R,R), n=2: bite threshold e0 and target threshold ebar
def thresholds_impl(pp, n=2):
    q00 = pp.q0
    e0 = -( pp.WL*(pp.q0*pp.alpha/n**2 - (1-q00)*(n-1)/(n**2*(1+pp.eps))) + pp.lam*pp.q0*pp.alpha*pp.ell/n )/pp.rho
    # benchmark: symmetric s maximising (1-Q)W/2 - Lam Q ell - kappa s^2/2  (one lab per state, common s, dQ/ds = -q0 alpha)
    sc = min(1.0, pp.q0*pp.alpha*(pp.W/2 + pp.Lam*pp.ell)/pp.kappa)
    q = pp.q0*(1-pp.alpha*sc)
    ebar = -( pp.WL*(pp.q0*pp.alpha/n**2 - (1-q)*(n-1)/(n**2*(1-sc+pp.eps))) - pp.kappa*sc + pp.lam*pp.q0*pp.alpha*pp.ell/n )/pp.rho
    return max(e0,0), sc, ebar
for kw in [dict(), dict(lam=0.3), dict(rho=8.0), dict(Lam=0.5), dict(ell=200), dict(kappa=3.0)]:
    pp=Par(**kw); e0, sc, eb = thresholds_impl(pp)
    print(kw, "bite e0=", round(e0,3), "s^c=", round(sc,3), "ebar=", round(eb,3), "feasible" if eb<=1 else "INFEASIBLE (>1)")
# --- risk dominance band for gamma=2
g=2.0
fR = lambda x: interstate(Par(ell=x,gamma=g))['DR']
fRD = lambda x: (lambda r: r['DR'] + r['DA'])(interstate(Par(ell=x,gamma=g)))
lR = brentq(fR, 10, 400); lrd = brentq(fRD, lR, 400); 
fA = lambda x: interstate(Par(ell=x,gamma=g))['DA']; lA = brentq(fA, lR, 800)
print("\ngamma=2: ell_R=%.2f  ell_rd (A,A risk-dom below)=%.2f  ell_A=%.2f" % (lR, lrd, lA))
r = interstate(Par(ell=60, gamma=g)); print("ell=60 gamma=2 cells:", {k:round(r[k],3) for k in ['uRR','uAA','uAR','uRA','DR','DA']}, "pstar=", round((-r['DA'])/((-r['DA'])+r['DR']),3))
r = interstate(Par(ell=120, gamma=g)); print("ell=120 gamma=2 cells:", {k:round(r[k],3) for k in ['uRR','uAA','uAR','uRA','DR','DA']})
# baseline table for the paper
print("\nBaseline gamma=1 table")
for ell in [10, 30, 33.5, 60, 150]:
    r = interstate(Par(ell=ell)); print(ell, {k:round(r[k],3) for k in ['uRR','uAA','uAR','uRA','DR','DA','Delta','psi','QAA','QAR','QRR']}, r['sAA'].round(3), r['sAR'].round(3), r['sRR'].round(3), r['kind'])
print("\nweak monitoring m=0.2 baseline")
for ell in [30, 60, 150]:
    r = interstate(Par(ell=ell, m=0.2)); print(ell, {k:round(r[k],3) for k in ['uRR','uAA','uAR','uRA','DR','DA']}, r['sRR'].round(3), r['kind'])
