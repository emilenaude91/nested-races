import numpy as np
from scipy.optimize import brentq

class Par:
    def __init__(self, **kw):
        d = dict(W=10.0, sigma=0.5, kappa=6.0, q0=0.2, alpha=0.8, eps=0.05, gamma=1.0,
                 lam=0.05, Lam=1.0, rho=4.0, m=1.0, ell=100.0, d=0.2)
        d.update(kw); self.__dict__.update(d)
    @property
    def WL(self): return self.sigma * self.W

def Q(s, p):
    sb = np.mean(s); return p.q0 * (1 - p.alpha * sb ** p.gamma)
def dQ(s, p, n):
    sb = np.mean(s)
    return -p.q0 * p.alpha * p.gamma * (sb ** (p.gamma - 1) if p.gamma != 1 else 1.0) / n
def pi_all(s, p):
    v = 1 - np.asarray(s) + p.eps; return v / v.sum()

def du(sj, s, p, j, e):
    s = np.array(s, float); s[j] = sj; n = len(s)
    v = 1 - s + p.eps; V = v.sum(); pij = v[j] / V; dpij = -(V - v[j]) / V**2
    q = Q(s, p); dq = dQ(s, p, n)
    return -dq * pij * p.WL + (1 - q) * dpij * p.WL - p.kappa * sj - p.lam * dq * p.ell + p.rho * e

def br(s, p, j, e):
    if du(0.0, s, p, j, e) <= 0: return 0.0
    if du(1.0, s, p, j, e) >= 0: return 1.0
    return brentq(lambda x: du(x, s, p, j, e), 0.0, 1.0, xtol=1e-13)

def lab_eq(p, e, s0=None):
    n = len(e); s = np.full(n, 0.3) if s0 is None else np.array(s0, float)
    for it in range(5000):
        s_new = np.array([br(s, p, j, e[j]) for j in range(n)])
        if np.max(np.abs(s_new - s)) < 1e-12: return s_new
        s = 0.5 * s + 0.5 * s_new
    raise RuntimeError("no convergence")

def state_payoffs(p, r):
    e = np.array([p.m * r[0], p.m * r[1]]); s = lab_eq(p, e); q = Q(s, p); P = pi_all(s, p)
    V = [(1 - q) * P[i] * p.W - p.Lam * q * p.ell - p.d * r[i] for i in range(2)]
    return V, s, q

def interstate(p):
    (uRR, _), sRR, QRR = state_payoffs(p, (1, 1))
    (uAA, _), sAA, QAA = state_payoffs(p, (0, 0))
    (uAR, uRA), sAR, QAR = state_payoffs(p, (0, 1))
    DR = uRR - uAR; DA = uRA - uAA; Delta = 2*QAR - QRR - QAA
    if DR > 0 and DA > 0: kind = "Regulate dominant"
    elif DR < 0 and DA < 0: kind = "PD" if uRR > uAA else "Deadlock"
    elif DR > 0 > DA: kind = "Stag Hunt"
    else: kind = "Chicken"
    psi = (sAR[0]-sAA[0]) + (sAR[1]-sRR[0])
    return dict(uRR=uRR,uAA=uAA,uAR=uAR,uRA=uRA,DR=DR,DA=DA,Delta=Delta,kind=kind,
                sRR=sRR,sAA=sAA,sAR=sAR,QRR=QRR,QAA=QAA,QAR=QAR,psi=psi,
                idchk=(DR-DA)-(p.W/2+p.Lam*p.ell)*Delta)

def thresholds(base, key='ell', lo=1, hi=400):
    fA = lambda x: interstate(Par(**{**base, key: x}))['DA']
    fR = lambda x: interstate(Par(**{**base, key: x}))['DR']
    return brentq(fA, lo, hi), brentq(fR, lo, hi)
