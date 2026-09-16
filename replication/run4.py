from fast import *
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
# threshold ell_R and ell_A as functions of m (gamma=1 and gamma=2)
print("m, ell_A, ell_R (gamma=1) | ell_R, ell_A (gamma=2)")
for m in [0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0]:
    out=[m]
    for g in [1.0, 2.0]:
        try:
            fA = lambda x: interstate(Par(ell=x, m=m, gamma=g))['DA']
            fR = lambda x: interstate(Par(ell=x, m=m, gamma=g))['DR']
            out += [round(brentq(fA, 1, 3000),1), round(brentq(fR, 1, 3000),1)]
        except Exception as ex: out += ["none","none"]
    print(out)
# implementation illustration with feasible benchmark
def impl(pp, n=2):
    e0 = -( pp.WL*(pp.q0*pp.alpha/n**2 - (1-pp.q0)*(n-1)/(n**2*(1+pp.eps))) + pp.lam*pp.q0*pp.alpha*pp.ell/n )/pp.rho
    sc = min(1.0, pp.q0*pp.alpha*(pp.W/2 + pp.Lam*pp.ell)/pp.kappa); q = pp.q0*(1-pp.alpha*sc)
    ebar = -( pp.WL*(pp.q0*pp.alpha/n**2 - (1-q)*(n-1)/(n**2*(1-sc+pp.eps))) - pp.kappa*sc + pp.lam*pp.q0*pp.alpha*pp.ell/n )/pp.rho
    return max(e0,0), sc, ebar
for kw in [dict(kappa=20, ell=60), dict(kappa=20, ell=60, rho=8), dict(kappa=20, ell=60, lam=0.3), dict(kappa=20, ell=60, Lam=0.5), dict(kappa=20, ell=60, sigma=1.0)]:
    pp=Par(**kw); e0, sc, eb = impl(pp); s_at_1 = lab_eq(pp, np.array([1.0,1.0]))[0]
    print(kw, "e0=%.3f sc=%.3f ebar=%.3f s*(m=1)=%.3f %s" % (e0, sc, eb, s_at_1, "feasible" if eb<=1 else "INFEASIBLE"))
# Figure 1: D_R, D_A vs ell for gamma=1 and gamma=2
fig, axes = plt.subplots(1,2, figsize=(9,3.4))
for ax, g, ttl in zip(axes, [1.0, 2.0], [r"linear risk ($\gamma=1$)", r"breadth-of-compliance risk ($\gamma=2$)"]):
    ells = np.linspace(2, 200, 120); DR=[]; DA=[]
    for l in ells:
        r = interstate(Par(ell=l, gamma=g)); DR.append(r['DR']); DA.append(r['DA'])
    ax.plot(ells, DR, 'k-', label=r"$D_R=u_{RR}-u_{AR}$"); ax.plot(ells, DA, 'k--', label=r"$D_A=u_{RA}-u_{AA}$")
    ax.axhline(0, color='grey', lw=0.7); ax.set_xlabel(r"shared catastrophe loss $\ell$"); ax.set_title(ttl, fontsize=10)
    ax.set_ylim(-2, 6); ax.legend(fontsize=8, loc="upper left")
axes[0].set_ylabel("gain from Regulate")
plt.tight_layout(); plt.savefig("fig_taxonomy.pdf")
# Figure 2: thresholds vs m (gamma=1): ell_R(m)
ms = np.linspace(0.06, 1.0, 30); lR=[]; lA=[]
for m in ms:
    fR = lambda x: interstate(Par(ell=x, m=m))['DR']; fA = lambda x: interstate(Par(ell=x, m=m))['DA']
    lR.append(brentq(fR, 1, 5000)); lA.append(brentq(fA, 1, 5000))
fig, ax = plt.subplots(figsize=(4.6,3.2))
ax.plot(ms, lR, 'k-', label=r"$\ell_R(m)$: Regulate dominant above"); ax.plot(ms, lA, 'k--', label=r"$\ell_A(m)$: PD below")
ax.set_xlabel(r"monitoring precision $m$"); ax.set_ylabel(r"threshold loss $\ell$"); ax.set_yscale("log"); ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig("fig_thresholds_m.pdf"); print("figs done")
