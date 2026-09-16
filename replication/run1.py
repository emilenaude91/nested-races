from fast import *
import time; t=time.time()
print("baseline thresholds (ell_A, ell_R):", thresholds({}))
r = interstate(Par(ell=34)); print("ell=34:", r['kind'], round(r['DR'],4), round(r['DA'],4), "psi", round(r['psi'],4), "idchk", r['idchk'])
print("time", time.time()-t)
rng = np.random.default_rng(0); found=0; tot=0; kinds={}
for t_ in range(400):
    kw = dict(W=10, sigma=rng.uniform(0.2,1), kappa=rng.uniform(1,15), q0=rng.uniform(0.02,0.6),
              alpha=rng.uniform(0.3,1), eps=rng.uniform(0.01,0.3), lam=rng.uniform(0,0.3),
              rho=rng.uniform(0.5,10), m=rng.uniform(0.1,1), ell=rng.uniform(1,300), d=rng.uniform(0,1), Lam=rng.uniform(0.2,1))
    try: r = interstate(Par(**kw))
    except RuntimeError: continue
    tot+=1; kinds[r['kind']]=kinds.get(r['kind'],0)+1
    if r['psi'] < -1e-9 or r['Delta']>1e-9: found+=1; print("psi<0:", {k:round(v,3) for k,v in kw.items()}, round(r['psi'],5), round(r['Delta'],6), r['kind'])
print("random cases", tot, "psi<0 count", found, kinds, "time", time.time()-t)
