from core4 import *
import time; t=time.time()
for tech in ['avg','wl','ll']:
    print(f"=== {tech} ===")
    for ell in [10,30,34,50,60,80,100,150,200]:
        c=cells(Par(tech=tech,ell=ell))
        print(f"ell={ell:4d} {c['kind']:8s} DR={c['DR']:+.3f} DA={c['DA']:+.3f} Delta={c['Delta']:+.4f} sAA={c['sAA'][0]:.3f} sRR={c['sRR'][0]:.3f} sAR=({c['sAR'][0]:.3f},{c['sAR'][1]:.3f}) Q=({c['QAA']:.3f},{c['QAR']:.3f},{c['QRR']:.3f}) id={c['idchk']:.1e} nloc={c['nloc']}")
    print("thresholds:",thr(dict(tech=tech)), "t=",round(time.time()-t))
