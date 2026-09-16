from fast import *
for g in [0.5, 1.0, 1.5, 2.0, 3.0]:
    print(f"\n=== gamma={g} ===")
    kinds=[]
    for ell in [5, 20, 30, 35, 40, 50, 60, 80, 100, 150, 200, 300, 500]:
        r = interstate(Par(ell=ell, gamma=g))
        kinds.append(r['kind'])
        print(f"ell={ell:4.0f} DR={r['DR']:+.3f} DA={r['DA']:+.3f} Delta={r['Delta']:+.4f} psi={r['psi']:+.4f} sAA={r['sAA'][0]:.3f} sRR={r['sRR'][0]:.3f} sAR={r['sAR'][0]:.3f},{r['sAR'][1]:.3f} QAA={r['QAA']:.3f} QAR={r['QAR']:.3f} QRR={r['QRR']:.3f} {r['kind']}")
    try:
        print("thresholds:", thresholds(dict(gamma=g), lo=1, hi=2000))
    except Exception as ex: print("thr err", ex)
