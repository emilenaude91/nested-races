import numpy as np, wmap
from wmap import Par, cells
# add contest exponent to wmap via monkeypatch
def pi_all_r(s,p):
    x=(1-np.asarray(s)+p.eps)**getattr(p,'r',1.0); return x/x.sum()
wmap.pi_all=pi_all_r
for r in [1.0,2.0,3.0]:
    for q0,W in [(0.19,13.0),(0.28,23.0),(0.19,7.0)]:
        p=Par(tech='ll',W=W,ell=100.0,q0=q0); p.r=r
        k,DR,DA,out=cells(p); pa=wmap.pi_all(out['AR'][1],p)[0]
        print(f"r={r} q0={q0} W/l={W/100:.2f}: form={['PD','Ch','SH','R'][k]} DR={DR:+.3f} DA={DA:+.3f} pi_A^AR={pa:.2f} sAR={out['AR'][1].round(3)}")
