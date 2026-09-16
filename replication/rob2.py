from robust import *
print("--- falling-behind disutility eta (gamma=1) ---")
for eta in [0.0, 2.0, 5.0, 10.0]:
    ks=[]
    for ell in [10,30,45,60,90,150]:
        c=cells(Par(eta=eta,ell=ell)); ks.append((ell,c['kind'],round(c['Delta'],4),round(c['psi'],4),tuple(c['sAR'].round(3)),round(c['sRR'][0],3)))
    print("eta=",eta,ks); print("   thr(DA,DR)=",thr(dict(eta=eta)))
