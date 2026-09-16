from rfront import *
for g,ells in [(1.0,[28,32,35,37,50,60]),(2.0,[50,60,80,100,150])]:
    for ell in ells:
        p=Par(gamma=g,ell=ell,m=1.0); sym,U,Qm=sustainable(p)
        # also record best deviation from r=1
        j=10; dev=R[int(np.argmax(U[:,j]))]
        print(f"gamma={g} ell={ell} m=1: sustainable r={[round(float(x),1) for x in sym]} | best reply to rival r=1 is r={dev:.1f}")
