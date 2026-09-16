"""Appendix E: asymmetric enforcement. States differ in monitoring precision m_U, m_C."""
import json, sys
from core4 import Par, lab_eq, Q, pis, qj
def prof_payoffs(p, mU, mC, rU, rC):
    s,_=lab_eq(p, mU*rU, mC*rC); qq=Q(s[0],s[1],p); pU=pis(s[0],s[1],p); pC=1-pU
    if p.tech=='ll': prU=pU*(1-qj(s[0],p)); prC=pC*(1-qj(s[1],p))
    else: prU=(1-qq)*pU; prC=(1-qq)*pC
    return dict(VU=prU*p.W-p.Lam*qq*p.ell-p.d*rU, VC=prC*p.W-p.Lam*qq*p.ell-p.d*rC, s=s, Q=qq, pU=pU)
def game(p,mU,mC):
    P={k:prof_payoffs(p,mU,mC,*v) for k,v in dict(RR=(1,1),RA=(1,0),AR=(0,1),AA=(0,0)).items()}
    # U's gains from regulating: vs regulator (RR-AR), vs racer (RA-AA); C's: vs regulator (RR-RA), vs racer (AR-AA)
    DU_R=P['RR']['VU']-P['AR']['VU']; DU_A=P['RA']['VU']-P['AA']['VU']
    DC_R=P['RR']['VC']-P['RA']['VC']; DC_A=P['AR']['VC']-P['AA']['VC']
    def disp(DR,DA):
        if DR>=0 and DA>=0: return 'restrains regardless'
        if DR<0 and DA<0: return 'races regardless'
        if DR>=0 and DA<0: return 'follows'
        return 'compensates'
    ne=[]
    if DU_R>=0 and DC_R>=0: ne.append('RR')
    if DU_A<=0 and DC_A<=0: ne.append('AA')
    if DU_A>=0 and DC_R<=0: ne.append('RA')
    if DU_R<=0 and DC_A>=0: ne.append('AR')
    return dict(P=P,DU=(DU_R,DU_A),DC=(DC_R,DC_A),dispU=disp(DU_R,DU_A),dispC=disp(DC_R,DC_A),ne=ne)
rows=[]
for tech in ['avg','wl','ll']:
    for ell in [30,60,120]:
        for mU,mC in [(1,1),(1,0.5),(1,0.1),(1,0.0)]:
            g=game(Par(tech=tech,ell=ell),mU,mC)
            rows.append(dict(tech=tech,ell=ell,mU=mU,mC=mC,DU=[round(x,2) for x in g['DU']],DC=[round(x,2) for x in g['DC']],dispU=g['dispU'],dispC=g['dispC'],ne=g['ne'],
                             sRR=[round(x,2) for x in g['P']['RR']['s']],sRA=[round(x,2) for x in g['P']['RA']['s']],QRR=round(g['P']['RR']['Q'],3),QRA=round(g['P']['RA']['Q'],3),QAA=round(g['P']['AA']['Q'],3)))
            r=rows[-1]; print(f"{tech} ell={ell:3d} mC={mC:<4} DU={r['DU']} DC={r['DC']} U:{r['dispU']:<20} C:{r['dispC']:<20} NE={r['ne']} sRR={r['sRR']} sRA={r['sRA']} Q: RR {r['QRR']} RA {r['QRA']} AA {r['QAA']}"); sys.stdout.flush()
json.dump(rows,open('asym.json','w'),indent=1)
