#!/usr/bin/env python3
"""Transferencia Xiamen <-> Caltech en ambos modos de procesamiento."""
import os, warnings
import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
warnings.filterwarnings("ignore")
SEED=42; PSEUDO=1e-6
BASE=os.path.expanduser("~/Documents/MicrobiomeRadar_Research")
TEA=f"{BASE}/01_TEA_Neurodivergent"
OUT=f"{BASE}/_crosscohort"
COH={"Xiamen":dict(dir=f"{TEA}/01g_Xiamen", nom="Xiamen", mapa="mapa_xiamen.tsv"),
     "Caltech":dict(dir=f"{TEA}/01h_Caltech_USA", nom="Caltech", mapa="mapa_caltech.tsv")}
def cargar(cf,modo):
    sub="02_processed" if modo=="single" else "02_processed_paired"
    X=pd.read_csv(f"{cf['dir']}/{sub}/{cf['nom']}_genus_counts.csv", index_col=0)
    mp=pd.read_csv(f"{cf['dir']}/00_screening/{cf['mapa']}", sep="\t", dtype=str)
    mp=mp[~mp["run"].duplicated(keep="first")].set_index("run")
    y,keep=[],[]
    for s in X.index:
        if s not in mp.index: continue
        g=mp.loc[s,"grupo"]
        if g=="ASD": y.append(1); keep.append(s)
        elif g=="Control": y.append(0); keep.append(s)
    return X.loc[keep], np.array(y)
def clr(df):
    rel=df.div(df.sum(axis=1),axis=0).replace(0,PSEUDO)
    lg=np.log(rel); return lg.sub(lg.mean(axis=1),axis=0)
def cv(X,y,nrep=5):
    o=[]
    for r in range(nrep):
        skf=StratifiedKFold(10,shuffle=True,random_state=SEED+r); pr=np.zeros(len(y))
        for tr,te in skf.split(X,y):
            m=RandomForestClassifier(n_estimators=500,max_features="sqrt",
                class_weight="balanced",random_state=SEED+r,n_jobs=-1)
            m.fit(X.iloc[tr],y[tr]); pr[te]=m.predict_proba(X.iloc[te])[:,1]
        o.append(roc_auc_score(y,pr))
    return float(np.mean(o))
print("="*76)
print("TRANSFERENCIA Xiamen <-> Caltech por modo de procesamiento")
print("="*76)
res=[]
for modo in ("single","paired"):
    d={n:cargar(cf,modo) for n,cf in COH.items()}
    com=sorted(set(d["Xiamen"][0].columns) & set(d["Caltech"][0].columns))
    feats=[g for g in com if min((d[n][0][g]>0).mean() for n in d)>=0.20]
    C={n:(clr(d[n][0][feats]),d[n][1]) for n in d}
    print(f"\n--- {modo.upper()} · {len(feats)} generos comunes ---")
    for a in C:
        Xa,ya=C[a]
        md=RandomForestClassifier(n_estimators=500,max_features="sqrt",
            class_weight="balanced",random_state=SEED,n_jobs=-1).fit(Xa,ya)
        for b in C:
            if a==b: continue
            Xb,yb=C[b]
            auc=roc_auc_score(yb, md.predict_proba(Xb[Xa.columns])[:,1])
            print(f"  {a:8} -> {b:8}  AUC = {auc:.3f}")
            res.append(dict(modo=modo, train=a, test=b, auc=round(auc,3)))
    # separabilidad de cohorte
    Xa,ya=C["Xiamen"]; Xb,yb=C["Caltech"]
    R=pd.concat([Xa,Xb],ignore_index=True); yy=np.array([0]*len(Xa)+[1]*len(Xb))
    sep=cv(R,yy)
    Ra=Xa[ya==0]; Rb=Xb[yb==0]
    Rc=pd.concat([Ra,Rb],ignore_index=True); yc=np.array([0]*len(Ra)+[1]*len(Rb))
    sepc=cv(Rc,yc)
    print(f"  separabilidad de cohorte: {sep:.3f}   solo controles: {sepc:.3f}")
    res.append(dict(modo=modo, train="separabilidad", test="todos", auc=round(sep,3)))
    res.append(dict(modo=modo, train="separabilidad", test="controles", auc=round(sepc,3)))
pd.DataFrame(res).to_csv(f"{OUT}/M3_transferencia_por_modo.csv", index=False)
print(f"\nGuardado en {OUT}/M3_transferencia_por_modo.csv")
print("="*76)
