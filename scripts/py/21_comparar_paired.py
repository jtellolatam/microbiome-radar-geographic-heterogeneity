#!/usr/bin/env python3
"""M3 - Compara single-end frente a paired-end en Xiamen y Caltech.
Cada cohorte se compara consigo misma: unica variable, el modo de procesamiento."""
import os, csv, warnings
import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
SEED=42; RNG=np.random.default_rng(42); PSEUDO=1e-6
BASE=os.path.expanduser("~/Documents/MicrobiomeRadar_Research")
TEA=os.path.join(BASE,"01_TEA_Neurodivergent")
OUT=os.path.join(BASE,"_crosscohort")
COH={
 "Xiamen": dict(dir=f"{TEA}/01g_Xiamen", mapa="mapa_xiamen.tsv", key="run"),
 "Caltech": dict(dir=f"{TEA}/01h_Caltech_USA", mapa="mapa_caltech.tsv", key="run"),
}
def cargar(cf, modo):
    sub = "02_processed" if modo=="single" else "02_processed_paired"
    nom = os.path.basename(cf["dir"]).split("_",1)[1].replace("_USA","")
    f = f"{cf['dir']}/{sub}/{nom}_genus_counts.csv"
    X = pd.read_csv(f, index_col=0)
    mp = pd.read_csv(f"{cf['dir']}/00_screening/{cf['mapa']}", sep="\t", dtype=str)
    mp = mp[~mp[cf["key"]].duplicated(keep="first")].set_index(cf["key"])
    y,keep=[],[]
    for s in X.index:
        base = s.split("_")[0] if s not in mp.index else s
        idx = s if s in mp.index else (base if base in mp.index else None)
        if idx is None: continue
        g = mp.loc[idx,"grupo"]
        if g=="ASD": y.append(1); keep.append(s)
        elif g=="Control": y.append(0); keep.append(s)
    return X.loc[keep], np.array(y)
def clr(df):
    rel=df.div(df.sum(axis=1),axis=0).replace(0,PSEUDO)
    lg=np.log(rel); return lg.sub(lg.mean(axis=1),axis=0)
def auc_cv(X,y,nrep=20):
    out=[]
    for r in range(nrep):
        skf=StratifiedKFold(10,shuffle=True,random_state=SEED+r); pr=np.zeros(len(y))
        for tr,te in skf.split(X,y):
            m=RandomForestClassifier(n_estimators=500,max_features="sqrt",
                class_weight="balanced",random_state=SEED+r,n_jobs=-1)
            m.fit(X.iloc[tr],y[tr]); pr[te]=m.predict_proba(X.iloc[te])[:,1]
        out.append(roc_auc_score(y,pr))
    return np.array(out), pr
def ic(yt,ys,n=1000):
    v=[];i=np.arange(len(yt))
    for _ in range(n):
        s=RNG.choice(i,len(i),True)
        if len(np.unique(yt[s]))<2: continue
        v.append(roc_auc_score(yt[s],ys[s]))
    return np.percentile(v,[2.5,97.5])
def perm(X,y,real,n=300):
    nul=[]
    for _ in range(n):
        yp=RNG.permutation(y)
        skf=StratifiedKFold(10,shuffle=True,random_state=SEED); pr=np.zeros(len(yp))
        for tr,te in skf.split(X,yp):
            m=RandomForestClassifier(n_estimators=100,max_features="sqrt",
                class_weight="balanced",random_state=SEED,n_jobs=-1)
            m.fit(X.iloc[tr],yp[tr]); pr[te]=m.predict_proba(X.iloc[te])[:,1]
        nul.append(roc_auc_score(yp,pr))
    return (np.sum(np.array(nul)>=real)+1)/(n+1)
print("="*84)
print("M3 · SINGLE-END frente a PAIRED-END")
print("="*84)
res=[]
for nom,cf in COH.items():
    print(f"\n--- {nom} ---")
    fila={"cohorte":nom}
    for modo in ("single","paired"):
        try:
            X,y = cargar(cf, modo)
        except Exception as e:
            print(f"  {modo}: no disponible ({e})"); continue
        feats=[c for c in X.columns if (X[c]>0).mean()>=0.20]
        Xc=clr(X[feats])
        a,pr = auc_cv(Xc,y)
        lo,hi = ic(y,pr); p = perm(Xc,y,a.mean())
        print(f"  {modo:7} n={len(y):3}  generos={len(feats):3}  "
              f"AUC={a.mean():.3f}  IC[{lo:.3f},{hi:.3f}]  p={p:.4f}  "
              f"profundidad={int(X.sum(axis=1).median()):,}")
        fila[f"auc_{modo}"]=round(a.mean(),3)
        fila[f"ic_{modo}"]=f"[{lo:.3f},{hi:.3f}]"
        fila[f"p_{modo}"]=round(p,4)
        fila[f"gen_{modo}"]=len(feats)
        fila[f"prof_{modo}"]=int(X.sum(axis=1).median())
    if "auc_single" in fila and "auc_paired" in fila:
        fila["dif"]=round(fila["auc_paired"]-fila["auc_single"],3)
        print(f"  {'diferencia':7} {fila['dif']:+.3f}")
    res.append(fila)
R=pd.DataFrame(res)
R.to_csv(f"{OUT}/M3_single_vs_paired.csv", index=False)
print("\n"+"="*84)
print(R.to_string(index=False))
print(f"\nGuardado en {OUT}/M3_single_vs_paired.csv")
