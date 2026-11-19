#!/usr/bin/env python3
"""Microbiome Radar - Geografia frente a condicion"""
import os, warnings
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
SEED=42; PSEUDO=1e-6; NREP=5
BASE=os.path.expanduser("~/Documents/MicrobiomeRadar_Research")
AD=os.path.join(BASE,"03_Alzheimer_Neurodegenerative")
TEA=os.path.join(BASE,"01_TEA_Neurodivergent")
PD=os.path.join(BASE,"02_Parkinson_Neurodegenerative")
OUT=os.path.join(BASE,"_crosscohort")
BRAZOS = {
 "Alzheimer": {
   "Shanghai": dict(m=f"{AD}/03d_Li2019_Shanghai_China/02_processed/Li_genus_counts.csv",
     p=f"{AD}/03d_Li2019_Shanghai_China/00_screening/mapa_heces.tsv",
     key="run_accession", col="grupo", modo="exacto", caso=["AD"], ctrl=["NC"]),
   "Lishui": dict(m=f"{AD}/03g_Ling2021_Lishui_China/02_processed/Ling_genus_counts.csv",
     p=f"{AD}/03g_Ling2021_Lishui_China/00_screening/ena_PRJNA633959.tsv",
     key="run_accession", col="library_name", modo="prefijo"),
   "Chongqing": dict(m=f"{AD}/03h_Zhuang2019_Chongqing_China/02_processed/Zhuang_genus_counts.csv",
     p=f"{AD}/03h_Zhuang2019_Chongqing_China/00_screening/meta_554111.tsv",
     key="run_accession", col="sample_alias", modo="prefijo"),
   "Estambul": dict(m=f"{AD}/03i_Yildirim2022_Turkey/02_processed/Turquia_genus_counts.csv",
     p=f"{AD}/03i_Yildirim2022_Turkey/00_screening/ena_PRJNA734525.tsv",
     key="run_accession", col="library_name", modo="turquia", centro="Estambul"),
   "Kayseri": dict(m=f"{AD}/03i_Yildirim2022_Turkey/02_processed/Turquia_genus_counts.csv",
     p=f"{AD}/03i_Yildirim2022_Turkey/00_screening/ena_PRJNA734525.tsv",
     key="run_accession", col="library_name", modo="turquia", centro="Kayseri"),
   "Brescia": dict(m=f"{AD}/03k_Italia_PRJEB89799/02_processed/Italia_genus_counts.csv",
     p=f"{AD}/03k_Italia_PRJEB89799/00_screening/mapa_baseline.tsv",
     key="run", col="grupo", modo="exacto", caso=["AD"], ctrl=["Control"]),
 },
 "TEA": {
   "Nanjing": dict(m=f"{TEA}/01f_Dan_Nanjing/02_processed/Dan_genus_counts.csv",
     p=f"{TEA}/01f_Dan_Nanjing/00_screening/mapa_dan.tsv",
     key="run", col="grupo", modo="exacto", caso=["ASD"], ctrl=["Control"]),
   "Xiamen": dict(m=f"{TEA}/01g_Xiamen/02_processed/Xiamen_genus_counts.csv",
     p=f"{TEA}/01g_Xiamen/00_screening/mapa_xiamen.tsv",
     key="run", col="grupo", modo="exacto", caso=["ASD"], ctrl=["Control"]),
   "Caltech": dict(m=f"{TEA}/01h_Caltech_USA/02_processed/Caltech_genus_counts.csv",
     p=f"{TEA}/01h_Caltech_USA/00_screening/mapa_caltech.tsv",
     key="run", col="grupo", modo="exacto", caso=["ASD"], ctrl=["Control"]),
 },
 "Parkinson": {
   "Finlandia": dict(m=f"{PD}/02e_Finlandia/02_processed/Finlandia_genus_counts.csv",
     p=f"{PD}/02e_Finlandia/00_screening/mapa_finlandia.tsv",
     key="run", col="grupo", modo="exacto", caso=["PD"], ctrl=["Control"]),
   "Rumania": dict(m=f"{PD}/02f_Rumania/02_processed/Rumania_genus_counts.csv",
     p=f"{PD}/02f_Rumania/00_screening/mapa_rumania.tsv",
     key="muestra", col="grupo", modo="exacto", caso=["PD"], ctrl=["Control"]),
 },
}
def log(*a): print(*a, flush=True)
def etiqueta(cf, v):
    s=str(v).upper()
    if cf["modo"]=="exacto":
        if v in cf["caso"]: return 1
        if v in cf["ctrl"]: return 0
        return None
    if cf["modo"]=="turquia":
        pref="".join(c for c in str(v) if not c.isdigit())
        if cf["centro"]=="Kayseri":
            if pref=="ADH": return 1
            if pref=="CK": return 0
        else:
            if pref=="AD": return 1
            if pref=="C": return 0
        return None
    if s.startswith("CON") or s.startswith("C"): return 0
    if s.startswith("AD"): return 1
    return None
def cargar(cf):
    X=pd.read_csv(cf["m"], index_col=0)
    mp=pd.read_csv(cf["p"], sep="\t", dtype=str)
    mp=mp[~mp[cf["key"]].duplicated(keep="first")].set_index(cf["key"])
    y,keep=[],[]
    for s in X.index:
        if s not in mp.index: continue
        lab=etiqueta(cf, mp.loc[s, cf["col"]])
        if lab is None: continue
        keep.append(s); y.append(lab)
    return X.loc[keep], np.array(y)
def clr(df):
    rel=df.div(df.sum(axis=1),axis=0).replace(0,PSEUDO)
    lg=np.log(rel); return lg.sub(lg.mean(axis=1),axis=0)
def auc_cv(X,y,nrep=NREP,k=10):
    if len(np.unique(y))<2 or min(np.bincount(y))<5: return np.nan
    out=[]
    for r in range(nrep):
        kk=min(k, min(np.bincount(y)))
        skf=StratifiedKFold(kk,shuffle=True,random_state=SEED+r); pred=np.zeros(len(y))
        for tr,te in skf.split(X,y):
            m=RandomForestClassifier(n_estimators=500,max_features="sqrt",
                class_weight="balanced",random_state=SEED+r,n_jobs=-1)
            m.fit(X.iloc[tr],y[tr]); pred[te]=m.predict_proba(X.iloc[te])[:,1]
        out.append(roc_auc_score(y,pred))
    return float(np.mean(out))
resumen=[]; detalle=[]
for brazo, COH in BRAZOS.items():
    log("\n"+"="*86); log(f"BRAZO: {brazo}"); log("="*86)
    crudo={}
    for n,cf in COH.items():
        X,y=cargar(cf); crudo[n]=(X,y)
        log(f"  {n:<11} n={len(y):>4}  casos={y.sum():>3}  ctrl={(y==0).sum():>3}")
    comunes=sorted(set.intersection(*[set(X.columns) for X,_ in crudo.values()]))
    feats=[g for g in comunes if min((X[g]>0).mean() for X,_ in crudo.values())>=0.20]
    log(f"\n  Descriptores comunes tras prevalencia >=20%: {len(feats)}")
    CLR={n:(clr(X[feats]),y) for n,(X,y) in crudo.items()}
    noms=list(CLR)
    log(f"\n  A. CASO vs CONTROL dentro de cada cohorte")
    ints=[]
    for n,(Xc,y) in CLR.items():
        a=auc_cv(Xc,y); ints.append(a)
        log(f"     {n:<11} {a:.3f}")
        detalle.append(dict(brazo=brazo, prueba="A_caso_vs_control", grupo=n, auc=a))
    log(f"     media: {np.nanmean(ints):.3f}")
    log(f"\n  B/C. ORIGEN de la muestra entre cohortes")
    log(f"     {'Par':<26}{'CONTROLES':>12}{'CASOS':>10}{'dif':>9}")
    bs, cs = [], []
    for i,a in enumerate(noms):
        for b in noms[i+1:]:
            Xa,ya=CLR[a]; Xb,yb=CLR[b]
            Ca=Xa[ya==0]; Cb=Xb[yb==0]
            R=pd.concat([Ca,Cb],ignore_index=True)
            yy=np.array([0]*len(Ca)+[1]*len(Cb))
            aB=auc_cv(R,yy)
            Ka=Xa[ya==1]; Kb=Xb[yb==1]
            R2=pd.concat([Ka,Kb],ignore_index=True)
            yy2=np.array([0]*len(Ka)+[1]*len(Kb))
            aC=auc_cv(R2,yy2)
            bs.append(aB); cs.append(aC)
            dif = aC-aB if (not np.isnan(aB) and not np.isnan(aC)) else np.nan
            log(f"     {a+' vs '+b:<26}{aB:>12.3f}{aC:>10.3f}{dif:>+9.3f}")
            detalle.append(dict(brazo=brazo, prueba="B_origen_controles", grupo=f"{a}|{b}", auc=aB))
            detalle.append(dict(brazo=brazo, prueba="C_origen_casos", grupo=f"{a}|{b}", auc=aC))
    mB, mC = np.nanmean(bs), np.nanmean(cs)
    log(f"     {'media':<26}{mB:>12.3f}{mC:>10.3f}{mC-mB:>+9.3f}")
    resumen.append(dict(brazo=brazo,
        A_caso_vs_control=round(float(np.nanmean(ints)),3),
        B_origen_controles=round(float(mB),3),
        C_origen_casos=round(float(mC),3),
        dif_C_menos_B=round(float(mC-mB),3)))
log("\n"+"="*86); log("RESUMEN: GEOGRAFIA FRENTE A CONDICION"); log("="*86)
R=pd.DataFrame(resumen)
log(R.to_string(index=False))
log("")
log("  A = distinguir caso de control DENTRO de una cohorte")
log("  B = distinguir el origen de dos sujetos SANOS de cohortes distintas")
log("  C = distinguir el origen de dos sujetos ENFERMOS de cohortes distintas")
R.to_csv(f"{OUT}/geografia_vs_condicion_resumen.csv", index=False)
pd.DataFrame(detalle).to_csv(f"{OUT}/geografia_vs_condicion_detalle.csv", index=False)
log(f"\nSalidas en {OUT}")
log("="*86)
