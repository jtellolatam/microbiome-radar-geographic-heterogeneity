#!/usr/bin/env python3
"""m2 - Sensibilidad al pseudoconteo en las once cohortes"""
import os, warnings
import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
SEED=42
BASE=os.path.expanduser("~/Documents/MicrobiomeRadar_Research")
AD=f"{BASE}/03_Alzheimer_Neurodegenerative"
TEA=f"{BASE}/01_TEA_Neurodivergent"
PD=f"{BASE}/02_Parkinson_Neurodegenerative"
OUT=f"{BASE}/_crosscohort"
VALORES=[1e-7, 1e-6, 1e-5, 1e-4]
COH = {
 "Shanghai": dict(m=f"{AD}/03d_Li2019_Shanghai_China/02_processed/Li_genus_counts.csv",
   p=f"{AD}/03d_Li2019_Shanghai_China/00_screening/mapa_heces.tsv",
   key="run_accession", col="grupo", modo="exacto", caso=["AD"], ctrl=["NC"], brazo="AD"),
 "Lishui": dict(m=f"{AD}/03g_Ling2021_Lishui_China/02_processed/Ling_genus_counts.csv",
   p=f"{AD}/03g_Ling2021_Lishui_China/00_screening/ena_PRJNA633959.tsv",
   key="run_accession", col="library_name", modo="prefijo", brazo="AD"),
 "Chongqing": dict(m=f"{AD}/03h_Zhuang2019_Chongqing_China/02_processed/Zhuang_genus_counts.csv",
   p=f"{AD}/03h_Zhuang2019_Chongqing_China/00_screening/meta_554111.tsv",
   key="run_accession", col="sample_alias", modo="prefijo", brazo="AD"),
 "Estambul": dict(m=f"{AD}/03i_Yildirim2022_Turkey/02_processed/Turquia_genus_counts.csv",
   p=f"{AD}/03i_Yildirim2022_Turkey/00_screening/ena_PRJNA734525.tsv",
   key="run_accession", col="library_name", modo="turquia", centro="Estambul", brazo="AD"),
 "Kayseri": dict(m=f"{AD}/03i_Yildirim2022_Turkey/02_processed/Turquia_genus_counts.csv",
   p=f"{AD}/03i_Yildirim2022_Turkey/00_screening/ena_PRJNA734525.tsv",
   key="run_accession", col="library_name", modo="turquia", centro="Kayseri", brazo="AD"),
 "Brescia": dict(m=f"{AD}/03k_Italia_PRJEB89799/02_processed/Italia_genus_counts.csv",
   p=f"{AD}/03k_Italia_PRJEB89799/00_screening/mapa_baseline.tsv",
   key="run", col="grupo", modo="exacto", caso=["AD"], ctrl=["Control"], brazo="AD"),
 "Nanjing": dict(m=f"{TEA}/01f_Dan_Nanjing/02_processed/Dan_genus_counts.csv",
   p=f"{TEA}/01f_Dan_Nanjing/00_screening/mapa_dan.tsv",
   key="run", col="grupo", modo="exacto", caso=["ASD"], ctrl=["Control"], brazo="TEA"),
 "Xiamen": dict(m=f"{TEA}/01g_Xiamen/02_processed/Xiamen_genus_counts.csv",
   p=f"{TEA}/01g_Xiamen/00_screening/mapa_xiamen.tsv",
   key="run", col="grupo", modo="exacto", caso=["ASD"], ctrl=["Control"], brazo="TEA"),
 "Caltech": dict(m=f"{TEA}/01h_Caltech_USA/02_processed/Caltech_genus_counts.csv",
   p=f"{TEA}/01h_Caltech_USA/00_screening/mapa_caltech.tsv",
   key="run", col="grupo", modo="exacto", caso=["ASD"], ctrl=["Control"], brazo="TEA"),
 "Finlandia": dict(m=f"{PD}/02e_Finlandia/02_processed/Finlandia_genus_counts.csv",
   p=f"{PD}/02e_Finlandia/00_screening/mapa_finlandia.tsv",
   key="run", col="grupo", modo="exacto", caso=["PD"], ctrl=["Control"], brazo="PD"),
 "Rumania": dict(m=f"{PD}/02f_Rumania/02_processed/Rumania_genus_counts.csv",
   p=f"{PD}/02f_Rumania/00_screening/mapa_rumania.tsv",
   key="muestra", col="grupo", modo="exacto", caso=["PD"], ctrl=["Control"], brazo="PD"),
}

def etiqueta(cf,v):
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

def cargar(n):
    cf=COH[n]; X=pd.read_csv(cf["m"], index_col=0)
    mp=pd.read_csv(cf["p"], sep="\t", dtype=str)
    mp=mp[~mp[cf["key"]].duplicated(keep="first")].set_index(cf["key"])
    y,keep=[],[]
    for s in X.index:
        if s not in mp.index: continue
        lab=etiqueta(cf, mp.loc[s, cf["col"]])
        if lab is None: continue
        keep.append(s); y.append(lab)
    return X.loc[keep], np.array(y)

def clr(df, pc):
    rel=df.div(df.sum(axis=1),axis=0).replace(0,pc)
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
    return float(np.mean(out))

print("="*80)
print("m2 · SENSIBILIDAD AL PSEUDOCONTEO")
print("="*80)
print(f"{'Cohorte':<12}{'Brazo':<7}" + "".join(f"{f'pc={v:g}':>12}" for v in VALORES) + f"{'rango':>9}")
print("-"*80)
filas=[]
for n,cf in COH.items():
    X,y = cargar(n)
    feats=[c for c in X.columns if (X[c]>0).mean()>=0.20]
    vals=[auc_cv(clr(X[feats],pc), y) for pc in VALORES]
    rng=max(vals)-min(vals)
    filas.append(dict(cohorte=n, brazo=cf["brazo"],
        **{f"auc_pc_{v:g}":round(a,3) for v,a in zip(VALORES,vals)}, rango=round(rng,3)))
    print(f"{n:<12}{cf['brazo']:<7}" + "".join(f"{a:>12.3f}" for a in vals) + f"{rng:>9.3f}")
R=pd.DataFrame(filas)
R.to_csv(f"{OUT}/m2_sensibilidad_escala.csv", index=False)
cols=[c for c in R.columns if c.startswith("auc_pc_")]
print("-"*80)
print(f"{'MEDIA':<19}" + "".join(f"{R[c].mean():>12.3f}" for c in cols) + f"{R['rango'].mean():>9.3f}")
print(f"\nRango maximo entre valores de pseudoconteo: {R['rango'].max():.3f}")
print(f"Guardado en {OUT}/m2_sensibilidad_escala.csv")
print("="*80)
