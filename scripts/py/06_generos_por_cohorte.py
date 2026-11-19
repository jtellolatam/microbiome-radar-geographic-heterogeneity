#!/usr/bin/env python3
"""Que generos sostienen la senal en cada cohorte (todos los generos, SHAP 20 semillas)"""
import os, warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
SEED = 42
BASE = os.path.expanduser("~/Documents/MicrobiomeRadar_Research")
AD = os.path.join(BASE, "03_Alzheimer_Neurodegenerative")
OUT = os.path.join(BASE, "_crosscohort")
GENEROS_TFM = ["Faecalibacterium","Bifidobacterium","Akkermansia","Lactobacillus",
               "Blautia","Clostridium","Sutterella","Bacteroides","Veillonella","Prevotella"]
TABLA_AD = ["Bacteroides","Bifidobacterium","Prevotella","Escherichia"]
COHORTES = {
 "Li_Shanghai": dict(
   matriz=f"{AD}/03d_Li2019_Shanghai_China/02_processed/Li_genus_counts.csv",
   mapa=f"{AD}/03d_Li2019_Shanghai_China/00_screening/mapa_heces.tsv",
   key="run_accession", grpcol="grupo", modo="exacto", caso=["AD"], ctrl=["NC"]),
 "Ling_Lishui": dict(
   matriz=f"{AD}/03g_Ling2021_Lishui_China/02_processed/Ling_genus_counts.csv",
   mapa=f"{AD}/03g_Ling2021_Lishui_China/00_screening/ena_PRJNA633959.tsv",
   key="run_accession", grpcol="library_name", modo="prefijo"),
 "Zhuang_Chongqing": dict(
   matriz=f"{AD}/03h_Zhuang2019_Chongqing_China/02_processed/Zhuang_genus_counts.csv",
   mapa=f"{AD}/03h_Zhuang2019_Chongqing_China/00_screening/meta_554111.tsv",
   key="run_accession", grpcol="sample_alias", modo="prefijo"),
}
PSEUDO = 1e-6
def log(*a): print(*a, flush=True)
def etiqueta(cf, v):
    s = str(v).upper()
    if cf["modo"]=="exacto":
        if v in cf["caso"]: return 1
        if v in cf["ctrl"]: return 0
        return None
    if s.startswith("CON") or s.startswith("C"): return 0
    if s.startswith("AD"): return 1
    return None
def cargar(nom):
    cf = COHORTES[nom]
    X = pd.read_csv(cf["matriz"], index_col=0)
    mp = pd.read_csv(cf["mapa"], sep="\t", dtype=str).set_index(cf["key"])
    y, keep = [], []
    for run in X.index:
        if run not in mp.index: continue
        lab = etiqueta(cf, mp.loc[run, cf["grpcol"]])
        if lab is None: continue
        keep.append(run); y.append(lab)
    X = X.loc[keep]; rel = X.div(X.sum(axis=1), axis=0)
    prev = (rel>0).mean(); cols = prev[prev>=0.20].index
    sub = rel[cols].replace(0,PSEUDO); lg = np.log(sub)
    return lg.sub(lg.mean(axis=1),axis=0), np.array(y), rel
def cv_auc(X, y, nrep=10, k=10):
    out=[]
    for r in range(nrep):
        skf = StratifiedKFold(k, shuffle=True, random_state=SEED+r); pred=np.zeros(len(y))
        for tr,te in skf.split(X,y):
            m = RandomForestClassifier(n_estimators=500, max_features="sqrt",
                class_weight="balanced", random_state=SEED+r, n_jobs=-1)
            m.fit(X.iloc[tr],y[tr]); pred[te]=m.predict_proba(X.iloc[te])[:,1]
        out.append(roc_auc_score(y,pred))
    return np.mean(out)
def cliff(a,b):
    a,b=np.asarray(a),np.asarray(b); n=len(a)*len(b)
    return (sum((x>b).sum() for x in a)-sum((x<b).sum() for x in a))/n if n else np.nan
def es_tfm(g): return any(g==t or g.startswith(t+"_") or g.startswith(t+"-") for t in GENEROS_TFM)
def es_tabla(g): return any(g==t or g.startswith(t+"_") or g.startswith(t+"-") for t in TABLA_AD)
try:
    import shap; USAR=True
except ImportError: USAR=False
log("="*88); log("QUE GENEROS SOSTIENEN LA SENAL EN CADA COHORTE"); log("="*88)
datos, IMP = {}, {}
for n in COHORTES:
    X,y,rel = cargar(n); datos[n]=(X,y,rel)
    log(f"  {n:<20} n={len(y):>4}  generos={X.shape[1]}")
    reps=[]
    for r in range(20):
        md = RandomForestClassifier(n_estimators=500, max_features="sqrt",
            class_weight="balanced", random_state=SEED+r, n_jobs=-1).fit(X,y)
        if USAR:
            sv = shap.TreeExplainer(md).shap_values(X)
            if isinstance(sv,list): sv=sv[1]
            if getattr(sv,"ndim",2)==3: sv=sv[:,:,1]
            reps.append(pd.Series(np.abs(sv).mean(axis=0), index=X.columns))
        else: reps.append(pd.Series(md.feature_importances_, index=X.columns))
    IMP[n]=pd.concat(reps,axis=1).mean(axis=1).sort_values(ascending=False)
for n in COHORTES:
    X,y,rel = datos[n]
    log("\n"+"="*88); log(f"TOP 20 - {n}"); log("="*88)
    log(f"{'#':<4}{'Genero':<34}{'SHAP':>9}{'delta':>8}{'p':>9}   origen")
    for i,(g,v) in enumerate(IMP[n].head(20).items(),1):
        a=rel.loc[y==1,g].values; c=rel.loc[y==0,g].values
        try: _,p = stats.mannwhitneyu(a,c,alternative="two-sided")
        except ValueError: p=np.nan
        d = cliff(a,c)
        org = "TABLA AD" if es_tabla(g) else ("set TFM" if es_tfm(g) else "NUEVO")
        log(f"{i:<4}{g[:33]:<34}{v:>9.4f}{d:>8.2f}{p:>9.4f}   {org}")
    nn = sum(1 for g in IMP[n].head(20).index if not es_tfm(g) and not es_tabla(g))
    log(f"\n  Nuevos en el top 20: {nn}/20")
log("\n"+"="*88); log("AUC SEGUN NUMERO DE FEATURES"); log("="*88)
log(f"{'Cohorte':<22}{'top5':>8}{'top10':>8}{'top20':>8}{'top30':>8}{'todos':>8}")
res={}
for n in COHORTES:
    X,y,_ = datos[n]; fila=[]
    for k in (5,10,20,30): fila.append(cv_auc(X[list(IMP[n].head(k).index)], y))
    fila.append(cv_auc(X,y)); res[n]=fila
    log(f"{n:<22}"+"".join(f"{v:>8.3f}" for v in fila))
log("\n"+"="*88); log("SOLAPAMIENTO DE LOS TOP 20"); log("="*88)
noms=list(COHORTES); tops={n:set(IMP[n].head(20).index) for n in noms}
for i,a in enumerate(noms):
    for b in noms[i+1:]:
        it = tops[a]&tops[b]
        log(f"{a} & {b}: {len(it)}/20")
        if it: log(f"    {', '.join(sorted(it))}")
cm = set.intersection(*tops.values())
log(f"\nComun a las TRES: {len(cm)}/20")
if cm: log(f"  {', '.join(sorted(cm))}")
pd.DataFrame(IMP).to_csv(f"{OUT}/shap_todos_generos_por_cohorte.csv")
pd.DataFrame(res, index=["top5","top10","top20","top30","todos"]).T.to_csv(
    f"{OUT}/auc_por_numero_features.csv")
log(f"\nSalidas en {OUT}")
