#!/usr/bin/env python3
"""Control de efecto de lote: son las cohortes tecnicamente distinguibles?"""
import os, warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
SEED = 42; RNG = np.random.default_rng(42)
BASE = os.path.expanduser("~/Documents/MicrobiomeRadar_Research")
AD = os.path.join(BASE, "03_Alzheimer_Neurodegenerative")
OUT = os.path.join(BASE, "_crosscohort")
COHORTES = {
 "Li_Shanghai": dict(
   matriz=f"{AD}/03d_Li2019_Shanghai_China/02_processed/Li_genus_counts.csv",
   mapa=f"{AD}/03d_Li2019_Shanghai_China/00_screening/mapa_heces.tsv",
   key="run_accession", grpcol="grupo", modo="exacto", caso=["AD"], ctrl=["NC"],
   track=f"{AD}/03d_Li2019_Shanghai_China/02_processed/Li_track.csv"),
 "Ling_Lishui": dict(
   matriz=f"{AD}/03g_Ling2021_Lishui_China/02_processed/Ling_genus_counts.csv",
   mapa=f"{AD}/03g_Ling2021_Lishui_China/00_screening/ena_PRJNA633959.tsv",
   key="run_accession", grpcol="library_name", modo="prefijo",
   track=f"{AD}/03g_Ling2021_Lishui_China/02_processed/Ling_track.csv"),
 "Zhuang_Chongqing": dict(
   matriz=f"{AD}/03h_Zhuang2019_Chongqing_China/02_processed/Zhuang_genus_counts.csv",
   mapa=f"{AD}/03h_Zhuang2019_Chongqing_China/00_screening/meta_554111.tsv",
   key="run_accession", grpcol="sample_alias", modo="prefijo",
   track=f"{AD}/03h_Zhuang2019_Chongqing_China/02_processed/Zhuang_track.csv"),
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
    return X.loc[keep], np.array(y)
def cv_auc(X, y, nrep=10, k=10):
    out=[]
    for r in range(nrep):
        skf = StratifiedKFold(k, shuffle=True, random_state=SEED+r); pred=np.zeros(len(y))
        for tr,te in skf.split(X,y):
            m = RandomForestClassifier(n_estimators=500, max_features="sqrt",
                class_weight="balanced", random_state=SEED+r, n_jobs=-1)
            m.fit(X.iloc[tr],y[tr]); pred[te]=m.predict_proba(X.iloc[te])[:,1]
        out.append(roc_auc_score(y,pred))
    return np.array(out)
def cliff(a,b):
    a,b=np.asarray(a),np.asarray(b); n=len(a)*len(b)
    return (sum((x>b).sum() for x in a)-sum((x<b).sum() for x in a))/n if n else np.nan

log("="*84); log("CONTROL DE EFECTO DE LOTE"); log("="*84)
crudo = {}
for n in COHORTES:
    X,y = cargar(n); crudo[n]=(X,y)
    log(f"  {n:<20} n={len(y):>4}  generos={X.shape[1]}")

# ---- 1. profundidad vs etiqueta DENTRO de cada cohorte ----
log("\n"+"="*84)
log("1. PROFUNDIDAD DE SECUENCIACION vs ETIQUETA (dentro de cada cohorte)")
log("="*84)
log(f"{'Cohorte':<22}{'med AD':>10}{'med ctrl':>10}{'p':>10}{'delta':>8}   veredicto")
for n,(X,y) in crudo.items():
    prof = X.sum(axis=1).values
    a, c = prof[y==1], prof[y==0]
    _, p = stats.mannwhitneyu(a, c, alternative="two-sided")
    d = cliff(a, c)
    v = ">>> CONFUSOR <<<" if p < 0.05 and abs(d) > 0.2 else "ok"
    log(f"{n:<22}{np.median(a):>10.0f}{np.median(c):>10.0f}{p:>10.4f}{d:>8.2f}   {v}")

# ---- 2. riqueza vs etiqueta ----
log("\n"+"="*84)
log("2. RIQUEZA (n generos observados) vs ETIQUETA")
log("="*84)
log(f"{'Cohorte':<22}{'AD':>10}{'ctrl':>10}{'p':>10}{'delta':>8}")
for n,(X,y) in crudo.items():
    riq = (X>0).sum(axis=1).values
    a, c = riq[y==1], riq[y==0]
    _, p = stats.mannwhitneyu(a, c, alternative="two-sided")
    log(f"{n:<22}{np.median(a):>10.0f}{np.median(c):>10.0f}{p:>10.4f}{cliff(a,c):>8.2f}")

# ---- 3. clasificador de COHORTE (ignorando diagnostico) ----
log("\n"+"="*84)
log("3. CLASIFICADOR DE COHORTE — son tecnicamente distinguibles?")
log("="*84)
comunes = sorted(set.intersection(*[set(X.columns) for X,_ in crudo.values()]))
log(f"  generos comunes: {len(comunes)}")
frames, lote = [], []
for i,(n,(X,y)) in enumerate(crudo.items()):
    rel = X[comunes].div(X[comunes].sum(axis=1), axis=0).replace(0,PSEUDO)
    lg = np.log(rel); frames.append(lg.sub(lg.mean(axis=1),axis=0))
    lote += [i]*len(X)
XA = pd.concat(frames, ignore_index=True); lote = np.array(lote)
noms = list(crudo)
log(f"\n{'Par de cohortes':<46}{'AUC':>8}   interpretacion")
for i in range(len(noms)):
    for j in range(i+1, len(noms)):
        m = (lote==i)|(lote==j)
        yy = (lote[m]==j).astype(int)
        a = cv_auc(XA[m].reset_index(drop=True), yy, nrep=5)
        itp = ("SEPARACION TOTAL" if a.mean()>0.95 else
               "muy distinguibles" if a.mean()>0.85 else
               "distinguibles" if a.mean()>0.70 else "similares")
        log(f"{noms[i]+' vs '+noms[j]:<46}{a.mean():>8.3f}   {itp}")

# ---- 4. comparacion: separacion por COHORTE vs por DIAGNOSTICO ----
log("\n"+"="*84)
log("4. QUE SEPARA MAS: la cohorte o el diagnostico?")
log("="*84)
log(f"{'Comparacion':<46}{'AUC':>8}")
for i,(n,(X,y)) in enumerate(crudo.items()):
    rel = X[comunes].div(X[comunes].sum(axis=1),axis=0).replace(0,PSEUDO)
    lg = np.log(rel); Xc = lg.sub(lg.mean(axis=1),axis=0)
    a = cv_auc(Xc, y, nrep=5)
    log(f"{'DIAGNOSTICO dentro de '+n:<46}{a.mean():>8.3f}")

# ---- 5. controles de una cohorte vs controles de otra ----
log("\n"+"="*84)
log("5. CONTROLES vs CONTROLES entre cohortes (deberian parecerse)")
log("="*84)
log(f"{'Par':<46}{'AUC ctrl-ctrl':>15}")
for i in range(len(noms)):
    for j in range(i+1, len(noms)):
        Xi, yi = crudo[noms[i]]; Xj, yj = crudo[noms[j]]
        ri = Xi[comunes][yi==0]; rj = Xj[comunes][yj==0]
        R = pd.concat([ri, rj], ignore_index=True)
        rel = R.div(R.sum(axis=1),axis=0).replace(0,PSEUDO)
        lg = np.log(rel); Rc = lg.sub(lg.mean(axis=1),axis=0)
        yy = np.array([0]*len(ri)+[1]*len(rj))
        a = cv_auc(Rc, yy, nrep=5)
        log(f"{noms[i]+' vs '+noms[j]:<46}{a.mean():>15.3f}")

log("\n"+"="*84)
log("LECTURA: si (5) da AUC > 0.9, los controles sanos de dos ciudades son")
log("mas distinguibles entre si que enfermos de sanos dentro de una ciudad.")
log("Eso significa que el efecto cohorte domina sobre el efecto enfermedad.")
log("="*84)
