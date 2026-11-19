#!/usr/bin/env python3
"""Microbiome Radar - Modelado y matriz de transferencia (bloque China)"""
import os, json, warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")
RNG = np.random.default_rng(42)
SEED = 42
BASE = os.path.expanduser("~/Documents/MicrobiomeRadar_Research")
AD = os.path.join(BASE, "03_Alzheimer_Neurodegenerative")
OUT = os.path.join(BASE, "_crosscohort")
os.makedirs(os.path.join(OUT, "transfer_matrix"), exist_ok=True)

COHORTES = {
 "Li_Shanghai": dict(
   matriz=f"{AD}/03d_Li2019_Shanghai_China/02_processed/Li_genus_counts.csv",
   mapa=f"{AD}/03d_Li2019_Shanghai_China/00_screening/mapa_heces.tsv",
   key="run_accession", grpcol="grupo", modo="exacto",
   caso=["AD"], ctrl=["NC"], edad=None),
 "Ling_Lishui": dict(
   matriz=f"{AD}/03g_Ling2021_Lishui_China/02_processed/Ling_genus_counts.csv",
   mapa=f"{AD}/03g_Ling2021_Lishui_China/00_screening/ena_PRJNA633959.tsv",
   key="run_accession", grpcol="library_name", modo="prefijo",
   caso=None, ctrl=None, edad=None),
 "Zhuang_Chongqing": dict(
   matriz=f"{AD}/03h_Zhuang2019_Chongqing_China/02_processed/Zhuang_genus_counts.csv",
   mapa=f"{AD}/03h_Zhuang2019_Chongqing_China/00_screening/meta_554111.tsv",
   key="run_accession", grpcol="sample_alias", modo="prefijo",
   caso=None, ctrl=None,
   edad=f"{AD}/03h_Zhuang2019_Chongqing_China/00_screening/meta_554111_full.tsv"),
}
PREV_MIN, PSEUDO, N_REPEATS, N_PERM, N_BOOT = 0.20, 1e-6, 20, 300, 1000
def log(*a): print(*a, flush=True)

def etiqueta(cf, valor):
    v = str(valor).upper()
    if cf["modo"] == "exacto":
        if valor in cf["caso"]: return 1
        if valor in cf["ctrl"]: return 0
        return None
    if v.startswith("CON") or v.startswith("C"): return 0
    if v.startswith("AD"): return 1
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
    X = X.loc[keep]; y = np.array(y)
    log(f"  {nom:<20} n={len(y):>4}  casos={y.sum():>3}  ctrl={(y==0).sum():>3}  gen={X.shape[1]}")
    return X, y

def clr(df):
    rel = df.div(df.sum(axis=1), axis=0).replace(0, PSEUDO)
    lg = np.log(rel)
    return lg.sub(lg.mean(axis=1), axis=0)

def cv_auc(X, y, nrep=N_REPEATS, ntree=500):
    out = []
    for r in range(nrep):
        skf = StratifiedKFold(5, shuffle=True, random_state=SEED+r)
        pred = np.zeros(len(y))
        for tr, te in skf.split(X, y):
            m = RandomForestClassifier(n_estimators=ntree, random_state=SEED+r,
                                       n_jobs=-1, class_weight="balanced")
            m.fit(X.iloc[tr], y[tr]); pred[te] = m.predict_proba(X.iloc[te])[:,1]
        out.append(roc_auc_score(y, pred))
    return np.array(out)

def perm_test(X, y, real, n=N_PERM):
    nul = []
    for _ in range(n):
        yp = RNG.permutation(y)
        skf = StratifiedKFold(5, shuffle=True, random_state=SEED)
        pred = np.zeros(len(yp))
        for tr, te in skf.split(X, yp):
            m = RandomForestClassifier(n_estimators=100, random_state=SEED,
                                       n_jobs=-1, class_weight="balanced")
            m.fit(X.iloc[tr], yp[tr]); pred[te] = m.predict_proba(X.iloc[te])[:,1]
        nul.append(roc_auc_score(yp, pred))
    nul = np.array(nul)
    return (np.sum(nul >= real) + 1) / (n + 1)

def ic_boot(yt, ys, n=N_BOOT):
    v, idx = [], np.arange(len(yt))
    for _ in range(n):
        s = RNG.choice(idx, len(idx), replace=True)
        if len(np.unique(yt[s])) < 2: continue
        v.append(roc_auc_score(yt[s], ys[s]))
    return np.percentile(v, [2.5, 97.5])

log("="*74); log("CARGA"); log("="*74)
datos = {}
for n in COHORTES:
    try: datos[n] = cargar(n)
    except Exception as e: log(f"  !! {n}: {e}")
if len(datos) < 2: raise SystemExit("Se necesitan >=2 cohortes")

comunes = set.intersection(*[set(X.columns) for X,_ in datos.values()])
log(f"\nGeneros comunes: {len(comunes)}")
feats = sorted([g for g in comunes
                if min((X[g]>0).mean() for X,_ in datos.values()) >= PREV_MIN])
log(f"Tras prevalencia >={PREV_MIN:.0%}: {len(feats)} generos")
CLR = {n: (clr(X[feats]), y) for n,(X,y) in datos.items()}

log("\n"+"="*74); log(f"1. AUC INTERNO (5-fold x {N_REPEATS})"); log("="*74)
log(f"{'Cohorte':<22}{'AUC':>8}{'DE':>8}{'IC95':>20}{'p_perm':>9}")
internos = {}
for n,(Xc,y) in CLR.items():
    a = cv_auc(Xc, y); m, sd = a.mean(), a.std()
    skf = StratifiedKFold(5, shuffle=True, random_state=SEED)
    pred = np.zeros(len(y))
    for tr,te in skf.split(Xc,y):
        md = RandomForestClassifier(n_estimators=500, random_state=SEED,
                                    n_jobs=-1, class_weight="balanced")
        md.fit(Xc.iloc[tr], y[tr]); pred[te] = md.predict_proba(Xc.iloc[te])[:,1]
    lo,hi = ic_boot(y, pred); p = perm_test(Xc, y, m)
    internos[n] = dict(auc=m, sd=sd, ic=[lo,hi], p=p)
    log(f"{n:<22}{m:>8.3f}{sd:>8.3f}    [{lo:.3f}, {hi:.3f}]{p:>9.4f}")

log("\n"+"="*74); log("2. LINEA BASE DEMOGRAFICA"); log("="*74)
for n, cf in COHORTES.items():
    if not cf.get("edad") or not os.path.exists(cf["edad"] or "") or n not in CLR:
        log(f"  {n:<22} sin edad -> NO EVALUABLE"); continue
    dem = pd.read_csv(cf["edad"], sep="\t", dtype=str).set_index("sample_name")
    mp = pd.read_csv(cf["mapa"], sep="\t", dtype=str).set_index(cf["key"])
    Xc, y = CLR[n]
    ed, sx, yy, idx = [], [], [], []
    for i, run in enumerate(Xc.index):
        al = mp.loc[run, cf["grpcol"]]
        if al in dem.index and str(dem.loc[al,"age"]).isdigit():
            ed.append(int(dem.loc[al,"age"]))
            sx.append(1 if str(dem.loc[al,"sex"]).lower().startswith("m") else 0)
            yy.append(y[i]); idx.append(i)
    if len(ed) < 20: log(f"  {n}: solo {len(ed)} con edad"); continue
    yy = np.array(yy); D = pd.DataFrame({"edad":ed,"sexo":sx})
    Xs = Xc.iloc[idx].reset_index(drop=True)
    log(f"\n  {n}  (n={len(yy)})")
    for k, M in [("solo demografia",D), ("solo microbioma",Xs),
                 ("ambos", pd.concat([Xs,D],axis=1))]:
        a = cv_auc(M, yy, nrep=10)
        log(f"    {k:<20} AUC = {a.mean():.3f} +/- {a.std():.3f}")

log("\n"+"="*74); log("3. MATRIZ DE TRANSFERENCIA (fila=train, col=test)"); log("="*74)
noms = list(CLR.keys())
M = pd.DataFrame(index=noms, columns=noms, dtype=float); det = []
for a in noms:
    Xa, ya = CLR[a]
    md = RandomForestClassifier(n_estimators=500, random_state=SEED,
                                n_jobs=-1, class_weight="balanced").fit(Xa, ya)
    for b in noms:
        if a == b: M.loc[a,b] = internos[a]["auc"]; continue
        Xb, yb = CLR[b]
        sc = md.predict_proba(Xb[Xa.columns])[:,1]
        auc = roc_auc_score(yb, sc); lo,hi = ic_boot(yb, sc)
        M.loc[a,b] = auc
        det.append(dict(train=a, test=b, auc=auc, ic_lo=lo, ic_hi=hi))
log("\n"+M.round(3).to_string())
M.to_csv(f"{OUT}/transfer_matrix/matriz_transferencia_china.csv")
pd.DataFrame(det).to_csv(f"{OUT}/transfer_matrix/transferencia_detalle.csv", index=False)
diag = np.mean([M.loc[n,n] for n in noms])
fu = [M.loc[a,b] for a in noms for b in noms if a!=b]
log(f"\n  AUC interno medio : {diag:.3f}")
log(f"  AUC cruzado medio : {np.mean(fu):.3f}")
log(f"  CAIDA             : {diag-np.mean(fu):.3f}")
log(f"  Pares < 0.5 (inversion): {sum(1 for v in fu if v<0.5)}/{len(fu)}")

log("\n"+"="*74); log("4. CONCORDANCIA DE RANKINGS (H2)"); log("="*74)
try:
    import shap; usar = True
except ImportError:
    log("  shap no instalado -> importancia de RF"); usar = False
imp = {}
for n,(Xc,y) in CLR.items():
    md = RandomForestClassifier(n_estimators=500, random_state=SEED,
                                n_jobs=-1, class_weight="balanced").fit(Xc, y)
    if usar:
        sv = shap.TreeExplainer(md).shap_values(Xc)
        if isinstance(sv, list): sv = sv[1]
        if getattr(sv,'ndim',2) == 3: sv = sv[:,:,1]
        imp[n] = pd.Series(np.abs(sv).mean(axis=0), index=Xc.columns)
    else:
        imp[n] = pd.Series(md.feature_importances_, index=Xc.columns)
log(f"\n{'Par':<46}{'Spearman':>10}{'p':>10}")
for i,a in enumerate(noms):
    for b in noms[i+1:]:
        rho,p = stats.spearmanr(imp[a], imp[b])
        log(f"{a+' vs '+b:<46}{rho:>10.3f}{p:>10.4f}")
top = pd.DataFrame(imp); top["media"] = top.mean(axis=1)
top = top.sort_values("media", ascending=False)
log(f"\nTop 15 por importancia media:\n{top.head(15).round(4).to_string()}")
top.to_csv(f"{OUT}/importancia_china.csv")
json.dump({k:{kk:(vv.tolist() if isinstance(vv,np.ndarray) else vv)
              for kk,vv in v.items()} for k,v in internos.items()},
          open(f"{OUT}/auc_internos_china.json","w"), indent=2, default=float)
log("\n"+"="*74); log("LISTO. Salidas en "+OUT); log("="*74)
