#!/usr/bin/env python3
"""EXPERIMENTO DECISIVO: Estambul vs Kayseri — dos ciudades, un solo estudio"""
import os, warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
SEED=42; RNG=np.random.default_rng(42)
BASE=os.path.expanduser("~/Documents/MicrobiomeRadar_Research")
TR=os.path.join(BASE,"03_Alzheimer_Neurodegenerative","03i_Yildirim2022_Turkey")
OUT=os.path.join(BASE,"_crosscohort"); PSEUDO=1e-6
def log(*a): print(*a, flush=True)
X = pd.read_csv(f"{TR}/02_processed/Turquia_genus_counts.csv", index_col=0)
mp = pd.read_csv(f"{TR}/00_screening/ena_PRJNA734525.tsv", sep="\t", dtype=str).set_index("run_accession")
filas=[]
for run in X.index:
    if run not in mp.index: continue
    lib=str(mp.loc[run,"library_name"]); pref="".join(c for c in lib if not c.isdigit())
    if pref in ("ADH","CK"): centro,grupo="Kayseri",("AD" if pref=="ADH" else "C")
    elif pref in ("AD","C"): centro,grupo="Estambul",pref
    elif pref=="MCI": centro,grupo="Estambul","MCI"
    else: continue
    filas.append((run,centro,grupo))
M=pd.DataFrame(filas,columns=["run","centro","grupo"]).set_index("run"); X=X.loc[M.index]
log("="*80); log("ESTAMBUL vs KAYSERI - dos ciudades, UN solo estudio"); log("="*80)
log("\n"+pd.crosstab(M.centro,M.grupo).to_string())
rel=X.div(X.sum(axis=1),axis=0); prev=(rel>0).mean(); cols=prev[prev>=0.20].index
sub=rel[cols].replace(0,PSEUDO); lg=np.log(sub); XC=lg.sub(lg.mean(axis=1),axis=0)
log(f"\nGeneros tras prevalencia >=20%: {XC.shape[1]}")
def cv_auc(Xd,y,nrep=10,k=10):
    out=[]
    for r in range(nrep):
        skf=StratifiedKFold(k,shuffle=True,random_state=SEED+r); pred=np.zeros(len(y))
        for tr,te in skf.split(Xd,y):
            m=RandomForestClassifier(n_estimators=500,max_features="sqrt",
                class_weight="balanced",random_state=SEED+r,n_jobs=-1)
            m.fit(Xd.iloc[tr],y[tr]); pred[te]=m.predict_proba(Xd.iloc[te])[:,1]
        out.append(roc_auc_score(y,pred))
    return np.array(out)
def ic_boot(yt,ys,n=1000):
    v,idx=[],np.arange(len(yt))
    for _ in range(n):
        s=RNG.choice(idx,len(idx),replace=True)
        if len(np.unique(yt[s]))<2: continue
        v.append(roc_auc_score(yt[s],ys[s]))
    return np.percentile(v,[2.5,97.5])
log("\n"+"="*80); log("1. CLASIFICADOR DE CENTRO (ignorando diagnostico)"); log("="*80)
yc=(M.centro=="Kayseri").astype(int).values; a=cv_auc(XC,yc)
log(f"  Estambul vs Kayseri, todos los sujetos: AUC = {a.mean():.3f} +/- {a.std():.3f}")
log("\n"+"="*80); log("2. CONTROLES SANOS de una ciudad vs de la otra"); log("="*80)
mc=(M.grupo=="C").values; Xc2=XC[mc]; yc2=(M.centro[mc]=="Kayseri").astype(int).values
log(f"  n: Estambul={sum(yc2==0)}  Kayseri={sum(yc2==1)}")
a2=cv_auc(Xc2,yc2); log(f"  AUC control-control: {a2.mean():.3f} +/- {a2.std():.3f}")
log("  REFERENCIA China (estudios distintos): 0.998 - 1.000")
log("\n"+"="*80); log("3. DIAGNOSTICO (AD vs control) DENTRO de cada centro"); log("="*80)
log(f"{'Centro':<16}{'n':>6}{'AUC':>9}{'DE':>8}{'IC95':>20}")
mods={}
for c in ("Estambul","Kayseri"):
    m=((M.centro==c)&M.grupo.isin(["AD","C"])).values
    Xs=XC[m]; ys=(M.grupo[m]=="AD").astype(int).values; a3=cv_auc(Xs,ys)
    skf=StratifiedKFold(10,shuffle=True,random_state=SEED); pred=np.zeros(len(ys))
    for tr,te in skf.split(Xs,ys):
        md=RandomForestClassifier(n_estimators=500,max_features="sqrt",
            class_weight="balanced",random_state=SEED,n_jobs=-1)
        md.fit(Xs.iloc[tr],ys[tr]); pred[te]=md.predict_proba(Xs.iloc[te])[:,1]
    lo,hi=ic_boot(ys,pred)
    log(f"{c:<16}{len(ys):>6}{a3.mean():>9.3f}{a3.std():>8.3f}   [{lo:.3f}, {hi:.3f}]")
    mods[c]=(Xs,ys)
log("\n"+"="*80); log("4. TRANSFERENCIA ENTRE CENTROS (mismo estudio)"); log("="*80)
for a_,b_ in (("Estambul","Kayseri"),("Kayseri","Estambul")):
    Xa,ya=mods[a_]; Xb,yb=mods[b_]
    md=RandomForestClassifier(n_estimators=500,max_features="sqrt",
        class_weight="balanced",random_state=SEED,n_jobs=-1).fit(Xa,ya)
    sc=md.predict_proba(Xb[Xa.columns])[:,1]; auc=roc_auc_score(yb,sc); lo,hi=ic_boot(yb,sc)
    log(f"  train {a_:<10} -> test {b_:<10}  AUC = {auc:.3f}   [{lo:.3f}, {hi:.3f}]")
log("\n"+"="*80); log("5. CONFUSOR: profundidad por centro"); log("="*80)
prof=X.sum(axis=1).values
e=prof[(M.centro=="Estambul").values]; k=prof[(M.centro=="Kayseri").values]
_,p=stats.mannwhitneyu(e,k,alternative="two-sided")
log(f"  Estambul mediana={np.median(e):.0f}  Kayseri mediana={np.median(k):.0f}  p={p:.4f}")
log("  (las dos ciudades usaron quimicas MiSeq distintas: 284 vs 234 pb)")
log("\n"+"="*80); log("6. CONCORDANCIA DE BIOMARCADORES ENTRE CENTROS"); log("="*80)
try:
    import shap; usar=True
except ImportError: usar=False
imp={}
for c in ("Estambul","Kayseri"):
    Xs,ys=mods[c]; reps=[]
    for r in range(20):
        md=RandomForestClassifier(n_estimators=500,max_features="sqrt",
            class_weight="balanced",random_state=SEED+r,n_jobs=-1).fit(Xs,ys)
        if usar:
            sv=shap.TreeExplainer(md).shap_values(Xs)
            if isinstance(sv,list): sv=sv[1]
            if getattr(sv,"ndim",2)==3: sv=sv[:,:,1]
            reps.append(pd.Series(np.abs(sv).mean(axis=0),index=Xs.columns))
        else: reps.append(pd.Series(md.feature_importances_,index=Xs.columns))
    imp[c]=pd.concat(reps,axis=1).mean(axis=1)
rho,p=stats.spearmanr(imp["Estambul"],imp["Kayseri"])
log(f"  Spearman rankings SHAP: rho = {rho:.3f}  p = {p:.4f}")
log("  REFERENCIA China: -0.23 a 0.48, ninguno significativo")
t={c:set(imp[c].sort_values(ascending=False).head(20).index) for c in imp}
it=t["Estambul"]&t["Kayseri"]
log(f"\n  Top-20 compartidos: {len(it)}/20")
if it: log(f"    {', '.join(sorted(it))}")
log("  REFERENCIA China: 0/20 comun a las tres ciudades")
pd.DataFrame(imp).to_csv(f"{OUT}/shap_turquia_por_centro.csv")
M.to_csv(f"{OUT}/turquia_centros_mapa.csv")
log("\n"+"="*80)
log("LECTURA: si (1) y (2) dan AUC alto como en China -> separacion GEOGRAFICA")
log("         si dan AUC bajo -> lo de China era EFECTO DE ESTUDIO")
log("="*80)
