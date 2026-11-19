#!/usr/bin/env python3
"""Microbiome Radar - Cinco cohortes: China, Turquia, Italia"""
import os, json, warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
SEED=42; RNG=np.random.default_rng(42)
BASE=os.path.expanduser("~/Documents/MicrobiomeRadar_Research")
AD=os.path.join(BASE,"03_Alzheimer_Neurodegenerative")
OUT=os.path.join(BASE,"_crosscohort")
os.makedirs(os.path.join(OUT,"transfer_matrix"), exist_ok=True)
PSEUDO=1e-6; N_REP=20; N_PERM=300; N_BOOT=1000
COH = {
 "Shanghai": dict(m=f"{AD}/03d_Li2019_Shanghai_China/02_processed/Li_genus_counts.csv",
   p=f"{AD}/03d_Li2019_Shanghai_China/00_screening/mapa_heces.tsv",
   key="run_accession", col="grupo", modo="exacto", caso=["AD"], ctrl=["NC"], pais="China"),
 "Lishui": dict(m=f"{AD}/03g_Ling2021_Lishui_China/02_processed/Ling_genus_counts.csv",
   p=f"{AD}/03g_Ling2021_Lishui_China/00_screening/ena_PRJNA633959.tsv",
   key="run_accession", col="library_name", modo="prefijo", pais="China"),
 "Chongqing": dict(m=f"{AD}/03h_Zhuang2019_Chongqing_China/02_processed/Zhuang_genus_counts.csv",
   p=f"{AD}/03h_Zhuang2019_Chongqing_China/00_screening/meta_554111.tsv",
   key="run_accession", col="sample_alias", modo="prefijo", pais="China"),
 "Estambul": dict(m=f"{AD}/03i_Yildirim2022_Turkey/02_processed/Turquia_genus_counts.csv",
   p=f"{AD}/03i_Yildirim2022_Turkey/00_screening/ena_PRJNA734525.tsv",
   key="run_accession", col="library_name", modo="turquia", centro="Estambul", pais="Turquia"),
 "Kayseri": dict(m=f"{AD}/03i_Yildirim2022_Turkey/02_processed/Turquia_genus_counts.csv",
   p=f"{AD}/03i_Yildirim2022_Turkey/00_screening/ena_PRJNA734525.tsv",
   key="run_accession", col="library_name", modo="turquia", centro="Kayseri", pais="Turquia"),
 "Brescia": dict(m=f"{AD}/03k_Italia_PRJEB89799/02_processed/Italia_genus_counts.csv",
   p=f"{AD}/03k_Italia_PRJEB89799/00_screening/mapa_baseline.tsv",
   key="run", col="grupo", modo="exacto", caso=["AD"], ctrl=["Control"], pais="Italia"),
}
def log(*a): print(*a, flush=True)
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
    for run in X.index:
        if run not in mp.index: continue
        lab=etiqueta(cf, mp.loc[run, cf["col"]])
        if lab is None: continue
        keep.append(run); y.append(lab)
    return X.loc[keep], np.array(y)
def clr(df):
    rel=df.div(df.sum(axis=1),axis=0).replace(0,PSEUDO)
    lg=np.log(rel); return lg.sub(lg.mean(axis=1),axis=0)
def cv_auc(X,y,nrep=N_REP,k=10,ntree=500):
    out=[]
    for r in range(nrep):
        skf=StratifiedKFold(k,shuffle=True,random_state=SEED+r); pred=np.zeros(len(y))
        for tr,te in skf.split(X,y):
            m=RandomForestClassifier(n_estimators=ntree,max_features="sqrt",
                class_weight="balanced",random_state=SEED+r,n_jobs=-1)
            m.fit(X.iloc[tr],y[tr]); pred[te]=m.predict_proba(X.iloc[te])[:,1]
        out.append(roc_auc_score(y,pred))
    return np.array(out)
def perm_test(X,y,real,n=N_PERM):
    nul=[]
    for _ in range(n):
        yp=RNG.permutation(y)
        skf=StratifiedKFold(10,shuffle=True,random_state=SEED); pred=np.zeros(len(yp))
        for tr,te in skf.split(X,yp):
            m=RandomForestClassifier(n_estimators=100,max_features="sqrt",
                class_weight="balanced",random_state=SEED,n_jobs=-1)
            m.fit(X.iloc[tr],yp[tr]); pred[te]=m.predict_proba(X.iloc[te])[:,1]
        nul.append(roc_auc_score(yp,pred))
    nul=np.array(nul); return (np.sum(nul>=real)+1)/(n+1)
def ic_boot(yt,ys,n=N_BOOT):
    v,idx=[],np.arange(len(yt))
    for _ in range(n):
        s=RNG.choice(idx,len(idx),replace=True)
        if len(np.unique(yt[s]))<2: continue
        v.append(roc_auc_score(yt[s],ys[s]))
    return np.percentile(v,[2.5,97.5])
def cliff(a,b):
    a,b=np.asarray(a),np.asarray(b); n=len(a)*len(b)
    return (sum((x>b).sum() for x in a)-sum((x<b).sum() for x in a))/n if n else np.nan
log("="*84); log("CINCO COHORTES - tres paises, seis ciudades"); log("="*84)
crudo={}
for n in COH:
    try:
        X,y=cargar(n); crudo[n]=(X,y)
        log(f"  {n:<12} {COH[n]['pais']:<9} n={len(y):>4}  AD={y.sum():>3}  ctrl={(y==0).sum():>3}  gen={X.shape[1]}")
    except Exception as e: log(f"  !! {n}: {e}")
comunes=sorted(set.intersection(*[set(X.columns) for X,_ in crudo.values()]))
log(f"\nGeneros comunes: {len(comunes)}")
feats=[g for g in comunes if min((X[g]>0).mean() for X,_ in crudo.values())>=0.20]
log(f"Tras prevalencia >=20% en TODAS: {len(feats)}")
CLR={n:(clr(X[feats]),y) for n,(X,y) in crudo.items()}; noms=list(CLR)
log("\n"+"="*84); log(f"1. AUC INTERNO (10-fold x {N_REP})"); log("="*84)
log(f"{'Ciudad':<14}{'Pais':<10}{'AUC':>8}{'DE':>8}{'IC95':>20}{'p':>9}")
internos={}
for n,(Xc,y) in CLR.items():
    a=cv_auc(Xc,y); m,sd=a.mean(),a.std()
    skf=StratifiedKFold(10,shuffle=True,random_state=SEED); pred=np.zeros(len(y))
    for tr,te in skf.split(Xc,y):
        md=RandomForestClassifier(n_estimators=500,max_features="sqrt",
            class_weight="balanced",random_state=SEED,n_jobs=-1)
        md.fit(Xc.iloc[tr],y[tr]); pred[te]=md.predict_proba(Xc.iloc[te])[:,1]
    lo,hi=ic_boot(y,pred); p=perm_test(Xc,y,m)
    internos[n]=dict(auc=m,sd=sd,ic=[lo,hi],p=p)
    log(f"{n:<14}{COH[n]['pais']:<10}{m:>8.3f}{sd:>8.3f}    [{lo:.3f}, {hi:.3f}]{p:>9.4f}")
log("\n"+"="*84); log("2. MATRIZ DE TRANSFERENCIA"); log("="*84)
M=pd.DataFrame(index=noms,columns=noms,dtype=float); det=[]
for a in noms:
    Xa,ya=CLR[a]
    md=RandomForestClassifier(n_estimators=500,max_features="sqrt",
        class_weight="balanced",random_state=SEED,n_jobs=-1).fit(Xa,ya)
    for b in noms:
        if a==b: M.loc[a,b]=internos[a]["auc"]; continue
        Xb,yb=CLR[b]; sc=md.predict_proba(Xb[Xa.columns])[:,1]
        auc=roc_auc_score(yb,sc); lo,hi=ic_boot(yb,sc)
        M.loc[a,b]=auc
        det.append(dict(train=a,test=b,pais_train=COH[a]["pais"],pais_test=COH[b]["pais"],
                        auc=auc,ic_lo=lo,ic_hi=hi))
log("\n"+M.round(3).to_string())
M.to_csv(f"{OUT}/transfer_matrix/matriz_5cohortes.csv")
D=pd.DataFrame(det); D.to_csv(f"{OUT}/transfer_matrix/transferencia_5cohortes.csv",index=False)
intra=D[D.pais_train==D.pais_test]; inter=D[D.pais_train!=D.pais_test]
log(f"\n  interno (misma ciudad)  : {np.mean([M.loc[n,n] for n in noms]):.3f}")
log(f"  intra-pais, otra ciudad : {intra.auc.mean():.3f}   (n={len(intra)})")
log(f"  inter-pais              : {inter.auc.mean():.3f}   (n={len(inter)})")
log(f"  pares < 0.5             : {(D.auc<0.5).sum()}/{len(D)}")
log("\n"+"="*84); log("3. SEPARABILIDAD ENTRE COHORTES"); log("="*84)
log(f"{'Par':<30}{'todos':>9}{'solo controles':>17}")
for i,a in enumerate(noms):
    for b in noms[i+1:]:
        Xa,ya=CLR[a]; Xb,yb=CLR[b]
        R=pd.concat([Xa,Xb],ignore_index=True); yy=np.array([0]*len(Xa)+[1]*len(Xb))
        t=cv_auc(R,yy,nrep=5).mean()
        Ra=Xa[ya==0]; Rb=Xb[yb==0]
        Rc=pd.concat([Ra,Rb],ignore_index=True); yc=np.array([0]*len(Ra)+[1]*len(Rb))
        c=cv_auc(Rc,yc,nrep=5).mean()
        log(f"{a+' vs '+b:<30}{t:>9.3f}{c:>17.3f}")
log("\n"+"="*84); log("4. DIRECCION DE BIOMARCADORES CLAVE"); log("="*84)
foco=["Bifidobacterium","Bacteroides","Akkermansia","Faecalibacterium",
      "Prevotella_9","Escherichia-Shigella","Lactobacillus","Sutterella",
      "Blautia","Streptococcus","Roseburia","Alistipes"]
def resolver(g,cols):
    if g in cols: return g
    c=[x for x in cols if x==g or x.startswith(g+"_") or x.startswith(g+"-")]
    return c[0] if c else None
log(f"{'Genero':<22}" + "".join(f"{n[:9]:>11}" for n in noms))
td={}
for g in foco:
    fila=[]
    for n in noms:
        X,y=crudo[n]; col=resolver(g,X.columns)
        if col is None: fila.append("-"); continue
        rel=X.div(X.sum(axis=1),axis=0)
        a=rel.loc[y==1,col].values; c=rel.loc[y==0,col].values
        d=cliff(a,c)
        fila.append("plano" if abs(d)<0.1 else f"{'+' if d>0 else '-'}{abs(d):.2f}")
    td[g]=fila
    log(f"{g:<22}" + "".join(f"{v:>11}" for v in fila))
pd.DataFrame(td,index=noms).T.to_csv(f"{OUT}/direcciones_5cohortes.csv")
json.dump({k:{kk:(vv.tolist() if isinstance(vv,np.ndarray) else vv)
    for kk,vv in v.items()} for k,v in internos.items()},
    open(f"{OUT}/auc_5cohortes.json","w"), indent=2, default=float)
log("\n"+"="*84); log("LISTO. Salidas en "+OUT); log("="*84)
