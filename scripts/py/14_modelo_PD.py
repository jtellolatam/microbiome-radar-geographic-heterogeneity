#!/usr/bin/env python3
"""Microbiome Radar - Modelado brazo Parkinson: Finlandia y Rumania"""
import os, json, warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
SEED=42; RNG=np.random.default_rng(42); PSEUDO=1e-6
N_REP=20; N_PERM=300; N_BOOT=1000
BASE=os.path.expanduser("~/Documents/MicrobiomeRadar_Research")
PD=os.path.join(BASE,"02_Parkinson_Neurodegenerative")
OUT=os.path.join(BASE,"_crosscohort")
os.makedirs(os.path.join(OUT,"transfer_matrix"), exist_ok=True)
COH = {
 "Finlandia": dict(m=f"{PD}/02e_Finlandia/02_processed/Finlandia_genus_counts.csv",
   p=f"{PD}/02e_Finlandia/00_screening/mapa_finlandia.tsv",
   key="run", col="grupo", caso="PD", ctrl="Control", pais="Finlandia"),
 "Rumania": dict(m=f"{PD}/02f_Rumania/02_processed/Rumania_genus_counts.csv",
   p=f"{PD}/02f_Rumania/00_screening/mapa_rumania.tsv",
   key="muestra", col="grupo", caso="PD", ctrl="Control", pais="Rumania"),
}
def log(*a): print(*a, flush=True)
def cargar(n):
    cf=COH[n]; X=pd.read_csv(cf["m"], index_col=0)
    mp=pd.read_csv(cf["p"], sep="\t", dtype=str)
    mp=mp[~mp[cf["key"]].duplicated(keep="first")].set_index(cf["key"])
    y,keep=[],[]
    for s in X.index:
        if s not in mp.index: continue
        g=mp.loc[s, cf["col"]]
        if g==cf["caso"]: y.append(1); keep.append(s)
        elif g==cf["ctrl"]: y.append(0); keep.append(s)
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
def resolver(g,cols):
    if g in cols: return g
    c=[x for x in cols if x==g or x.startswith(g+"_") or x.startswith(g+"-")]
    return c[0] if c else None
log("="*84); log("BRAZO PARKINSON - dos cohortes europeas"); log("="*84)
crudo={}
for n in COH:
    X,y=cargar(n); crudo[n]=(X,y)
    log(f"  {n:<11} {COH[n]['pais']:<10} n={len(y):>4}  PD={y.sum():>3}  ctrl={(y==0).sum():>3}  gen={X.shape[1]}")
comunes=sorted(set.intersection(*[set(X.columns) for X,_ in crudo.values()]))
feats=[g for g in comunes if min((X[g]>0).mean() for X,_ in crudo.values())>=0.20]
log(f"\nGeneros comunes: {len(comunes)}   tras prevalencia >=20%: {len(feats)}")
CLR={n:(clr(X[feats]),y) for n,(X,y) in crudo.items()}
noms=list(CLR)
log("\n"+"="*84); log(f"1. AUC INTERNO (10-fold x {N_REP})"); log("="*84)
log(f"{'Cohorte':<13}{'AUC':>8}{'DE':>8}{'IC95':>20}{'p':>9}")
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
    log(f"{n:<13}{m:>8.3f}{sd:>8.3f}    [{lo:.3f}, {hi:.3f}]{p:>9.4f}")
log("\n"+"="*84); log("2. TRANSFERENCIA"); log("="*84)
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
        det.append(dict(train=a,test=b,auc=auc,ic_lo=lo,ic_hi=hi))
log("\n"+M.round(3).to_string())
for d in det:
    log(f"  {d['train']:>10} -> {d['test']:<10} AUC={d['auc']:.3f}  [{d['ic_lo']:.3f}, {d['ic_hi']:.3f}]")
M.to_csv(f"{OUT}/transfer_matrix/matriz_PD.csv")
pd.DataFrame(det).to_csv(f"{OUT}/transfer_matrix/transferencia_PD.csv",index=False)
log("\n"+"="*84); log("3. CONTROL DE LOTE"); log("="*84)
Xa,ya=CLR["Finlandia"]; Xb,yb=CLR["Rumania"]
R=pd.concat([Xa,Xb],ignore_index=True); yy=np.array([0]*len(Xa)+[1]*len(Xb))
log(f"  Clasificar cohorte, todos los sujetos : {cv_auc(R,yy,nrep=5).mean():.3f}")
Ra=Xa[ya==0]; Rb=Xb[yb==0]
Rc=pd.concat([Ra,Rb],ignore_index=True); yc=np.array([0]*len(Ra)+[1]*len(Rb))
log(f"  Solo controles sanos                  : {cv_auc(Rc,yc,nrep=5).mean():.3f}")
log("  REFERENCIA: Alzheimer 0.87-1.00   TEA 0.99-1.00")
log("\n"+"="*84); log("4. SHAP"); log("="*84)
try:
    import shap; usar=True
except ImportError: log("  shap no instalado"); usar=False
imp_rep={n:[] for n in noms}
for n,(Xc,y) in CLR.items():
    for r in range(20):
        md=RandomForestClassifier(n_estimators=500,max_features="sqrt",
            class_weight="balanced",random_state=SEED+r,n_jobs=-1).fit(Xc,y)
        if usar:
            sv=shap.TreeExplainer(md).shap_values(Xc)
            if isinstance(sv,list): sv=sv[1]
            if getattr(sv,"ndim",2)==3: sv=sv[:,:,1]
            imp_rep[n].append(pd.Series(np.abs(sv).mean(axis=0), index=Xc.columns))
        else:
            imp_rep[n].append(pd.Series(md.feature_importances_, index=Xc.columns))
log("\nEstabilidad INTRA-cohorte:")
for n in noms:
    rs=[stats.spearmanr(imp_rep[n][i],imp_rep[n][j])[0] for i in range(20) for j in range(i+1,20)]
    log(f"  {n:<13} rho medio = {np.mean(rs):.3f}")
imp={n: pd.concat(imp_rep[n],axis=1).mean(axis=1) for n in noms}
rho,p=stats.spearmanr(imp["Finlandia"],imp["Rumania"])
log(f"\nConcordancia Finlandia vs Rumania: rho = {rho:.3f}  p = {p:.4f}")
tops={n:set(imp[n].sort_values(ascending=False).head(20).index) for n in noms}
it=tops["Finlandia"]&tops["Rumania"]
log(f"Top-20 compartidos: {len(it)}/20")
if it: log(f"  {', '.join(sorted(it))}")
for n in noms:
    log(f"\n--- {n} ---")
    log(imp[n].sort_values(ascending=False).head(20).round(5).to_string())
pd.DataFrame(imp).to_csv(f"{OUT}/shap_PD_por_cohorte.csv")
log("\n"+"="*84); log("5. DIRECCION DE BIOMARCADORES"); log("="*84)
cont={}
for n in COH:
    for g in crudo[n][0].columns:
        base=g.split("_")[0] if not g.startswith("[") else g
        cont[base]=cont.get(base,0)+1
foco=sorted([g for g,c in cont.items() if c>=2])
res=[]
for g in foco:
    fila={}; ok=True
    for n in noms:
        X,y=crudo[n]; rel=X.div(X.sum(axis=1),axis=0)
        col=resolver(g,rel.columns)
        if col is None: ok=False; break
        a=rel.loc[y==1,col].values; c=rel.loc[y==0,col].values
        if a.mean()+c.mean()<1e-5: ok=False; break
        fila[f"d_{n}"]=cliff(a,c)
        _,pv=stats.mannwhitneyu(a,c,alternative="two-sided"); fila[f"p_{n}"]=pv
    if not ok: continue
    df_,dr_=fila["d_Finlandia"],fila["d_Rumania"]
    fila["genero"]=g
    fila["patron"]=("INVIERTE" if np.sign(df_)!=np.sign(dr_) else "consistente") \
        if (abs(df_)>=0.1 and abs(dr_)>=0.1) else \
        ("solo Finlandia" if abs(df_)>=0.1 else "solo Rumania" if abs(dr_)>=0.1 else "sin efecto")
    fila["transversal"]=(np.sign(df_)==np.sign(dr_) and abs(df_)>=0.1 and abs(dr_)>=0.1
                         and fila["p_Finlandia"]<0.05 and fila["p_Rumania"]<0.05)
    res.append(fila)
R2=pd.DataFrame(res)
R2["_a"]=R2[["d_Finlandia","d_Rumania"]].abs().max(axis=1)
R2=R2.sort_values("_a",ascending=False).drop(columns="_a")
cols=["genero","d_Finlandia","p_Finlandia","d_Rumania","p_Rumania","patron"]
log(f"\nGeneros evaluados: {len(R2)}")
log(R2.patron.value_counts().to_string())
log("\n>>> TRANSVERSALES:")
TR=R2[R2.transversal]
log(TR[cols].round(4).to_string(index=False) if len(TR) else "   ninguno")
log("\n>>> INVIERTEN:")
INV=R2[R2.patron=="INVIERTE"]
log(INV[cols].round(4).head(15).to_string(index=False) if len(INV) else "   ninguno")
log("\n>>> TOP 20 por magnitud:")
log(R2[cols].head(20).round(4).to_string(index=False))
R2.to_csv(f"{OUT}/direcciones_PD.csv", index=False)
json.dump({k:{kk:(vv.tolist() if isinstance(vv,np.ndarray) else vv)
    for kk,vv in v.items()} for k,v in internos.items()},
    open(f"{OUT}/auc_PD.json","w"), indent=2, default=float)
log("\n"+"="*84); log("LISTO"); log("="*84)
