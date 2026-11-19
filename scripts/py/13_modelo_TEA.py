#!/usr/bin/env python3
"""Microbiome Radar - Modelado brazo TEA: Dan, Xiamen, Caltech"""
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
TEA=os.path.join(BASE,"01_TEA_Neurodivergent")
OUT=os.path.join(BASE,"_crosscohort")
os.makedirs(os.path.join(OUT,"transfer_matrix"), exist_ok=True)
COH = {
 "Nanjing": dict(m=f"{TEA}/01f_Dan_Nanjing/02_processed/Dan_genus_counts.csv",
   p=f"{TEA}/01f_Dan_Nanjing/00_screening/mapa_dan.tsv", key="run", col="grupo", pais="China"),
 "Xiamen": dict(m=f"{TEA}/01g_Xiamen/02_processed/Xiamen_genus_counts.csv",
   p=f"{TEA}/01g_Xiamen/00_screening/mapa_xiamen.tsv", key="run", col="grupo", pais="China"),
 "Caltech": dict(m=f"{TEA}/01h_Caltech_USA/02_processed/Caltech_genus_counts.csv",
   p=f"{TEA}/01h_Caltech_USA/00_screening/mapa_caltech.tsv", key="run", col="grupo", pais="EEUU"),
}
def log(*a): print(*a, flush=True)
def cargar(n):
    cf=COH[n]; X=pd.read_csv(cf["m"], index_col=0)
    mp=pd.read_csv(cf["p"], sep="\t", dtype=str)
    mp=mp[~mp[cf["key"]].duplicated(keep="first")].set_index(cf["key"])
    y,keep=[],[]
    for run in X.index:
        if run not in mp.index: continue
        g=mp.loc[run, cf["col"]]
        if g=="ASD": y.append(1); keep.append(run)
        elif g=="Control": y.append(0); keep.append(run)
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
def cliff_ic(a,b,nboot=500):
    v=[]; ia,ib=np.arange(len(a)),np.arange(len(b))
    for _ in range(nboot):
        v.append(cliff(a[RNG.choice(ia,len(a),True)], b[RNG.choice(ib,len(b),True)]))
    return np.percentile(v,[2.5,97.5])
def resolver(g,cols):
    if g in cols: return g
    c=[x for x in cols if x==g or x.startswith(g+"_") or x.startswith(g+"-")]
    return c[0] if c else None
log("="*88); log("BRAZO TEA - tres cohortes, dos paises"); log("="*88)
crudo={}
for n in COH:
    X,y=cargar(n); crudo[n]=(X,y)
    log(f"  {n:<10} {COH[n]['pais']:<7} n={len(y):>4}  ASD={y.sum():>3}  ctrl={(y==0).sum():>3}  gen={X.shape[1]}")
comunes=sorted(set.intersection(*[set(X.columns) for X,_ in crudo.values()]))
feats=[g for g in comunes if min((X[g]>0).mean() for X,_ in crudo.values())>=0.20]
log(f"\nGeneros comunes: {len(comunes)}   tras prevalencia >=20%: {len(feats)}")
CLR={n:(clr(X[feats]),y) for n,(X,y) in crudo.items()}
noms=list(CLR)
log("\n"+"="*88); log(f"1. AUC INTERNO (10-fold x {N_REP})"); log("="*88)
log(f"{'Cohorte':<12}{'Pais':<8}{'AUC':>8}{'DE':>8}{'IC95':>20}{'p':>9}")
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
    log(f"{n:<12}{COH[n]['pais']:<8}{m:>8.3f}{sd:>8.3f}    [{lo:.3f}, {hi:.3f}]{p:>9.4f}")
log("\n"+"="*88); log("2. MATRIZ DE TRANSFERENCIA"); log("="*88)
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
M.to_csv(f"{OUT}/transfer_matrix/matriz_TEA.csv")
D=pd.DataFrame(det); D.to_csv(f"{OUT}/transfer_matrix/transferencia_TEA.csv",index=False)
intra=D[D.pais_train==D.pais_test]; inter=D[D.pais_train!=D.pais_test]
log(f"\n  interno medio      : {np.mean([M.loc[n,n] for n in noms]):.3f}")
if len(intra): log(f"  intra-pais (China) : {intra.auc.mean():.3f}   (n={len(intra)})")
log(f"  inter-pais         : {inter.auc.mean():.3f}   (n={len(inter)})")
log(f"  pares < 0.5        : {(D.auc<0.5).sum()}/{len(D)}")
log("\n"+"="*88); log("3. CONTROL DE LOTE"); log("="*88)
log(f"{'Par':<26}{'todos':>9}{'solo controles':>17}")
for i,a in enumerate(noms):
    for b in noms[i+1:]:
        Xa,ya=CLR[a]; Xb,yb=CLR[b]
        R=pd.concat([Xa,Xb],ignore_index=True); yy=np.array([0]*len(Xa)+[1]*len(Xb))
        t=cv_auc(R,yy,nrep=5).mean()
        Ra=Xa[ya==0]; Rb=Xb[yb==0]
        Rc=pd.concat([Ra,Rb],ignore_index=True); yc=np.array([0]*len(Ra)+[1]*len(Rb))
        c=cv_auc(Rc,yc,nrep=5).mean()
        log(f"{a+' vs '+b:<26}{t:>9.3f}{c:>17.3f}")
log("\n  REFERENCIA Alzheimer: cohorte 0.93-0.99   controles 0.87-1.00")
log("\n"+"="*88); log("4. SHAP"); log("="*88)
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
    log(f"  {n:<12} rho medio = {np.mean(rs):.3f}")
imp={n: pd.concat(imp_rep[n],axis=1).mean(axis=1) for n in noms}
log("\nConcordancia ENTRE cohortes:")
log(f"{'Par':<28}{'Spearman':>10}{'p':>10}")
for i,a in enumerate(noms):
    for b in noms[i+1:]:
        rho,p=stats.spearmanr(imp[a],imp[b])
        log(f"{a+' vs '+b:<28}{rho:>10.3f}{p:>10.4f}")
tops={n:set(imp[n].sort_values(ascending=False).head(20).index) for n in noms}
log("\nSolapamiento top-20:")
for i,a in enumerate(noms):
    for b in noms[i+1:]:
        it=tops[a]&tops[b]
        log(f"  {a} & {b}: {len(it)}/20  {', '.join(sorted(it)) if it else ''}")
cm=set.intersection(*tops.values())
log(f"  Comun a las TRES: {len(cm)}/20  {', '.join(sorted(cm)) if cm else ''}")
pd.DataFrame(imp).to_csv(f"{OUT}/shap_TEA_por_cohorte.csv")
log("\n"+"="*88); log("5. CONTRASTE AGRUPADO CHINA vs EEUU"); log("="*88)
chi=[n for n in COH if COH[n]["pais"]=="China"]
usa=[n for n in COH if COH[n]["pais"]=="EEUU"]
log(f"  China: {sum(len(crudo[n][1]) for n in chi)}   EEUU: {sum(len(crudo[n][1]) for n in usa)}\n")
cont={}
for n in COH:
    for g in crudo[n][0].columns:
        base=g.split("_")[0] if not g.startswith("[") else g
        cont[base]=cont.get(base,0)+1
foco=sorted([g for g,c in cont.items() if c>=2])
def agrupar(nombres,g):
    A,C=[],[]
    for n in nombres:
        X,y=crudo[n]; rel=X.div(X.sum(axis=1),axis=0)
        col=resolver(g,rel.columns)
        if col is None: continue
        A.append(rel.loc[y==1,col].values); C.append(rel.loc[y==0,col].values)
    if not A: return None,None
    return np.concatenate(A), np.concatenate(C)
res=[]
for g in foco:
    ac,cc=agrupar(chi,g); au,cu=agrupar(usa,g)
    if ac is None or au is None: continue
    if len(ac)<20 or len(au)<20: continue
    if (ac.mean()+cc.mean())<1e-5 or (au.mean()+cu.mean())<1e-5: continue
    dc=cliff(ac,cc); du=cliff(au,cu)
    lo_c,hi_c=cliff_ic(ac,cc); lo_u,hi_u=cliff_ic(au,cu)
    sol = not (hi_c<lo_u or hi_u<lo_c)
    _,pc=stats.mannwhitneyu(ac,cc,alternative="two-sided")
    _,pu=stats.mannwhitneyu(au,cu,alternative="two-sided")
    res.append(dict(genero=g, d_china=dc, ic_china=f"[{lo_c:.2f},{hi_c:.2f}]", p_china=pc,
        d_eeuu=du, ic_eeuu=f"[{lo_u:.2f},{hi_u:.2f}]", p_eeuu=pu, solapan=sol,
        inversion=(np.sign(dc)!=np.sign(du) and abs(dc)>=0.1 and abs(du)>=0.1 and not sol),
        transversal=(np.sign(dc)==np.sign(du) and abs(dc)>=0.1 and abs(du)>=0.1 and pc<0.05 and pu<0.05)))
R=pd.DataFrame(res)
R["_a"]=R[["d_china","d_eeuu"]].abs().max(axis=1)
R=R.sort_values("_a",ascending=False).drop(columns="_a")
log(">>> INVERSIONES CONFIRMADAS:")
INV=R[R.inversion]
log(INV[["genero","d_china","ic_china","p_china","d_eeuu","ic_eeuu","p_eeuu"]].round(4).to_string(index=False) if len(INV) else "   ninguna")
log("\n>>> TRANSVERSALES:")
TR=R[R.transversal]
log(TR[["genero","d_china","p_china","d_eeuu","p_eeuu"]].round(4).to_string(index=False) if len(TR) else "   ninguno")
log("\n>>> TOP 20 por magnitud:")
log(R[["genero","d_china","p_china","d_eeuu","p_eeuu","solapan","inversion"]].head(20).round(4).to_string(index=False))
R.to_csv(f"{OUT}/contraste_TEA_china_eeuu.csv", index=False)
json.dump({k:{kk:(vv.tolist() if isinstance(vv,np.ndarray) else vv)
    for kk,vv in v.items()} for k,v in internos.items()},
    open(f"{OUT}/auc_TEA.json","w"), indent=2, default=float)
log(f"\nGeneros: {len(R)}   inversiones: {R.inversion.sum()}   transversales: {R.transversal.sum()}")
log("="*88)
