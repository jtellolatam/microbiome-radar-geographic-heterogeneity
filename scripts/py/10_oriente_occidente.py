#!/usr/bin/env python3
"""Contraste ORIENTE (China) vs OCCIDENTE (Turquia + Italia)"""
import os, warnings
import numpy as np
import pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
BASE=os.path.expanduser("~/Documents/MicrobiomeRadar_Research")
AD=os.path.join(BASE,"03_Alzheimer_Neurodegenerative")
OUT=os.path.join(BASE,"_crosscohort")
RNG=np.random.default_rng(42); PSEUDO=1e-6
COH = {
 "Shanghai": dict(m=f"{AD}/03d_Li2019_Shanghai_China/02_processed/Li_genus_counts.csv",
   p=f"{AD}/03d_Li2019_Shanghai_China/00_screening/mapa_heces.tsv",
   key="run_accession", col="grupo", modo="exacto", caso=["AD"], ctrl=["NC"],
   region="Oriente", pais="China"),
 "Lishui": dict(m=f"{AD}/03g_Ling2021_Lishui_China/02_processed/Ling_genus_counts.csv",
   p=f"{AD}/03g_Ling2021_Lishui_China/00_screening/ena_PRJNA633959.tsv",
   key="run_accession", col="library_name", modo="prefijo", region="Oriente", pais="China"),
 "Chongqing": dict(m=f"{AD}/03h_Zhuang2019_Chongqing_China/02_processed/Zhuang_genus_counts.csv",
   p=f"{AD}/03h_Zhuang2019_Chongqing_China/00_screening/meta_554111.tsv",
   key="run_accession", col="sample_alias", modo="prefijo", region="Oriente", pais="China"),
 "Estambul": dict(m=f"{AD}/03i_Yildirim2022_Turkey/02_processed/Turquia_genus_counts.csv",
   p=f"{AD}/03i_Yildirim2022_Turkey/00_screening/ena_PRJNA734525.tsv",
   key="run_accession", col="library_name", modo="turquia", centro="Estambul",
   region="Occidente", pais="Turquia"),
 "Kayseri": dict(m=f"{AD}/03i_Yildirim2022_Turkey/02_processed/Turquia_genus_counts.csv",
   p=f"{AD}/03i_Yildirim2022_Turkey/00_screening/ena_PRJNA734525.tsv",
   key="run_accession", col="library_name", modo="turquia", centro="Kayseri",
   region="Occidente", pais="Turquia"),
 "Brescia": dict(m=f"{AD}/03k_Italia_PRJEB89799/02_processed/Italia_genus_counts.csv",
   p=f"{AD}/03k_Italia_PRJEB89799/00_screening/mapa_baseline.tsv",
   key="run", col="grupo", modo="exacto", caso=["AD"], ctrl=["Control"],
   region="Occidente", pais="Italia"),
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
    X=X.loc[keep]
    return X.div(X.sum(axis=1),axis=0), np.array(y)
def cliff(a,b):
    a,b=np.asarray(a),np.asarray(b); n=len(a)*len(b)
    if n==0: return np.nan
    return (sum((x>b).sum() for x in a)-sum((x<b).sum() for x in a))/n
def cliff_ic(a,b,nboot=500):
    v=[]; ia,ib=np.arange(len(a)),np.arange(len(b))
    for _ in range(nboot):
        v.append(cliff(a[RNG.choice(ia,len(a),True)], b[RNG.choice(ib,len(b),True)]))
    return np.percentile(v,[2.5,97.5])
def resolver(g,cols):
    if g in cols: return g
    c=[x for x in cols if x==g or x.startswith(g+"_") or x.startswith(g+"-")]
    return c[0] if c else None
log("="*96); log("CONTRASTE ORIENTE vs OCCIDENTE"); log("="*96)
D={}
for n in COH:
    rel,y=cargar(n); D[n]=(rel,y)
    log(f"  {n:<11} {COH[n]['region']:<10} {COH[n]['pais']:<9} n={len(y):>4}  AD={y.sum():>3}  ctrl={(y==0).sum():>3}  gen={rel.shape[1]}")
ori=[n for n in COH if COH[n]["region"]=="Oriente"]
occ=[n for n in COH if COH[n]["region"]=="Occidente"]
log(f"\n  ORIENTE: {sum(len(D[n][1]) for n in ori)} muestras   OCCIDENTE: {sum(len(D[n][1]) for n in occ)} muestras")
cont={}
for n in COH:
    for g in D[n][0].columns:
        base=g.split("_")[0] if not g.startswith("[") else g
        cont[base]=cont.get(base,0)+1
foco=sorted([g for g,c in cont.items() if c>=4])
log("\n"+"="*96); log("A. DIRECCION POR COHORTE (cada una con TODOS sus generos)"); log("="*96)
log(f"Generos presentes en >=4 cohortes: {len(foco)}")
filas=[]
for g in foco:
    fila={"genero":g}; ds={}
    for n in COH:
        rel,y=D[n]; col=resolver(g,rel.columns)
        if col is None: fila[n]=np.nan; continue
        a=rel.loc[y==1,col].values; c=rel.loc[y==0,col].values
        if a.mean()+c.mean()<1e-5: fila[n]=np.nan; continue
        d=cliff(a,c); fila[n]=d; ds[n]=d
    do=[ds[n] for n in ori if n in ds and abs(ds[n])>=0.10]
    dw=[ds[n] for n in occ if n in ds and abs(ds[n])>=0.10]
    fila["n_ori"]=len(do); fila["n_occ"]=len(dw)
    if do and dw:
        so={np.sign(x) for x in do}; sw={np.sign(x) for x in dw}
        fila["patron"]=("INVIERTE" if so!=sw else "consistente") if (len(so)==1 and len(sw)==1) else "mixto"
    elif len(do)>=2 and not dw: fila["patron"]="solo Oriente"
    elif len(dw)>=2 and not do: fila["patron"]="solo Occidente"
    else: fila["patron"]="sin efecto"
    filas.append(fila)
T=pd.DataFrame(filas).set_index("genero")
inv=T[T.patron=="INVIERTE"]
log(f"\n>>> INVIERTEN entre Oriente y Occidente: {len(inv)}")
if len(inv): log(T.loc[inv.index, list(COH)+["patron"]].round(2).to_string())
cons=T[T.patron=="consistente"]
log(f"\n>>> CONSISTENTES en ambas regiones: {len(cons)}")
if len(cons): log(T.loc[cons.index, list(COH)+["patron"]].round(2).to_string())
log(f"\n  solo Oriente: {(T.patron=='solo Oriente').sum()}   solo Occidente: {(T.patron=='solo Occidente').sum()}   mixto: {(T.patron=='mixto').sum()}   sin efecto: {(T.patron=='sin efecto').sum()}")
T.round(3).to_csv(f"{OUT}/direcciones_oriente_occidente.csv")
log("\n"+"="*96); log("B. CONTRASTE AGRUPADO - interaccion region x taxon"); log("="*96)
def agrupar(nombres,g):
    A,C=[],[]
    for n in nombres:
        rel,y=D[n]; col=resolver(g,rel.columns)
        if col is None: continue
        A.append(rel.loc[y==1,col].values); C.append(rel.loc[y==0,col].values)
    if not A: return None,None
    return np.concatenate(A), np.concatenate(C)
res=[]
for g in foco:
    ao,co=agrupar(ori,g); aw,cw=agrupar(occ,g)
    if ao is None or aw is None: continue
    if len(ao)<20 or len(aw)<20: continue
    if (ao.mean()+co.mean())<1e-5 or (aw.mean()+cw.mean())<1e-5: continue
    do=cliff(ao,co); dw=cliff(aw,cw)
    lo_o,hi_o=cliff_ic(ao,co); lo_w,hi_w=cliff_ic(aw,cw)
    solapan = not (hi_o<lo_w or hi_w<lo_o)
    _,po=stats.mannwhitneyu(ao,co,alternative="two-sided")
    _,pw=stats.mannwhitneyu(aw,cw,alternative="two-sided")
    res.append(dict(genero=g, delta_ori=do, ic_ori=f"[{lo_o:.2f},{hi_o:.2f}]", p_ori=po,
        delta_occ=dw, ic_occ=f"[{lo_w:.2f},{hi_w:.2f}]", p_occ=pw, solapan_IC=solapan,
        inversion=(np.sign(do)!=np.sign(dw) and abs(do)>=0.1 and abs(dw)>=0.1 and not solapan)))
R=pd.DataFrame(res)
R["_a"]=R[["delta_ori","delta_occ"]].abs().max(axis=1)
R=R.sort_values("_a",ascending=False).drop(columns="_a")
log("\n>>> INVERSIONES CONFIRMADAS (signo opuesto + IC no solapados):")
INV=R[R.inversion]
if len(INV): log(INV[["genero","delta_ori","ic_ori","p_ori","delta_occ","ic_occ","p_occ"]].round(4).to_string(index=False))
else: log("   ninguna supera el criterio")
log("\n>>> TOP 25 por magnitud:")
log(R[["genero","delta_ori","p_ori","delta_occ","p_occ","solapan_IC","inversion"]].head(25).round(4).to_string(index=False))
R.to_csv(f"{OUT}/contraste_oriente_occidente.csv", index=False)
log("\n"+"="*96)
log(f"Generos evaluados: {len(R)}   Inversiones confirmadas: {R.inversion.sum()}")
log("="*96)
