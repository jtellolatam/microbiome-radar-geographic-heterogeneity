#!/usr/bin/env python3
"""Microbiome Radar - Auditoria de direccion de biomarcadores (AD, China)
N1 = tabla de literatura | N2 = reportado por autores | N3 = nuestro pipeline"""
import os, warnings
import numpy as np
import pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")

BASE = os.path.expanduser("~/Documents/MicrobiomeRadar_Research")
AD = os.path.join(BASE, "03_Alzheimer_Neurodegenerative")
OUT = os.path.join(BASE, "_crosscohort")
os.makedirs(OUT, exist_ok=True)

N1_ALZHEIMER = {"Bacteroides":"ALTO","Bifidobacterium":"BAJO",
                "Prevotella":"BAJO","Escherichia":"ALTO"}
GENEROS_TFM = ["Faecalibacterium","Bifidobacterium","Akkermansia","Lactobacillus",
               "Blautia","Clostridium","Sutterella","Bacteroides","Veillonella","Prevotella"]

N2 = {
 "Li_Shanghai": {"Bifidobacterium":"ALTO","Escherichia":"ALTO","Akkermansia":"ALTO",
                 "Bacteroides":"BAJO","Streptococcus":"ALTO","Sutterella":"ALTO",
                 "Blautia":"ALTO"},
 "Ling_Lishui": {"Bifidobacterium":"ALTO","Akkermansia":"ALTO","Faecalibacterium":"BAJO",
                 "Roseburia":"BAJO","Lactobacillus":"ALTO","Enterococcus":"ALTO"},
 "Zhuang_Chongqing": {"Bifidobacterium":"ALTO","Bacteroides":"BAJO"},
}

COHORTES = {
 "Li_Shanghai": dict(
   matriz=f"{AD}/03d_Li2019_Shanghai_China/02_processed/Li_genus_counts.csv",
   mapa=f"{AD}/03d_Li2019_Shanghai_China/00_screening/mapa_heces.tsv",
   key="run_accession", grpcol="grupo", modo="exacto",
   caso=["AD"], ctrl=["NC"], ciudad="Shanghai"),
 "Ling_Lishui": dict(
   matriz=f"{AD}/03g_Ling2021_Lishui_China/02_processed/Ling_genus_counts.csv",
   mapa=f"{AD}/03g_Ling2021_Lishui_China/00_screening/ena_PRJNA633959.tsv",
   key="run_accession", grpcol="library_name", modo="prefijo", ciudad="Lishui"),
 "Zhuang_Chongqing": dict(
   matriz=f"{AD}/03h_Zhuang2019_Chongqing_China/02_processed/Zhuang_genus_counts.csv",
   mapa=f"{AD}/03h_Zhuang2019_Chongqing_China/00_screening/meta_554111.tsv",
   key="run_accession", grpcol="sample_alias", modo="prefijo", ciudad="Chongqing"),
}

def log(*a): print(*a, flush=True)

def etiqueta(cf, v):
    s = str(v).upper()
    if cf["modo"] == "exacto":
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
    X = X.loc[keep]
    return X.div(X.sum(axis=1), axis=0), np.array(y)

def resolver(g, cols):
    if g in cols: return g
    c = [x for x in cols if x == g or x.startswith(g+"_") or x.startswith(g+"-")]
    return c[0] if c else None

def cliff(a, b):
    a, b = np.asarray(a), np.asarray(b)
    n = len(a)*len(b)
    if n == 0: return np.nan
    return (sum((x > b).sum() for x in a) - sum((x < b).sum() for x in a)) / n

log("="*92)
log("AUDITORIA DE DIRECCION DE BIOMARCADORES - Alzheimer, cohortes chinas")
log("="*92)
datos = {}
for n in COHORTES:
    try:
        rel, y = cargar(n); datos[n] = (rel, y)
        log(f"  {n:<20} n={len(y):>4}  AD={y.sum():>3}  ctrl={(y==0).sum():>3}  gen={rel.shape[1]}")
    except Exception as e:
        log(f"  !! {n}: {e}")

generos = sorted(set(GENEROS_TFM) | set(N1_ALZHEIMER) | {g for d in N2.values() for g in d})
ciudades = [COHORTES[n]["ciudad"] for n in datos]

filas = []
for g in generos:
    f = {"genero": g, "N1_literatura": N1_ALZHEIMER.get(g, "-")}
    for nom, (rel, y) in datos.items():
        cd = COHORTES[nom]["ciudad"]
        col = resolver(g, rel.columns)
        if col is None:
            f[f"{cd}_N3"]="ausente"; f[f"{cd}_p"]=np.nan
            f[f"{cd}_delta"]=np.nan; f[f"{cd}_ratio"]=np.nan
        else:
            a = rel.loc[y==1, col].values; c = rel.loc[y==0, col].values
            ma, mc = a.mean(), c.mean()
            try: _, p = stats.mannwhitneyu(a, c, alternative="two-sided")
            except ValueError: p = np.nan
            d = cliff(a, c)
            f[f"{cd}_N3"] = "plano" if abs(d) < 0.1 else ("ALTO" if ma > mc else "BAJO")
            f[f"{cd}_p"]=p; f[f"{cd}_delta"]=d
            f[f"{cd}_ratio"] = ma/mc if mc>0 else np.inf
        f[f"{cd}_N2"] = N2.get(nom, {}).get(g, "")
    filas.append(f)
T = pd.DataFrame(filas).set_index("genero")

log("\n"+"="*92)
log("a) VALIDEZ DEL PIPELINE - reproducimos lo publicado por los autores?")
log("="*92)
log(f"{'Cohorte':<14}{'Genero':<20}{'N2 autor':<11}{'N3 nuestro':<12}{'p':>9}{'delta':>8}   veredicto")
ok = mal = 0
for nom in datos:
    cd = COHORTES[nom]["ciudad"]
    for g, rep in N2.get(nom, {}).items():
        if not rep or g not in T.index: continue
        nu = T.loc[g, f"{cd}_N3"]; p = T.loc[g, f"{cd}_p"]; d = T.loc[g, f"{cd}_delta"]
        if nu in ("ausente","plano"): v = "sin efecto"
        elif nu == rep: v = "COINCIDE"; ok += 1
        else: v = ">>> DISCREPA <<<"; mal += 1
        log(f"{cd:<14}{g:<20}{rep:<11}{str(nu):<12}{p:>9.4f}{d:>8.2f}   {v}")
log(f"\n  Coinciden: {ok}   Discrepan: {mal}")
log("  !! Pipeline NO reproduce lo publicado" if mal > ok else "  Pipeline consistente")

log("\n"+"="*92)
log("b) LITERATURA OCCIDENTAL vs COHORTES CHINAS")
log("="*92)
log(f"{'Genero':<20}{'N1 Occid':<11}" + "".join(f"{c:<13}" for c in ciudades) + "  veredicto")
for g in sorted(N1_ALZHEIMER):
    if g not in T.index: continue
    n1 = N1_ALZHEIMER[g]
    ds = [T.loc[g, f"{c}_N3"] for c in ciudades]
    re = [d for d in ds if d in ("ALTO","BAJO")]
    v = ("sin datos" if not re else
         ">>> INVIERTE vs OCCIDENTE <<<" if all(d != n1 for d in re) else
         "coincide" if all(d == n1 for d in re) else "mixto entre ciudades")
    log(f"{g:<20}{n1:<11}" + "".join(f"{str(d):<13}" for d in ds) + f"  {v}")

log("\n"+"="*92)
log("c) GENEROS CON EFECTO EN >=2 COHORTES, AUSENTES EN LA TABLA")
log("="*92)
todos = set.intersection(*[set(r.columns) for r,_ in datos.values()])
cand = []
for g in sorted(todos):
    ds, ps = [], []
    for nom,(rel,y) in datos.items():
        a = rel.loc[y==1, g].values; c = rel.loc[y==0, g].values
        if a.mean()+c.mean() < 1e-4: continue
        try: _, p = stats.mannwhitneyu(a, c, alternative="two-sided")
        except ValueError: continue
        ds.append(cliff(a,c)); ps.append(p)
    fu = [i for i in range(len(ds)) if abs(ds[i])>=0.3 and ps[i]<0.05]
    if len(fu) >= 2:
        sg = {np.sign(ds[i]) for i in fu}
        cand.append(dict(genero=g, n_cohortes=len(fu),
            direccion="consistente" if len(sg)==1 else "DISCORDANTE",
            delta_medio=np.mean([ds[i] for i in fu]),
            en_tabla_AD=g in N1_ALZHEIMER, en_set_TFM=g in GENEROS_TFM))
C = pd.DataFrame(cand)
if len(C):
    C = C.sort_values("delta_medio", key=abs, ascending=False)
    nv = C[~C.en_tabla_AD & ~C.en_set_TFM]
    log(f"\nCon efecto fuerte en >=2 cohortes: {len(C)}")
    log(f"NO estan en tabla AD ni en set TFM: {len(nv)}")
    log("\nTop 20 candidatos nuevos:")
    log(nv.head(20).to_string(index=False))
    C.to_csv(f"{OUT}/candidatos_nuevos_biomarcadores.csv", index=False)
else:
    log("  ninguno supera el umbral")

T.round(4).to_csv(f"{OUT}/auditoria_direccion_biomarcadores.csv")
log("\n"+"="*92); log("TABLA COMPLETA"); log("="*92)
cols = ["N1_literatura"]+[f"{c}_{s}" for c in ciudades for s in ("N2","N3","p","delta")]
log(T[[c for c in cols if c in T.columns]].round(3).to_string())
log(f"\nSalidas en {OUT}")
