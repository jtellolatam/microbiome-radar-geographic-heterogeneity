#!/usr/bin/env python3
"""Microbiome Radar - Figuras 2 a 6, version corregida"""
import os, json, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats as st
warnings.filterwarnings("ignore")
BASE = os.path.expanduser("~/Documents/MicrobiomeRadar_Research")
CC   = os.path.join(BASE, "_crosscohort")
TM   = os.path.join(CC, "transfer_matrix")
OUT  = os.path.join(BASE, "_figuras")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({
    "font.family":"sans-serif",
    "font.sans-serif":["Helvetica","Arial","DejaVu Sans"],
    "font.size":9, "axes.linewidth":0.8, "axes.edgecolor":"#333333",
    "axes.labelcolor":"#1A1A1A", "xtick.color":"#333333", "ytick.color":"#333333",
    "axes.spines.top":False, "axes.spines.right":False, "figure.dpi":110})
AZUL,VERDE,NARANJA = "#2E5A9C","#2F6B4F","#B5651D"
GRIS = "#8A8A8A"
COL = {"Alzheimer":AZUL, "TEA":VERDE, "Parkinson":NARANJA}
REN = {"Suzhou":"Xiamen"}
def fix(s): return REN.get(str(s), str(s))
def guardar(fig,n):
    for e in ("png","pdf"):
        fig.savefig(os.path.join(OUT,f"{n}.{e}"), dpi=300,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig); print(f"  {n}")
def cj(f):
    p=os.path.join(CC,f)
    return json.load(open(p)) if os.path.exists(p) else {}
print("Regenerando figuras...\n")
print("Figura 2")
filas=[]
for brazo,f in [("Alzheimer","auc_5cohortes.json"),("TEA","auc_TEA.json"),
                ("Parkinson","auc_PD.json")]:
    for coh,v in cj(f).items():
        ic=v.get("ic",[np.nan,np.nan])
        filas.append(dict(brazo=brazo,cohorte=fix(coh),auc=v.get("auc",np.nan),
                          lo=ic[0],hi=ic[1],p=v.get("p",np.nan)))
D=pd.DataFrame(filas)
if len(D):
    fig,ax=plt.subplots(figsize=(7.6,5.4))
    y=np.arange(len(D))[::-1]
    for i,(_,r) in enumerate(D.iterrows()):
        yy=y[i]; c=COL[r.brazo]
        sig=(not np.isnan(r.p)) and r.p<0.05
        ax.plot([r.lo,r.hi],[yy,yy],color=c,lw=1.8,
                alpha=0.9 if sig else 0.42, solid_capstyle="round", zorder=2)
        ax.scatter(r.auc,yy,s=52,color=c if sig else "white",
                   edgecolor=c,zorder=4,linewidth=1.5)
        et=f"{r.auc:.3f}" + ("" if sig else "   n.s.")
        ax.text(1.015,yy,et,va="center",fontsize=8,
                color="#1A1A1A" if sig else GRIS, family="monospace")
    ax.axvline(0.5,color=GRIS,ls="--",lw=0.9,zorder=0)
    lim_izq = -0.42
    for brazo,c in COL.items():
        idx=[list(D.index).index(i) for i in D[D.brazo==brazo].index]
        if not idx: continue
        ys=[y[i] for i in idx]
        ax.text(lim_izq, np.mean(ys), brazo, transform=ax.get_yaxis_transform(),
                ha="center", va="center", rotation=90, fontsize=9,
                color=c, fontweight="bold", clip_on=False)
        ax.plot([lim_izq+0.06]*2, [min(ys)-0.34, max(ys)+0.34],
                transform=ax.get_yaxis_transform(), color=c, lw=2.4,
                solid_capstyle="round", clip_on=False)
    ax.set_yticks(y); ax.set_yticklabels(D.cohorte,fontsize=9)
    ax.set_xlim(0.26,1.12); ax.set_xticks([0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
    ax.set_xlabel("AUC — validación cruzada 10×20, IC 95 % bootstrap", fontsize=9)
    ax.set_title("Separación entre casos y controles dentro de cada cohorte",
                 fontsize=11,pad=14,loc="left",color="#1A1A1A")
    fig.subplots_adjust(left=0.30)
    guardar(fig,"Figura_2_rendimiento_interno")
print("Figura 3")
cmap=LinearSegmentedColormap.from_list("t",["#B03030","#F5F5F5","#2E5A9C"])
disp=[(n,f) for n,f in [("Alzheimer","matriz_5cohortes.csv"),
      ("TEA","matriz_TEA.csv"),("Parkinson","matriz_PD.csv")]
      if os.path.exists(os.path.join(TM,f))]
if disp:
    anchos=[pd.read_csv(os.path.join(TM,f),index_col=0).shape[1] for _,f in disp]
    fig,axes=plt.subplots(1,len(disp),
        figsize=(sum(a*0.80 for a in anchos)+3.0, 4.8),
        gridspec_kw={"width_ratios":anchos, "wspace":0.55})
    if len(disp)==1: axes=[axes]
    for ax,(nombre,f) in zip(axes,disp):
        M=pd.read_csv(os.path.join(TM,f),index_col=0)
        M.index=[fix(i) for i in M.index]; M.columns=[fix(c) for c in M.columns]
        V=M.values.astype(float)
        im=ax.imshow(V,cmap=cmap,vmin=0.25,vmax=0.95)
        for i in range(V.shape[0]):
            for j in range(V.shape[1]):
                if np.isnan(V[i,j]): continue
                ax.text(j,i,f"{V[i,j]:.2f}",ha="center",va="center",fontsize=8.2,
                        color="white" if (V[i,j]<0.42 or V[i,j]>0.84) else "#1A1A1A",
                        fontweight="bold" if i==j else "normal")
                if i==j:
                    ax.add_patch(plt.Rectangle((j-.5,i-.5),1,1,fill=False,
                                 edgecolor="#1A1A1A",lw=1.6))
        ax.set_xticks(range(V.shape[1])); ax.set_yticks(range(V.shape[0]))
        ax.set_xticklabels(M.columns,rotation=42,ha="right",fontsize=8)
        ax.set_yticklabels(M.index,fontsize=8)
        ax.set_title(nombre,fontsize=10.5,color=COL[nombre],fontweight="bold",pad=10)
        ax.set_xlabel("modelo aplicado a",fontsize=8.6)
        ax.set_ylabel("modelo derivado en",fontsize=8.6)
        for s in ax.spines.values(): s.set_visible(False)
        ax.tick_params(length=0)
    cb=fig.colorbar(im,ax=axes,fraction=0.016,pad=0.03)
    cb.set_label("AUC",fontsize=8.6); cb.ax.tick_params(labelsize=8)
    fig.suptitle("Transferencia de modelos entre poblaciones",
                 fontsize=11,x=0.015,ha="left",y=1.04,color="#1A1A1A")
    guardar(fig,"Figura_3_transferencia")
print("Figura 4")
paneles=[]
for brazo,f in [("Alzheimer","shap_todos_generos_por_cohorte.csv"),
                ("TEA","shap_TEA_por_cohorte.csv"),
                ("Parkinson","shap_PD_por_cohorte.csv")]:
    p=os.path.join(CC,f)
    if not os.path.exists(p): continue
    S=pd.read_csv(p,index_col=0)
    S.columns=[fix(c) for c in S.columns]
    cols=list(S.columns); pares=[]
    for i,a in enumerate(cols):
        for b in cols[i+1:]:
            rho,_=st.spearmanr(S[a],S[b])
            ta=set(S[a].sort_values(ascending=False).head(20).index)
            tb=set(S[b].sort_values(ascending=False).head(20).index)
            pares.append((f"{a} – {b}",rho,len(ta&tb)))
    if pares: paneles.append((brazo,pares))
if paneles:
    alturas=[len(p) for _,p in paneles]
    fig,axes=plt.subplots(len(paneles),1,
        figsize=(7.4, sum(0.26*a for a in alturas)+2.4),
        gridspec_kw={"height_ratios":alturas, "hspace":0.42})
    if len(paneles)==1: axes=[axes]
    for ax,(brazo,pares) in zip(axes,paneles):
        nom=[x[0] for x in pares]; rhos=[x[1] for x in pares]; ovl=[x[2] for x in pares]
        y=np.arange(len(pares))[::-1]; c=COL[brazo]
        ax.barh(y,rhos,color=[c if r>0 else "#B03030" for r in rhos],
                alpha=0.75,height=0.62)
        for i,o in enumerate(ovl):
            ax.text(0.60,y[i],f"{o}/20",va="center",fontsize=7.6,
                    color="#555555",family="monospace")
        ax.axvline(0,color="#333333",lw=0.9)
        ax.set_yticks(y); ax.set_yticklabels(nom,fontsize=8)
        ax.set_xlim(-0.55,0.78); ax.set_xticks([-0.4,-0.2,0,0.2,0.4])
        ax.set_title(brazo,fontsize=10,color=c,fontweight="bold",pad=6,loc="left")
        if ax is axes[-1]:
            ax.set_xlabel("ρ de Spearman entre rankings de importancia",fontsize=8.6)
    fig.suptitle("Concordancia de descriptores entre cohortes  ·  solapamiento del top-20 a la derecha",
                 fontsize=11,x=0.015,ha="left",y=1.0,color="#1A1A1A")
    guardar(fig,"Figura_4_concordancia")
print("Figura 5")
disp=[]
for brazo,f,ca,cb,la,lb in [
    ("Alzheimer","contraste_oriente_occidente.csv","d_oriente","d_occidente","Oriente","Occidente"),
    ("TEA","contraste_TEA_china_eeuu.csv","d_china","d_eeuu","China","EE. UU."),
    ("Parkinson","direcciones_PD.csv","d_Finlandia","d_Rumania","Finlandia","Rumanía")]:
    p=os.path.join(CC,f)
    if not os.path.exists(p): continue
    D2=pd.read_csv(p)
    if ca in D2.columns and cb in D2.columns: disp.append((brazo,D2,ca,cb,la,lb))
if disp:
    fig,axes=plt.subplots(1,len(disp),figsize=(3.9*len(disp),4.3))
    if len(disp)==1: axes=[axes]
    for k,(ax,(brazo,D2,ca,cb,la,lb)) in enumerate(zip(axes,disp)):
        x=D2[ca].astype(float); y=D2[cb].astype(float); c=COL[brazo]
        tr=D2["transversal"] if "transversal" in D2.columns else pd.Series(False,index=D2.index)
        inv=D2["inversion"] if "inversion" in D2.columns else (
            D2["patron"]=="INVIERTE" if "patron" in D2.columns else pd.Series(False,index=D2.index))
        base=~(tr|inv)
        ax.scatter(x[base],y[base],s=13,color=GRIS,alpha=0.30,linewidth=0)
        ax.scatter(x[tr],y[tr],s=46,color=c,alpha=0.9,linewidth=0,
                   label=f"concordante ({int(tr.sum())})")
        ax.scatter(x[inv],y[inv],s=46,facecolor="white",edgecolor="#B03030",
                   linewidth=1.4,label=f"inversión ({int(inv.sum())})")
        lim=max(abs(np.nanmin([x.min(),y.min()])),abs(np.nanmax([x.max(),y.max()])))*1.14
        ax.plot([-lim,lim],[-lim,lim],color=GRIS,ls="--",lw=0.8,zorder=0)
        ax.axhline(0,color="#333333",lw=0.7); ax.axvline(0,color="#333333",lw=0.7)
        ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
        ax.set_xlabel(f"delta de Cliff · {la}",fontsize=8.6)
        ax.set_ylabel(f"delta de Cliff · {lb}",fontsize=8.6)
        ax.set_title(f"{brazo}   n = {len(D2)} géneros",fontsize=10,
                     color=c,fontweight="bold",pad=10,loc="left")
        ax.legend(fontsize=7.4,frameon=False,loc="lower right",
                  handletextpad=0.5, borderpad=0.2)
        ax.set_aspect("equal")
    fig.suptitle("Tamaño de efecto por región  ·  la diagonal marca concordancia perfecta",
                 fontsize=11,x=0.015,ha="left",y=1.03,color="#1A1A1A")
    fig.subplots_adjust(wspace=0.38)
    guardar(fig,"Figura_5_contraste_regional")
print("Figura 6")
p=os.path.join(CC,"geografia_vs_condicion_resumen.csv")
if os.path.exists(p):
    G=pd.read_csv(p)
    fig,ax=plt.subplots(figsize=(7.4,4.8)); x=np.arange(len(G)); w=0.26
    b1=ax.bar(x-w,G.A_caso_vs_control,w,color="#C9CDD4",edgecolor="#333333",lw=0.7)
    b2=ax.bar(x,G.B_origen_controles,w,color=AZUL,alpha=0.85,edgecolor="#333333",lw=0.7)
    b3=ax.bar(x+w,G.C_origen_casos,w,color=NARANJA,alpha=0.85,edgecolor="#333333",lw=0.7)
    for i,r in G.iterrows():
        for off,v in [(-w,r.A_caso_vs_control),(0,r.B_origen_controles),
                      (w,r.C_origen_casos)]:
            ax.text(i+off,v+0.014,f"{v:.3f}",ha="center",fontsize=7.6,family="monospace")
    ax.axhline(0.5,color=GRIS,ls="--",lw=0.9,zorder=0)
    ax.text(len(G)-0.42, 0.515, "azar", fontsize=7.4, color=GRIS, va="bottom")
    ax.set_xticks(x); ax.set_xticklabels(G.brazo,fontsize=9.6)
    ax.set_ylim(0,1.04); ax.set_ylabel("AUC",fontsize=9)
    ax.legend([b1,b2,b3],
              ["caso frente a control, dentro de una cohorte",
               "cohorte de origen · solo controles",
               "cohorte de origen · solo casos"],
              fontsize=8.2, frameon=False, ncol=3,
              loc="upper center", bbox_to_anchor=(0.5,-0.10),
              handlelength=1.4, columnspacing=1.4)
    ax.set_title("El origen de la muestra predomina sobre la condición",
                 fontsize=11,pad=14,loc="left",color="#1A1A1A")
    guardar(fig,"Figura_6_origen_vs_condicion")
print(f"\nFiguras en: {OUT}")
