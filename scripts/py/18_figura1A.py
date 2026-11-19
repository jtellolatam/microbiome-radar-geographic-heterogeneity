#!/usr/bin/env python3
"""Microbiome Radar - Figura 1A: diagrama de flujo del cribado"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
BASE = os.path.expanduser("~/Documents/MicrobiomeRadar_Research")
OUT  = os.path.join(BASE, "_figuras")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Helvetica","Arial","DejaVu Sans"], "font.size":9})
AZUL="#2E5A9C"; AZULOSC="#1F3864"; ROJO="#B03030"; FONDO="#F2F4F7"
fig, ax = plt.subplots(figsize=(7.8, 8.6))
ax.set_xlim(0, 10); ax.set_ylim(0, 15.6); ax.axis("off")
def caja(x,y,w,h,texto,sub=None,fill="white",borde=AZULOSC,lw=1.3,fs=9.6,bold=False):
    ax.add_patch(FancyBboxPatch((x,y),w,h,
        boxstyle="round,pad=0.10,rounding_size=0.14",
        facecolor=fill, edgecolor=borde, linewidth=lw, zorder=2))
    if sub:
        ax.text(x+w/2, y+h*0.64, texto, ha="center", va="center", fontsize=fs,
                color="#1A1A1A", fontweight="bold" if bold else "normal", zorder=3)
        ax.text(x+w/2, y+h*0.26, sub, ha="center", va="center",
                fontsize=fs-1.4, color="#555555", zorder=3)
    else:
        ax.text(x+w/2, y+h/2, texto, ha="center", va="center", fontsize=fs,
                color="#1A1A1A", fontweight="bold" if bold else "normal", zorder=3)
def flecha(x,y0,y1):
    ax.add_patch(FancyArrowPatch((x,y0),(x,y1), arrowstyle="-|>",
        mutation_scale=15, linewidth=1.3, color=AZULOSC, zorder=1))
def lateral(y,texto,detalle):
    ax.add_patch(FancyArrowPatch((2.6,y),(4.72,y), arrowstyle="-|>",
        mutation_scale=13, linewidth=1.1, color=ROJO,
        linestyle=(0,(4,2)), zorder=1))
    n=len(detalle); alto=0.62+n*0.36
    ax.add_patch(FancyBboxPatch((4.78, y-alto/2), 5.0, alto,
        boxstyle="round,pad=0.10,rounding_size=0.12",
        facecolor="#FDF3F3", edgecolor=ROJO, linewidth=1.1, zorder=2))
    ax.text(4.98, y+alto/2-0.34, texto, ha="left", va="center",
            fontsize=8.8, color=ROJO, fontweight="bold", zorder=3)
    for i,d in enumerate(detalle):
        ax.text(5.12, y+alto/2-0.72-i*0.35, d, ha="left", va="center",
                fontsize=7.8, color="#555555", zorder=3)
ax.text(0.0, 15.2, "A", fontsize=15, fontweight="bold", color=AZULOSC)
ax.text(0.62, 15.2, "Cribado de cohortes", fontsize=11.5,
        color="#1A1A1A", va="center")
caja(0.3, 13.2, 4.6, 1.30, "Depósitos identificados",
     "ENA · SRA · QIITA · Zenodo        n ≈ 30", fill=FONDO, bold=True)
flecha(2.6, 13.15, 10.8)
lateral(11.98, "Excluidos por incompatibilidad técnica    n = 14",
        ["región 16S distinta a la del brazo   4",
         "shotgun o RNA-Seq   4",
         "plataforma obsoleta o formato propietario   2",
         "sin datos crudos accesibles   2",
         "organismo o gen distinto al declarado   2"])
caja(0.3, 9.5, 4.6, 1.30, "Compatibles técnicamente", "n = 16", bold=True)
flecha(2.6, 9.45, 7.1)
lateral(8.28, "Excluidos por diseño del estudio    n = 6",
        ["intervención o varios puntos temporales",
         "pseudorreplicación",
         "etiquetas de grupo no disponibles"])
caja(0.3, 5.8, 4.6, 1.30, "Cohortes admitidas",
     "n = 10        1.185 muestras", fill="#EAF0F8", borde=AZUL, lw=1.8, bold=True)
brazos = [("Alzheimer","5 cohortes · 6 ciudades\n563 muestras","#2E5A9C"),
          ("TEA","3 cohortes · 3 ciudades\n450 muestras","#2F6B4F"),
          ("Parkinson","2 cohortes · 2 ciudades\n172 muestras","#B5651D")]
ancho=2.95; ytop=4.70; ybox=2.60
xs=[0.35+i*3.22+ancho/2 for i in range(3)]
ax.plot([2.6,2.6],[5.80,ytop], color=AZULOSC, lw=1.3, zorder=1)
ax.plot([min(xs),max(xs)],[ytop,ytop], color=AZULOSC, lw=1.3, zorder=1)
for i,(nom,det,col) in enumerate(brazos):
    x = 0.35 + i*3.22
    ax.add_patch(FancyArrowPatch((xs[i],ytop),(xs[i],ybox+1.72),
        arrowstyle="-|>", mutation_scale=13, linewidth=1.2,
        color=AZULOSC, zorder=1))
    ax.add_patch(FancyBboxPatch((x,ybox), ancho, 1.68,
        boxstyle="round,pad=0.10,rounding_size=0.12",
        facecolor="white", edgecolor=col, linewidth=1.5, zorder=2))
    ax.text(x+ancho/2, ybox+1.24, nom, ha="center", va="center",
            fontsize=9.8, color=col, fontweight="bold", zorder=3)
    ax.text(x+ancho/2, ybox+0.52, det, ha="center", va="center",
            fontsize=8.2, color="#555555", linespacing=1.5, zorder=3)
ax.text(0.3, 1.95,
        "El protocolo de admisión se aplicó sobre los metadatos públicos antes de\n"
        "descargar ficheros de secuencias. Algunas cohortes incumplieron más de un\n"
        "criterio; se indica el motivo principal de exclusión.",
        fontsize=7.8, color="#666666", va="top", linespacing=1.6)
fig.savefig(os.path.join(OUT,"Figura_1A_flujo_cribado.png"),
            dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(OUT,"Figura_1A_flujo_cribado.pdf"),
            bbox_inches="tight", facecolor="white")
plt.close(fig)
print("  Figura_1A_flujo_cribado.png / .pdf")
