#!/usr/bin/env python3
"""Microbiome Radar - Figura 1B: mapa en proyeccion Equal Earth"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
BASE = os.path.expanduser("~/Documents/MicrobiomeRadar_Research")
OUT  = os.path.join(BASE, "_figuras")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Helvetica","Arial","DejaVu Sans"], "font.size":9})
AZUL,VERDE,NARANJA = "#2E5A9C","#2F6B4F","#B5651D"
COL = {"Alzheimer":AZUL,"TEA":VERDE,"Parkinson":NARANJA}
AZULOSC = "#1F3864"
SITIOS = [
    ("Shanghái",   31.23, 121.47, "Alzheimer",  90, False, ( 20,   2)),
    ("Lishui",     28.47, 119.92, "Alzheimer", 171, False, ( 20, -14)),
    ("Chongqing",  29.56, 106.55, "Alzheimer",  86, False, (-18,   6)),
    ("Estambul",   41.01,  28.98, "Alzheimer",  62, False, ( 12, -16)),
    ("Kayseri",    38.73,  35.49, "Alzheimer",  63, False, ( 20,   4)),
    ("Brescia",    45.54,  10.22, "Alzheimer",  91, False, (-18,  -8)),
    ("Nanjing",    32.06, 118.80, "TEA",       285, False, ( 18,  24)),
    ("Xiamen",     24.48, 118.09, "TEA",        83, False, (  4, -24)),
    ("Pasadena",   34.15,-118.14, "TEA",        82, True,  (-18,   8)),
    ("Helsinki",   60.17,  24.94, "Parkinson", 133, True,  (-18,  10)),
    ("Iași",       47.16,  27.59, "Parkinson",  39, True,  ( 18,  14)),
]
proj = ccrs.EqualEarth(central_longitude=40)
fig = plt.figure(figsize=(10.6, 5.6))
ax = plt.axes(projection=proj)
ax.set_global()
ax.add_feature(cfeature.LAND.with_scale("110m"),
               facecolor="#E7E9EC", edgecolor="#C4C9D0", linewidth=0.45, zorder=1)
ax.add_feature(cfeature.OCEAN.with_scale("110m"), facecolor="#FAFBFC", zorder=0)
ax.add_feature(cfeature.BORDERS.with_scale("110m"),
               edgecolor="#D3D7DD", linewidth=0.32, zorder=2)
ax.gridlines(draw_labels=False, linewidth=0.35, color="#DFE3E8", alpha=0.8, zorder=1)
ax.spines["geo"].set_edgecolor("#C4C9D0"); ax.spines["geo"].set_linewidth(0.7)
for nombre, lat, lon, brazo, n, inf, (dx, dy) in SITIOS:
    c = COL[brazo]
    ax.scatter(lon, lat, s=28+n*0.44, color=c, alpha=0.80,
               edgecolor="white", linewidth=1.0, zorder=6,
               transform=ccrs.PlateCarree())
    ax.annotate(nombre + ("*" if inf else ""),
                xy=(lon, lat), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax),
                xytext=(dx, dy), textcoords="offset points",
                fontsize=8.4, color="#1A1A1A", zorder=7,
                ha="left" if dx > 0 else "right", va="center",
                arrowprops=dict(arrowstyle="-", color="#9AA0A8",
                                lw=0.65, shrinkA=0, shrinkB=4))
fig.text(0.055, 0.955, "B", fontsize=15, fontweight="bold", color=AZULOSC)
fig.text(0.088, 0.955, "Localizaciones de reclutamiento", fontsize=11.5,
         color="#1A1A1A", va="center")
lg = fig.add_axes([0.055, 0.02, 0.90, 0.10])
lg.set_xlim(0, 100); lg.set_ylim(0, 10); lg.axis("off")
for i, (brazo, c) in enumerate(COL.items()):
    x = 2 + i*14
    lg.scatter(x, 6.6, s=64, color=c, alpha=0.80,
               edgecolor="white", linewidth=1.0, zorder=5)
    lg.text(x+2.2, 6.6, brazo, fontsize=8.8, va="center", color="#1A1A1A")
lg.text(58, 9.2, "el tamaño del punto es proporcional al número de muestras",
        fontsize=8, color="#666666", ha="left", va="center")
for i, (n, lab) in enumerate([(39,"39"), (133,"133"), (285,"285")]):
    x = 62 + i*11
    lg.scatter(x, 5.4, s=28+n*0.44, color="#B8BEC6", alpha=0.75,
               edgecolor="white", linewidth=0.9, zorder=5)
    lg.text(x, 1.4, lab, fontsize=7.6, ha="center", color="#666666")
fig.text(0.055, -0.035,
         "* Ciudad no declarada en los metadatos; se indica la sede de la institución responsable.",
         fontsize=7.6, color="#666666", va="center")
fig.text(0.055, -0.065,
         "Proyección Equal Earth (Šavrič, Patterson y Jenny, 2018). Cartografía: Natural Earth, dominio público.",
         fontsize=7.2, color="#8A8A8A", va="center")
fig.savefig(os.path.join(OUT,"Figura_1B_mapa.png"),
            dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(OUT,"Figura_1B_mapa.pdf"),
            bbox_inches="tight", facecolor="white")
plt.close(fig)
print("  Figura_1B_mapa.png / .pdf")
print(f"\nGuardada en: {OUT}")
