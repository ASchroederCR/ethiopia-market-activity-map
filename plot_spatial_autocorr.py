"""
Visualize the spatial-autocorrelation analysis (spatial_autocorr.py):

Left:  quarterly global Moran's I of the activity-index change -- consistently
       positive => nearby markets' changes move together (spatial clustering),
       with total nearby-conflict per quarter overlaid.
Right: Local Moran (LISA) hotspot map for the highest-conflict quarter, showing
       significant decline-clusters (LL) vs rise-clusters (HH) and the rarer
       substitution cells (LH/HL), with ACLED conflict events overlaid.
"""
import json

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

moran = pd.read_csv("spatial_moran_by_quarter.csv")
lisa = pd.read_csv("spatial_lisa_hotspot.csv")
meta = json.load(open("spatial_hotspot_meta.json"))
hot_q = meta["quarter"]
adm1 = gpd.read_file("ethiopia_adm1.geojson").to_crs(4326)
conf_events = json.load(open("ethiopia_conflict_events.json", encoding="utf-8"))
ev = pd.DataFrame(conf_events["events"])
ev_hot = ev[ev["q"] == hot_q]

fig = plt.figure(figsize=(14, 6.2))
gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1])

# ----------------------------------------------------------------- Panel A
axA = fig.add_subplot(gs[0, 0])
x = np.arange(len(moran))
sig = moran["p_perm"] < 0.05
colors = np.where(sig & (moran.morans_i > 0), "#b2182b",
                  np.where(sig & (moran.morans_i < 0), "#2166ac", "#bbbbbb"))
axA.bar(x, moran["morans_i"], color=colors, zorder=3, width=0.7)
axA.axhline(0, color="#666", lw=1)
axA.set_ylabel("Global Moran's I of Δactivity", fontsize=10.5)
axA.set_title("Market-activity changes cluster in space, every quarter",
              fontsize=12, fontweight="bold")
axA.set_xticks([i for i, q in enumerate(moran["quarter"]) if q.endswith("Q1")])
axA.set_xticklabels([q[:4] for q in moran["quarter"] if q.endswith("Q1")])
for sp in ("top",):
    axA.spines[sp].set_visible(False)

# conflict intensity overlay
axA2 = axA.twinx()
axA2.plot(x, moran["total_conf"], color="#333", lw=1.4, marker="o",
          markersize=3, zorder=4, label="ACLED events near markets")
axA2.set_ylabel("Total nearby ACLED events (line)", fontsize=9.5, color="#333")
axA2.spines["top"].set_visible(False)
axA.legend(handles=[
    Line2D([0], [0], marker="s", color="none", markerfacecolor="#b2182b", markersize=9,
           label="Moran's I > 0, p<.05 (clustering)"),
    Line2D([0], [0], color="#333", lw=1.4, marker="o", markersize=3,
           label="Total nearby ACLED events"),
], loc="upper left", fontsize=8.5, framealpha=0.95)

# ----------------------------------------------------------------- Panel B
axB = fig.add_subplot(gs[0, 1])
adm1.boundary.plot(ax=axB, color="#cccccc", linewidth=0.6, zorder=1)

# non-significant markets: faint grey
ns = lisa[lisa["lisa_sig"] == "ns"]
axB.scatter(ns["lon"], ns["lat"], s=4, c="#dddddd", zorder=2, linewidths=0)

cls_style = {
    "LL": ("#b2182b", "Decline cluster (market ↓, neighbors ↓) — contagion"),
    "HH": ("#2166ac", "Rise cluster (market ↑, neighbors ↑)"),
    "LH": ("#f4a582", "Market ↓ amid rising neighbors — substitution"),
    "HL": ("#92c5de", "Market ↑ amid declining neighbors"),
}
for cls, (color, _) in cls_style.items():
    sub = lisa[lisa["lisa_sig"] == cls]
    axB.scatter(sub["lon"], sub["lat"], s=26, c=color, edgecolors="white",
                linewidths=0.4, zorder=4 if cls in ("LL", "HH") else 3)

# conflict events overlaid as small dark squares
axB.scatter(ev_hot["lon"], ev_hot["lat"], s=13, marker="s", c="#000000",
            alpha=0.35, zorder=5, linewidths=0)

axB.set_title(f"Local hotspots of correlated change, {hot_q} (peak conflict)",
              fontsize=12, fontweight="bold")
axB.set_xlabel("Longitude", fontsize=9.5)
axB.set_ylabel("Latitude", fontsize=9.5)
axB.set_aspect(1 / np.cos(np.radians(9)))
axB.set_xlim(32.5, 48)
axB.set_ylim(3, 15.2)

legend_items = [Line2D([0], [0], marker="o", color="none", markerfacecolor=c,
                       markeredgecolor="white", markersize=9, label=desc)
                for cls, (c, desc) in cls_style.items()]
legend_items.append(Line2D([0], [0], marker="s", color="none", markerfacecolor="#000000",
                           alpha=0.4, markersize=8, label="ACLED conflict event"))
legend_items.append(Line2D([0], [0], marker="o", color="none", markerfacecolor="#dddddd",
                           markersize=7, label="market, no significant cluster"))
axB.legend(handles=legend_items, loc="lower left", fontsize=7.6, framealpha=0.95)

fig.text(0.5, -0.02,
         "Spatial weights: markets within 40 km, row-standardized. Moran's I via 999 permutations. "
         "LISA quadrants shown only where locally significant (p<.05, 999 permutations). Δactivity is the "
         "quarter-over-quarter change in the activity index; declines cluster together near conflict rather "
         "than diverting trade to neighbors (which would show as LH substitution cells).",
         ha="center", fontsize=8, color="#666", wrap=True)
fig.suptitle("Conflict-linked market declines are spatially contagious, not substitutive",
             fontsize=13.5, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig("spatial_autocorr.png", dpi=150, bbox_inches="tight", facecolor="white")
print("Wrote spatial_autocorr.png")
