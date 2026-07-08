"""
Visualize the bivariate spatial association (own conflict x neighbors' activity
change) from spatial_bivariate.py.

Left:  global bivariate Moran's I per quarter (negative => conflict amid
       declining neighbors), with nearby-conflict intensity overlaid.
Right: local bivariate LISA map for the peak-conflict quarter, highlighting
       significant "conflict hotspot, neighbors declining" cells (the spillover)
       against the other quadrants.
"""
import json

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

biv = pd.read_csv("spatial_bivariate_by_quarter.csv")
hot = pd.read_csv("spatial_bivariate_hotspot.csv")
meta = json.load(open("spatial_bivariate_meta.json"))
hot_q = meta["quarter"]
adm1 = gpd.read_file("ethiopia_adm1.geojson").to_crs(4326)
ev = pd.DataFrame(json.load(open("ethiopia_conflict_events.json", encoding="utf-8"))["events"])
ev_hot = ev[ev["q"] == hot_q]

fig = plt.figure(figsize=(14, 6.2))
gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1])

# ----------------------------------------------------------------- Panel A
axA = fig.add_subplot(gs[0, 0])
x = np.arange(len(biv))
sig = biv["p_perm"] < 0.05
colors = np.where(sig & (biv.biv_moran_i < 0), "#6a1b9a", "#cbb8dd")
axA.bar(x, biv["biv_moran_i"], color=colors, zorder=3, width=0.7)
axA.axhline(0, color="#666", lw=1)
axA.set_ylabel("Bivariate Moran's I\n(own conflict × neighbors' Δactivity)", fontsize=10)
axA.set_title("Conflict sits amid declining neighbors (negative = spillover)",
              fontsize=11.5, fontweight="bold")
axA.set_xticks([i for i, q in enumerate(biv["quarter"]) if q.endswith("Q1")])
axA.set_xticklabels([q[:4] for q in biv["quarter"] if q.endswith("Q1")])
axA.spines["top"].set_visible(False)

axA2 = axA.twinx()
axA2.plot(x, biv["total_conf"], color="#333", lw=1.3, marker="o", markersize=3, zorder=4)
axA2.set_ylabel("Total nearby ACLED events (line)", fontsize=9.5, color="#333")
axA2.spines["top"].set_visible(False)
axA.legend(handles=[
    Line2D([0], [0], marker="s", color="none", markerfacecolor="#6a1b9a", markersize=9,
           label="I_B < 0, p<.05 (spillover)"),
    Line2D([0], [0], marker="s", color="none", markerfacecolor="#cbb8dd", markersize=9,
           label="not significant"),
    Line2D([0], [0], color="#333", lw=1.3, marker="o", markersize=3,
           label="total nearby ACLED events"),
], loc="lower left", fontsize=8.3, framealpha=0.95)

# ----------------------------------------------------------------- Panel B
axB = fig.add_subplot(gs[0, 1])
adm1.boundary.plot(ax=axB, color="#cccccc", linewidth=0.6, zorder=1)

ns = hot[hot["quad_sig"] == "ns"]
axB.scatter(ns["lon"], ns["lat"], s=4, c="#dddddd", zorder=2, linewidths=0)

style = {
    "conf_hi_nbr_down": ("#b2182b", "Conflict hotspot, neighbors declining (spillover)"),
    "conf_hi_nbr_up":   ("#f4a582", "Conflict hotspot, neighbors rising"),
    "conf_lo_nbr_down": ("#92c5de", "Low conflict, neighbors declining"),
    "conf_lo_nbr_up":   ("#2166ac", "Low conflict, neighbors rising"),
}
for k, (c, _) in style.items():
    sub = hot[hot["quad_sig"] == k]
    z = 5 if k == "conf_hi_nbr_down" else 3
    axB.scatter(sub["lon"], sub["lat"], s=28 if k == "conf_hi_nbr_down" else 22,
                c=c, edgecolors="white", linewidths=0.4, zorder=z)

axB.scatter(ev_hot["lon"], ev_hot["lat"], s=9, marker="s", c="#000000",
            alpha=0.22, zorder=4, linewidths=0)

axB.set_title(f"Where conflict drags down neighbors, {hot_q}", fontsize=11.5, fontweight="bold")
axB.set_xlabel("Longitude", fontsize=9.5)
axB.set_ylabel("Latitude", fontsize=9.5)
axB.set_aspect(1 / np.cos(np.radians(9)))
axB.set_xlim(32.5, 48)
axB.set_ylim(3, 15.2)

items = [Line2D([0], [0], marker="o", color="none", markerfacecolor=c,
                markeredgecolor="white", markersize=9, label=desc)
         for k, (c, desc) in style.items()]
items += [Line2D([0], [0], marker="s", color="none", markerfacecolor="#000000",
                 alpha=0.3, markersize=8, label="ACLED conflict event"),
          Line2D([0], [0], marker="o", color="none", markerfacecolor="#dddddd",
                 markersize=7, label="no significant association")]
axB.legend(handles=items, loc="lower left", fontsize=7.4, framealpha=0.95)

fig.text(0.5, -0.02,
         "Bivariate local Moran's I: each market's own nearby-conflict count vs the mean activity change of "
         "its neighbors (within 40 km, row-standardized), both z-scored per quarter; significance from 999 "
         "permutations of the neighbor variable. In the peak-conflict quarter, significant conflict hotspots "
         "with declining neighbors (113) outnumber those with rising neighbors (8) ~14 to 1.",
         ha="center", fontsize=8, color="#666", wrap=True)
fig.suptitle("Bivariate LISA: a market's conflict predicts its neighbors' decline",
             fontsize=13, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig("spatial_bivariate.png", dpi=150, bbox_inches="tight", facecolor="white")
print("Wrote spatial_bivariate.png")
