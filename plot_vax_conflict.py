"""
The counter-intuitive headline: low vaccination and 2018-23 conflict/markets sit
in largely DIFFERENT places. Two panels make the point and the confound explicit.
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

reg = pd.read_csv("vax_region_summary.csv")
zt = pd.read_csv("zone_integration_summary.csv")
z = zt[zt["mcv1_smooth"].notna()].copy()

fig, (axA, axB) = plt.subplots(1, 2, figsize=(14, 6))

# ---- Panel A: regional scatter MCV1 vs nearby conflict --------------------
axA.scatter(reg["conflict_near"], reg["mcv1"] * 100, s=reg["n_clusters"] * 4,
            c="#6a51a3", alpha=0.75, edgecolors="white", zorder=3)
for _, r in reg.iterrows():
    axA.annotate(r["region"], (r["conflict_near"], r["mcv1"] * 100),
                 fontsize=8, xytext=(4, 4), textcoords="offset points")
# trend line
b, a = np.polyfit(reg["conflict_near"], reg["mcv1"] * 100, 1)
xs = np.linspace(reg["conflict_near"].min(), reg["conflict_near"].max(), 50)
axA.plot(xs, a + b * xs, color="#b2182b", lw=1.5, ls="--", zorder=2)
axA.set_xlabel("Mean ACLED conflict events near clusters (2018-2023)", fontsize=10.5)
axA.set_ylabel("MCV1 coverage (%)", fontsize=10.5)
axA.set_title("More conflict → HIGHER vaccination (a confound, not a cause)",
              fontsize=11.5, fontweight="bold")
axA.text(0.03, 0.06,
         "Raw corr = +0.57;  partial = +0.37\nafter adjusting for population density.\n"
         "Deficit is in low-conflict pastoralist\nlowlands (Somali, Afar), not conflict zones.",
         transform=axA.transAxes, fontsize=8.5, color="#444",
         bbox=dict(boxstyle="round", fc="#f7f7f7", ec="#ddd"))
for sp in ("top", "right"):
    axA.spines[sp].set_visible(False)

# ---- Panel B: zones sorted by MCV1, conflict/market overlaid --------------
zz = z.sort_values("mcv1_smooth").reset_index(drop=True)
y = np.arange(len(zz))
axB.barh(y, zz["mcv1_smooth"] * 100, color="#9e9ac8", zorder=3, height=0.8,
         label="MCV1 coverage (%)")
axB.set_yticks([])
axB.set_xlabel("MCV1 coverage (%)  —  74 admin-2 zones, low to high", fontsize=10.5)
axB.set_ylabel("zones (each bar = one zone)", fontsize=9.5)
axB.set_title("Where vaccination is lowest, markets & conflict are scarce",
              fontsize=11.5, fontweight="bold")
# overlay: mark zones that have any market and any conflict
for i, r in zz.iterrows():
    if r["n_markets"] > 0:
        axB.plot(101, i, marker="o", ms=3, color="#1a6ec1", zorder=4)
    if r["conflict_events"] > zz["conflict_events"].median():
        axB.plot(105, i, marker="s", ms=3.5, color="#b2182b", alpha=0.6, zorder=4)
axB.set_xlim(0, 112)
axB.axhspan(-0.5, 9.5, color="#3f007d", alpha=0.06, zorder=0)
axB.text(40, 4.5, "lowest-coverage zones\n(Somali lowlands):\nfew/no markets, little conflict",
         fontsize=8.2, color="#3f007d", va="center")
axB.text(101, len(zz) + 0.5, "market", fontsize=7.5, color="#1a6ec1", rotation=90, va="bottom", ha="center")
axB.text(105, len(zz) + 0.5, "conflict", fontsize=7.5, color="#b2182b", rotation=90, va="bottom", ha="center")
from matplotlib.lines import Line2D
axB.legend(handles=[
    Line2D([0], [0], marker="o", color="none", markerfacecolor="#1a6ec1", markersize=6, label="zone has ≥1 market"),
    Line2D([0], [0], marker="s", color="none", markerfacecolor="#b2182b", markersize=6, label="zone above-median conflict"),
], loc="center right", fontsize=8, framealpha=0.95)
for sp in ("top", "right"):
    axB.spines[sp].set_visible(False)

fig.text(0.5, -0.02,
         "Vaccination: EDHS 2024-25 MCV1 (children 12-23mo), child-weighted to admin-2 zones (>=3 DHS clusters). "
         "Market activity & ACLED conflict: 2018-2023. Cross-sectional and temporally offset — associations are "
         "ecological, not causal. See accompanying document for full caveats.",
         ha="center", fontsize=8, color="#666", wrap=True)
fig.suptitle("Low measles vaccination does NOT coincide with conflict or market decline in Ethiopia",
             fontsize=13.5, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig("vax_conflict_market.png", dpi=150, bbox_inches="tight", facecolor="white")
print("Wrote vax_conflict_market.png")
