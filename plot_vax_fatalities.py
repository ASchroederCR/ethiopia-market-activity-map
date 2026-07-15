"""
Fatalities robustness figure. Left: MCV1 vs conflict fatalities by region — the
positive (confounded) link is weaker than with event counts. Right: MCV1 vs
fatalities-PER-EVENT (lethality) by zone — the one measure that hints at the
"expected" negative direction, though marginal and still confounded.
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

reg = pd.read_csv("vax_region_fatalities.csv")
zt = pd.read_csv("zone_integration_summary.csv")
z = zt[zt["mcv1_smooth"].notna()].copy()
z["fatal_per_event"] = np.where(z["conflict_events"] > 0, z["fatalities"] / z["conflict_events"], np.nan)
z = z.dropna(subset=["fatal_per_event"])

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5.6))

# Panel A: region MCV1 vs fatalities
axA.scatter(reg["fatalities"], reg["mcv1"] * 100, s=90, c="#8c2d04",
            alpha=0.8, edgecolors="white", zorder=3)
for _, r in reg.iterrows():
    axA.annotate(r["region"], (r["fatalities"], r["mcv1"] * 100), fontsize=8,
                 xytext=(4, 4), textcoords="offset points")
b, a = np.polyfit(reg["fatalities"], reg["mcv1"] * 100, 1)
xs = np.linspace(reg["fatalities"].min(), reg["fatalities"].max(), 50)
axA.plot(xs, a + b * xs, color="#b2182b", lw=1.4, ls="--", zorder=2)
axA.set_xlabel("Mean conflict FATALITIES near clusters (2018-2023)", fontsize=10.5)
axA.set_ylabel("MCV1 coverage (%)", fontsize=10.5)
axA.set_title("MCV1 vs fatalities: positive link is weaker than with event counts",
              fontsize=11, fontweight="bold")
axA.text(0.97, 0.06,
         "cluster raw r: events +0.57 → fatalities +0.21\n"
         "after density + region FE: fatalities β≈0 (p=0.60)",
         transform=axA.transAxes, ha="right", fontsize=8.3, color="#444",
         bbox=dict(boxstyle="round", fc="#f7f7f7", ec="#ddd"))
for sp in ("top", "right"):
    axA.spines[sp].set_visible(False)

# Panel B: zone MCV1 vs fatalities per event (lethality)
axB.scatter(z["fatal_per_event"], z["mcv1_smooth"] * 100, s=45,
            c="#6a51a3", alpha=0.75, edgecolors="white", zorder=3)
b2, a2 = np.polyfit(z["fatal_per_event"], z["mcv1_smooth"] * 100, 1)
xs2 = np.linspace(z["fatal_per_event"].min(), z["fatal_per_event"].max(), 50)
axB.plot(xs2, a2 + b2 * xs2, color="#b2182b", lw=1.4, ls="--", zorder=2)
axB.set_xlabel("Fatalities per conflict event (lethality)", fontsize=10.5)
axB.set_ylabel("MCV1 coverage (%)", fontsize=10.5)
axB.set_title("Higher lethality-per-event → slightly LOWER vaccination",
              fontsize=11, fontweight="bold")
axB.text(0.97, 0.93, "zone r = −0.24 (p=0.06)\nthe one hint toward the\n\"expected\" direction — but\nmarginal & still confounded",
         transform=axB.transAxes, ha="right", va="top", fontsize=8.3, color="#444",
         bbox=dict(boxstyle="round", fc="#f7f7f7", ec="#ddd"))
for sp in ("top", "right"):
    axB.spines[sp].set_visible(False)

fig.text(0.5, -0.02,
         "Fatalities and event counts are only moderately correlated across regions (r=0.42): urban areas "
         "(e.g. Addis) log many low-lethality events, while the northern war produced high fatalities. Neither "
         "measure supports 'conflict → low vaccination' once population density and region are accounted for.",
         ha="center", fontsize=8, color="#666", wrap=True)
fig.suptitle("Robustness: conflict FATALITIES instead of event counts",
             fontsize=13, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig("vax_fatalities.png", dpi=150, bbox_inches="tight", facecolor="white")
print("Wrote vax_fatalities.png")
