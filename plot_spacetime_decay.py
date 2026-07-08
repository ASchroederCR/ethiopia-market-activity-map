"""
Visualize the space x time decay of the conflict -> market activity effect
(finite distributed-lag ring model from spacetime_decay.py).

Left:  per-event effect by time lag (t, t-1, t-2), one series per distance ring
       -> shows how the effect fades over BOTH distance and time.
Right: contemporaneous vs cumulative (long-run) effect per ring, showing that
       summing the lags gives a larger total effect (conflict persists).

Reads spacetime_decay_rings.csv.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("spacetime_decay_rings.csv")
rings = ["0-10km", "10-20km", "20-50km"]
ring_colors = {"0-10km": "#b2182b", "10-20km": "#ef8a62", "20-50km": "#2166ac"}
lag_labels = {0: "t\n(same qtr)", 1: "t-1\n(+1 qtr)", 2: "t-2\n(+2 qtr)"}

fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.5, 5.2),
                               gridspec_kw={"width_ratios": [1.25, 1]})

# ---------------------------------------------------------------- Panel A
lags = [0, 1, 2]
dodge = {"0-10km": -0.12, "10-20km": 0.0, "20-50km": 0.12}
axL.axhline(0, color="#999", lw=1, zorder=1)
for ring in rings:
    sub = df[(df["ring"] == ring) & (df["lag"].astype(str).isin(["0", "1", "2"]))]
    sub = sub.assign(lag=sub["lag"].astype(int)).sort_values("lag")
    x = [l + dodge[ring] for l in sub["lag"]]
    y = sub["coef"].to_numpy()
    ci = 1.96 * sub["se"].to_numpy()
    axL.errorbar(x, y, yerr=ci, fmt="o-", color=ring_colors[ring],
                 ecolor=ring_colors[ring], elinewidth=1.5, capsize=4,
                 markersize=7, lw=1.4, label=ring, zorder=3)

axL.set_xticks(lags)
axL.set_xticklabels([lag_labels[l] for l in lags])
axL.set_xlim(-0.5, 2.5)
axL.set_xlabel("Time lag of the conflict change", fontsize=10.5)
axL.set_ylabel("Change in activity index per conflict event", fontsize=10.5)
axL.set_title("Effect fades over time as well as distance", fontsize=12, fontweight="bold")
axL.legend(title="Distance ring", fontsize=9, title_fontsize=9, loc="lower right")
axL.grid(axis="y", color="#f2f2f2", zorder=0)
for sp in ("top", "right"):
    axL.spines[sp].set_visible(False)

# ---------------------------------------------------------------- Panel B
# contemporaneous (lag 0) vs cumulative (long-run) per ring
xr = np.arange(len(rings))
bw = 0.36
contemp = [df[(df.ring == r) & (df.lag.astype(str) == "0")]["coef"].iloc[0] for r in rings]
contemp_ci = [1.96 * df[(df.ring == r) & (df.lag.astype(str) == "0")]["se"].iloc[0] for r in rings]
cum = [df[(df.ring == r) & (df.lag.astype(str) == "cumulative")]["coef"].iloc[0] for r in rings]
cum_ci = [1.96 * df[(df.ring == r) & (df.lag.astype(str) == "cumulative")]["se"].iloc[0] for r in rings]

axR.axhline(0, color="#999", lw=1, zorder=1)
axR.bar(xr - bw / 2, contemp, bw, yerr=contemp_ci, capsize=4,
        color="#9ecae1", edgecolor="#3182bd", label="Contemporaneous (t only)", zorder=3)
axR.bar(xr + bw / 2, cum, bw, yerr=cum_ci, capsize=4,
        color="#fc9272", edgecolor="#b2182b", label="Cumulative (t + t-1 + t-2)", zorder=3)
axR.set_xticks(xr)
axR.set_xticklabels(rings)
axR.set_xlabel("Distance ring", fontsize=10.5)
axR.set_ylabel("Effect on activity index per (sustained) event", fontsize=10.5)
axR.set_title("Long-run effect exceeds the same-quarter effect", fontsize=12, fontweight="bold")
axR.legend(fontsize=9, loc="lower right")
axR.grid(axis="y", color="#f2f2f2", zorder=0)
for sp in ("top", "right"):
    axR.spines[sp].set_visible(False)

fig.text(0.5, -0.02,
         "Finite distributed-lag model, first-differenced with quarter FE and market-clustered SEs "
         "(N=28,945 market-quarter changes). Bars/points are per-event marginal effects (±95% CI). "
         "Joint Wald test that all lag terms = 0: F=26.7, p=1.6e-4 (lags jointly significant). "
         "Effect is negative: conflict lowers market activity.",
         ha="center", fontsize=8, color="#666", wrap=True)
fig.suptitle("Conflict's effect on market activity decays over space AND time",
             fontsize=13.5, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig("spacetime_decay.png", dpi=150, bbox_inches="tight", facecolor="white")
print("Wrote spacetime_decay.png")
