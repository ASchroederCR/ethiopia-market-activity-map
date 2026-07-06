"""
Visualize the spatial decay of the conflict -> market activity relationship.

Left panel: a concentric-ring "bullseye" centered on a schematic market, each
ring shaded by the estimated per-event effect on the market activity index
(from the joint ring regression in distance_decay.py). Darker = larger hit.

Right panel: the same ring coefficients plotted as a marginal-effect gradient
vs distance, with 95% confidence intervals, plus the best-fitting exponential
decay kernel (lambda = 20 km) as an illustrative guide to the functional form.

Reads distance_decay_rings.csv and distance_decay_results.csv.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

rings = pd.read_csv("distance_decay_rings.csv")
decay = pd.read_csv("distance_decay_results.csv")

# ring geometry: (inner, outer) km and midpoint
ring_bounds = [(0, 10), (10, 20), (20, 50)]
mids = [(a + b) / 2 for a, b in ring_bounds]
coefs = rings["coef"].to_numpy()          # per-event effect (negative)
ses = rings["se"].to_numpy()
ci = 1.96 * ses

# color scale on the magnitude of the (negative) effect
mag = np.abs(coefs)
norm = Normalize(vmin=0, vmax=mag.max() * 1.15)
cmap = plt.cm.Reds
sm = ScalarMappable(norm=norm, cmap=cmap)

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5.6),
                               gridspec_kw={"width_ratios": [1, 1.15]})

# ------------------------------------------------------------------ Panel A
# bullseye: outer ring first so inner rings draw on top
axL.set_aspect("equal")
lim = 52
for (inner, outer), c in sorted(zip(ring_bounds, coefs), key=lambda z: -z[0][1]):
    color = cmap(norm(abs(c)))
    axL.add_patch(Circle((0, 0), outer, facecolor=color, edgecolor="white", lw=1.5, zorder=1))
# redraw ring boundaries as thin guides + labels
for (inner, outer), mid, c in zip(ring_bounds, mids, coefs):
    axL.add_patch(Circle((0, 0), outer, facecolor="none", edgecolor="#00000022", lw=0.8, zorder=3))
    # label placed in the ring band, upper area
    r_label = (inner + outer) / 2
    axL.text(0, r_label, f"{c:+.2f}", ha="center", va="center",
             fontsize=11, fontweight="bold",
             color="white" if abs(c) > mag.max() * 0.55 else "#7a1010", zorder=4)
    axL.text(0, -r_label, f"{inner}–{outer} km", ha="center", va="center",
             fontsize=8.5, color="white" if abs(c) > mag.max() * 0.55 else "#555", zorder=4)

# central market marker
axL.add_patch(Circle((0, 0), 1.6, facecolor="#1a6ec1", edgecolor="white", lw=1.2, zorder=5))
axL.annotate("market", (0, 0), (0, -lim * 0.92), ha="center", fontsize=9,
             color="#1a6ec1", fontweight="bold",
             arrowprops=dict(arrowstyle="-", color="#1a6ec1", lw=0.8))
axL.set_xlim(-lim, lim)
axL.set_ylim(-lim, lim)
axL.axis("off")
axL.set_title("Per-event effect on market activity, by distance",
              fontsize=12, fontweight="bold", pad=10)

cbar = fig.colorbar(sm, ax=axL, fraction=0.046, pad=0.04)
cbar.set_label("|effect| per conflict event\n(activity-index points)", fontsize=9)

# ------------------------------------------------------------------ Panel B
axR.axhline(0, color="#999", lw=1, ls="-", zorder=1)
axR.errorbar(mids, coefs, yerr=ci, fmt="o", color="#b2182b", ecolor="#b2182b",
             elinewidth=1.6, capsize=5, markersize=8, zorder=4,
             label="Ring estimate (±95% CI)")

# illustrative exponential decay kernel, lambda = 20 km (best AIC), scaled so
# the curve matches the innermost ring's point estimate at its midpoint
lam = 20.0
x = np.linspace(0, 52, 200)
scale = coefs[0] / np.exp(-mids[0] / lam)
axR.plot(x, scale * np.exp(-x / lam), color="#1a6ec1", lw=2, ls="--", zorder=3,
         label=f"Exp. decay kernel (λ={lam:.0f} km, best fit)")

# annotate ring bands
for (inner, outer) in ring_bounds:
    axR.axvspan(inner, outer, color="#00000006", zorder=0)
    axR.axvline(outer, color="#e0e0e0", lw=0.8, zorder=0)

axR.set_xlim(0, 52)
axR.set_xlabel("Distance from market (km)", fontsize=10.5)
axR.set_ylabel("Change in activity index per conflict event", fontsize=10.5)
axR.set_title("Spatial decay of the conflict effect", fontsize=12, fontweight="bold", pad=10)
axR.legend(loc="lower right", fontsize=9, framealpha=0.95)
axR.grid(axis="y", color="#f0f0f0", zorder=0)
for spine in ("top", "right"):
    axR.spines[spine].set_visible(False)

# footnote
fig.text(0.5, -0.02,
         "First-differenced panel (2018Q1–2023Q4, 1,764 markets, N=31,772 market-quarter changes); "
         "quarter fixed effects, SEs clustered by market. "
         "Ring coefficients from a single joint regression. Effect is negative: conflict lowers market activity.",
         ha="center", fontsize=8, color="#666", wrap=True)

fig.suptitle("Conflict events depress nearby rural market activity, with a steep distance gradient",
             fontsize=13.5, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig("distance_decay.png", dpi=150, bbox_inches="tight", facecolor="white")
print("Wrote distance_decay.png")
