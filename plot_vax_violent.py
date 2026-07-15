"""
One figure to summarise the conflict-definition robustness: the MCV1-conflict
correlation under every conflict definition and adjustment we tried. The naive
hypothesis predicts a NEGATIVE correlation (conflict -> low vaccination); across
all definitions the correlation is POSITIVE, weakening (never reversing) as we
strip out urban protests and adjust for population density.
"""
import numpy as np
import matplotlib.pyplot as plt

# cluster-level Pearson r, from analyze_vax_conflict / _confounder / _fatalities / _violent
rows = [
    ("Events, all types — raw",              0.57, "#b2182b"),
    ("Events, all types — adj. pop density", 0.37, "#b2182b"),
    ("Events, armed violence only — raw",    0.36, "#ef8a62"),
    ("Events, armed violence only — adj.",   0.32, "#ef8a62"),
    ("Fatalities, all types — raw",          0.21, "#6a51a3"),
    ("Fatalities, armed violence — raw",     0.16, "#6a51a3"),
]
labels = [r[0] for r in rows][::-1]
vals   = [r[1] for r in rows][::-1]
colors = [r[2] for r in rows][::-1]
y = np.arange(len(labels))

fig, ax = plt.subplots(figsize=(11, 5))
# shade the region the hypothesis predicts
ax.axvspan(-0.7, 0, color="#1a6ec1", alpha=0.06, zorder=0)
ax.text(-0.35, len(labels) - 0.4, 'hypothesis predicts\nnegative (conflict → low MCV1)',
        ha="center", va="top", fontsize=9, color="#1a6ec1", style="italic")
ax.axvline(0, color="#333", lw=1.2, zorder=2)
ax.barh(y, vals, color=colors, zorder=3, height=0.66)
for yi, v in zip(y, vals):
    ax.text(v + 0.012, yi, f"+{v:.2f}", va="center", fontsize=9.5, fontweight="bold")
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9.5)
ax.set_xlim(-0.7, 0.7)
ax.set_xlabel("Correlation of MCV1 coverage with nearby conflict (cluster level, n=722)", fontsize=10.5)
ax.set_title("Every conflict definition gives a POSITIVE MCV1–conflict correlation",
             fontsize=12.5, fontweight="bold", pad=12)
for sp in ("top", "right", "left"):
    ax.spines[sp].set_visible(False)
ax.tick_params(left=False)

fig.text(0.5, -0.04,
         "Restricting ACLED to armed violence (Battles, Violence against civilians, Explosions/Remote "
         "violence — 73% of events; dropping Protests, Riots, Strategic developments) weakens the confounded "
         "positive correlation but does not reverse it. The vaccination deficit remains in the low-violence "
         "pastoralist lowlands. Associations are ecological and temporally offset (see document).",
         ha="center", fontsize=8, color="#666", wrap=True)
fig.tight_layout()
fig.savefig("vax_conflict_definitions.png", dpi=150, bbox_inches="tight", facecolor="white")
print("Wrote vax_conflict_definitions.png")
