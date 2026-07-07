"""
Per-region distance-decay of the conflict -> market activity relationship.

Re-estimates the concentric-ring model (0-10 / 10-20 / 20-50 km, first
differenced, quarter FE, SEs clustered by market) separately for each admin-1
region with enough markets for reliable clustered inference, and draws a
small-multiples grid of the per-event marginal-effect gradient so the spatial
decay can be compared across regions. A national panel is included as baseline.

Regions with too few markets (clusters) for market-clustered SEs are skipped
and listed in the figure footnote.

Reads ethiopia_market_activity.json + the raw ACLED CSV; writes
distance_decay_by_region.png and distance_decay_by_region.csv.
"""
import json

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt

REPO = "../repo"
RING_BOUNDS = [(0, 10), (10, 20), (20, 50)]
MIDS = [(a + b) / 2 for a, b in RING_BOUNDS]
MIN_MARKETS = 150  # need enough clusters for market-clustered SEs to be credible;
                   # below this the ring CIs are too wide to interpret (e.g.
                   # Beneshangul Gumu at 39 markets has CIs of +/-4 to +/-8)

# --- activity panel ----------------------------------------------------------
with open("ethiopia_market_activity.json", encoding="utf-8") as f:
    data = json.load(f)

quarters = data["quarters"]
qindex = {q: i for i, q in enumerate(quarters)}
market_ids = list(data["markets"].keys())
mkt_lat = np.array([data["markets"][m]["lat"] for m in market_ids])
mkt_lon = np.array([data["markets"][m]["lon"] for m in market_ids])
mkt_region = np.array([data["markets"][m]["region"] for m in market_ids])

panel = pd.DataFrame([
    {"market": mkt_id, "region": data["markets"][mkt_id]["region"],
     "quarter": q, "t": qindex[q],
     "act": (data["markets"][mkt_id]["q"].get(q) or [np.nan])[0]}
    for mkt_id in market_ids for q in quarters
])

# --- ACLED events ------------------------------------------------------------
conflict = pd.read_csv(f"{REPO}/datasets/conflict/2012-07-01-2025-07-01-Ethiopia.csv")
conflict["event_date"] = pd.to_datetime(conflict["event_date"], format="mixed")
conflict["quarter"] = conflict["event_date"].dt.to_period("Q").astype(str)
conflict = conflict[conflict["quarter"].isin(quarters)].reset_index(drop=True)
ev_lat = conflict["latitude"].to_numpy()
ev_lon = conflict["longitude"].to_numpy()
ev_q = conflict["quarter"].to_numpy()
ev_idx_by_q = {q: np.nonzero(ev_q == q)[0] for q in quarters}

def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0088
    lat1r, lon1r, lat2r, lon2r = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat, dlon = lat2r - lat1r, lon2r - lon1r
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1r) * np.cos(lat2r) * np.sin(dlon / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))

# full market x event distance matrix (events can be in any region; we count
# whatever falls within the ring distance of each market)
dist_km = haversine_km(mkt_lat[:, None], mkt_lon[:, None],
                       ev_lat[None, :], ev_lon[None, :])

# precompute ring membership masks once
ring_masks = [((dist_km > lo) & (dist_km <= hi)) for lo, hi in RING_BOUNDS]

def ring_counts_for(mask_rows):
    """Return long DataFrame market x quarter with a count column per ring,
    for the market rows selected by boolean mask_rows."""
    sel = np.nonzero(mask_rows)[0]
    sel_ids = [market_ids[i] for i in sel]
    frames = []
    for q in quarters:
        eidx = ev_idx_by_q[q]
        row = {"market": sel_ids, "quarter": q}
        for r, m in enumerate(ring_masks):
            if len(eidx):
                row[f"r{r}"] = m[np.ix_(sel, eidx)].sum(axis=1)
            else:
                row[f"r{r}"] = np.zeros(len(sel))
        frames.append(pd.DataFrame(row))
    return pd.concat(frames, ignore_index=True)

def estimate(mask_rows):
    """Fit the joint ring model on the selected markets; return coef/CI per ring."""
    counts = ring_counts_for(mask_rows)
    df = panel.merge(counts, on=["market", "quarter"]).sort_values(["market", "t"])
    g = df.groupby("market")
    df["d_act"] = g["act"].diff()
    for r in range(len(RING_BOUNDS)):
        df[f"d_r{r}"] = g[f"r{r}"].diff()
    est = df.dropna(subset=["d_act"] + [f"d_r{r}" for r in range(len(RING_BOUNDS))])
    # drop rings with no within-region variation to avoid singular design
    terms = [f"d_r{r}" for r in range(len(RING_BOUNDS)) if est[f"d_r{r}"].std() > 1e-9]
    if not terms:
        return None
    res = smf.ols("d_act ~ " + " + ".join(terms) + " + C(quarter)",
                  data=est).fit(cov_type="cluster",
                                cov_kwds={"groups": est["market"]})
    out = {"n_markets": int(mask_rows.sum()), "n_obs": int(res.nobs)}
    for r in range(len(RING_BOUNDS)):
        k = f"d_r{r}"
        if k in terms:
            out[f"coef{r}"] = res.params[k]
            out[f"ci{r}"] = 1.96 * res.bse[k]
        else:
            out[f"coef{r}"] = np.nan
            out[f"ci{r}"] = np.nan
    return out

# --- estimate national + qualifying regions ---------------------------------
region_counts = pd.Series(mkt_region).value_counts()
regions = [r for r in region_counts.index if region_counts[r] >= MIN_MARKETS]
excluded = [f"{r} ({region_counts[r]})" for r in region_counts.index
            if region_counts[r] < MIN_MARKETS]

results = {}
results["All Ethiopia"] = estimate(np.ones(len(market_ids), dtype=bool))
for r in regions:
    results[r] = estimate(mkt_region == r)

# order panels: national first, then regions by market count (desc)
panel_order = ["All Ethiopia"] + sorted(regions, key=lambda r: -region_counts[r])

# save table
rows = []
for name in panel_order:
    res = results[name]
    rows.append({"region": name, "n_markets": res["n_markets"], "n_obs": res["n_obs"],
                 **{f"coef_{lo}_{hi}km": res[f"coef{i}"]
                    for i, (lo, hi) in enumerate(RING_BOUNDS)},
                 **{f"ci95_{lo}_{hi}km": res[f"ci{i}"]
                    for i, (lo, hi) in enumerate(RING_BOUNDS)}})
pd.DataFrame(rows).to_csv("distance_decay_by_region.csv", index=False)

# --- figure: small multiples -------------------------------------------------
# fixed, readable y-window; a few CIs for the smaller/noisier regions extend
# beyond it and are clipped at the panel edge (noted in the footnote).
ylo, yhi = -3.5, 1.5
nat = results["All Ethiopia"]
nat_coefs = [nat[f"coef{r}"] for r in range(len(RING_BOUNDS))]

ncol = 3
nrow = int(np.ceil(len(panel_order) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(4.3 * ncol, 3.5 * nrow),
                         sharey=True, sharex=True)
axes = np.atleast_1d(axes).ravel()

for ax, name in zip(axes, panel_order):
    res = results[name]
    coefs = [res[f"coef{r}"] for r in range(len(RING_BOUNDS))]
    cis = [res[f"ci{r}"] for r in range(len(RING_BOUNDS))]
    is_national = name == "All Ethiopia"
    color = "#333" if is_national else "#b2182b"

    for lo, hi in RING_BOUNDS:
        ax.axvspan(lo, hi, color="#00000005", zorder=0)
        ax.axvline(hi, color="#eee", lw=0.8, zorder=0)
    ax.axhline(0, color="#999", lw=1, zorder=1)
    # faint national gradient as a common reference in every regional panel
    if not is_national:
        ax.plot(MIDS, nat_coefs, "o-", color="#bbb", lw=1, markersize=3,
                zorder=2, label="National")
    ax.errorbar(MIDS, coefs, yerr=cis, fmt="o-", color=color, ecolor=color,
                elinewidth=1.4, capsize=4, markersize=6, lw=1.3, zorder=3)

    ax.set_title(f"{name}\n{res['n_markets']} markets, N={res['n_obs']:,}",
                 fontsize=10.5, fontweight="bold" if is_national else "normal")
    ax.set_xlim(0, 52)
    ax.set_ylim(ylo, yhi)
    if not is_national:
        ax.legend(loc="upper left", fontsize=7.5, framealpha=0.9,
                  handlelength=1.2, borderpad=0.3)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(axis="y", color="#f4f4f4", zorder=0)

# hide any unused axes
for ax in axes[len(panel_order):]:
    ax.axis("off")

# shared axis labels
fig.supxlabel("Distance from market (km)", fontsize=11, y=0.04)
fig.supylabel("Change in activity index per conflict event", fontsize=11, x=0.005)

foot = ("Concentric-ring model per admin-1 region (0-10 / 10-20 / 20-50 km), first-differenced, "
        "quarter FE, SEs clustered by market. Points are per-event marginal effects (±95% CI); grey line "
        "repeats the national gradient for reference. Some CIs for the smaller regions extend beyond the "
        f"axis. Regions with <{MIN_MARKETS} markets omitted for unreliable clustered inference: "
        + ", ".join(excluded) + ".")
fig.text(0.5, -0.01, foot, ha="center", fontsize=8, color="#666", wrap=True)

fig.suptitle("Conflict–market-activity distance decay varies by region",
             fontsize=14, fontweight="bold", y=1.02)
fig.tight_layout(rect=[0.02, 0.05, 1, 1])
fig.savefig("distance_decay_by_region.png", dpi=150, bbox_inches="tight", facecolor="white")
print("Wrote distance_decay_by_region.png and distance_decay_by_region.csv")
print("Regions estimated:", panel_order)
print("Excluded (too few markets):", excluded)
