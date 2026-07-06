"""
Distance-decay estimation of the conflict -> market activity relationship.

Two specifications, both first-differenced with quarter FE and market-clustered
SEs (same design as analyze_conflict_activity.py / radius_sensitivity.py):

(A) Concentric rings, one model:
      d_act ~ d_conf[0-10km] + d_conf[10-20km] + d_conf[20-50km] + quarter FE
    The ring counts are mutually exclusive, so the coefficients trace the
    spatial gradient directly and can be compared within a single regression.

(B) Continuous exponential decay:
      exposure[i,t] = sum over events j in quarter t of exp(-dist_ij / lambda)
    estimated separately for several decay lengths lambda; model fit (R2/AIC)
    indicates which spatial scale best describes the data. An event at
    distance lambda counts ~0.37 of an adjacent event.
"""
import json

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

REPO = "../repo"
RINGS = [(0, 10), (10, 20), (20, 50)]
LAMBDAS_KM = [5, 10, 20, 40]

# --- activity panel ----------------------------------------------------------
with open("ethiopia_market_activity.json", encoding="utf-8") as f:
    data = json.load(f)

quarters = data["quarters"]
qindex = {q: i for i, q in enumerate(quarters)}
market_ids = list(data["markets"].keys())
mkt_lat = np.array([data["markets"][m]["lat"] for m in market_ids])
mkt_lon = np.array([data["markets"][m]["lon"] for m in market_ids])

panel = pd.DataFrame([
    {"market": mkt_id, "quarter": q, "t": qindex[q],
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

def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0088
    lat1r, lon1r, lat2r, lon2r = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat, dlon = lat2r - lat1r, lon2r - lon1r
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1r) * np.cos(lat2r) * np.sin(dlon / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))

dist_km = haversine_km(mkt_lat[:, None], mkt_lon[:, None],
                       ev_lat[None, :], ev_lon[None, :])
ev_idx_by_q = {q: np.nonzero(ev_q == q)[0] for q in quarters}

def quarterly_series(weight_matrix):
    """Sum a market-x-event weight matrix over events, per quarter.
    Returns a DataFrame market x quarter -> value (long format)."""
    out = []
    for q in quarters:
        idx = ev_idx_by_q[q]
        vals = weight_matrix[:, idx].sum(axis=1) if len(idx) else np.zeros(len(market_ids))
        out.append(pd.DataFrame({"market": market_ids, "quarter": q, "v": vals}))
    return pd.concat(out, ignore_index=True)

def fit(formula, df):
    return smf.ols(formula, data=df).fit(
        cov_type="cluster", cov_kwds={"groups": df["market"]})

def first_diff(df, cols):
    df = df.sort_values(["market", "t"]).copy()
    g = df.groupby("market")
    for c in cols:
        df["d_" + c] = g[c].diff()
    return df

print(f"Panel: {len(market_ids)} markets x {len(quarters)} quarters; "
      f"{len(conflict)} ACLED events\n")

# =============================================================================
# (A) Ring specification
# =============================================================================
df = panel.copy()
ring_cols = []
for lo, hi in RINGS:
    col = f"conf_{lo}_{hi}"
    ring_cols.append(col)
    w = ((dist_km > lo) & (dist_km <= hi)).astype(float)
    s = quarterly_series(w).rename(columns={"v": col})
    df = df.merge(s, on=["market", "quarter"])

df = first_diff(df, ["act"] + ring_cols)
est = df.dropna(subset=["d_act"] + ["d_" + c for c in ring_cols])

formula = "d_act ~ " + " + ".join("d_" + c for c in ring_cols) + " + C(quarter)"
resA = fit(formula, est)

print("=" * 78)
print("(A) Concentric rings, single model: d_act ~ ring counts + quarter FE")
print("-" * 78)
for c in ring_cols:
    k = "d_" + c
    lo, hi = c.split("_")[1:]
    print(f"  {lo:>3}-{hi:<3} km    coef = {resA.params[k]:+.4f}   "
          f"se = {resA.bse[k]:.4f}   t = {resA.tvalues[k]:>6.2f}   "
          f"p = {resA.pvalues[k]:.2e}")
print(f"  N = {int(resA.nobs)}   R2 = {resA.rsquared:.4f}\n")

# =============================================================================
# (B) Exponential decay exposure at several decay lengths
# =============================================================================
print("=" * 78)
print("(B) Exponential-decay exposure: d_act ~ d_exposure(lambda) + quarter FE")
print("-" * 78)
rowsB = []
for lam in LAMBDAS_KM:
    w = np.exp(-dist_km / lam)
    s = quarterly_series(w).rename(columns={"v": "expo"})
    dfl = first_diff(panel.merge(s, on=["market", "quarter"]), ["act", "expo"])
    estl = dfl.dropna(subset=["d_act", "d_expo"])
    res = fit("d_act ~ d_expo + C(quarter)", estl)
    b, se, p = res.params["d_expo"], res.bse["d_expo"], res.pvalues["d_expo"]
    sd = estl["d_expo"].std()
    rowsB.append({"lambda_km": lam, "coef": b, "se": se, "p": p,
                   "effect_1sd": b * sd, "r2": res.rsquared, "aic": res.aic,
                   "N": int(res.nobs)})
    print(f"  lambda = {lam:>2} km   coef = {b:+.4f} (se {se:.4f}, p {p:.2e})   "
          f"1SD effect = {b*sd:+.2f} pts   R2 = {res.rsquared:.4f}   "
          f"AIC = {res.aic:,.0f}")

best = min(rowsB, key=lambda r: r["aic"])
print(f"\nBest fit by AIC: lambda = {best['lambda_km']} km")

pd.DataFrame(rowsB).to_csv("distance_decay_results.csv", index=False)

ringA = pd.DataFrame({
    "ring": [f"{lo}-{hi}km" for lo, hi in RINGS],
    "coef": [resA.params["d_" + c] for c in ring_cols],
    "se": [resA.bse["d_" + c] for c in ring_cols],
    "p": [resA.pvalues["d_" + c] for c in ring_cols],
})
ringA.to_csv("distance_decay_rings.csv", index=False)
print("Wrote distance_decay_results.csv and distance_decay_rings.csv")
