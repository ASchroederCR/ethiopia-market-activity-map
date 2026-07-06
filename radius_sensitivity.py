"""
Radius sensitivity for the conflict -> market activity first-difference model.

Recomputes the market<->ACLED linkage at 10, 20, and 50 km catchments and
re-runs the main specification for each:

    d_act[i,t] ~ d_conf[i,t] + quarter FE,  SEs clustered by market

Activity data comes from ethiopia_market_activity.json (unchanged across
radii); conflict counts are rebuilt from the raw ACLED CSV per radius.
"""
import json

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

REPO = "../repo"
RADII_KM = [10, 20, 50]

# --- activity panel (radius-independent) -----------------------------------
with open("ethiopia_market_activity.json", encoding="utf-8") as f:
    data = json.load(f)

quarters = data["quarters"]
qindex = {q: i for i, q in enumerate(quarters)}
market_ids = list(data["markets"].keys())
mkt_lat = np.array([data["markets"][m]["lat"] for m in market_ids])
mkt_lon = np.array([data["markets"][m]["lon"] for m in market_ids])

rows = []
for mkt_id in market_ids:
    m = data["markets"][mkt_id]
    for q in quarters:
        act = m["q"].get(q)
        rows.append({
            "market": mkt_id,
            "quarter": q,
            "t": qindex[q],
            "act": act[0] if act else np.nan,
        })
panel = pd.DataFrame(rows)

# --- ACLED events -----------------------------------------------------------
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

# one distance matrix, thresholded per radius
dist_km = haversine_km(mkt_lat[:, None], mkt_lon[:, None],
                       ev_lat[None, :], ev_lon[None, :])

def fit(formula, df):
    return smf.ols(formula, data=df).fit(
        cov_type="cluster", cov_kwds={"groups": df["market"]})

print(f"Panel: {len(market_ids)} markets x {len(quarters)} quarters; "
      f"{len(conflict)} ACLED events in range\n")

results = []
for radius in RADII_KM:
    within = dist_km <= radius

    # conflict count per market-quarter at this radius
    conf_rows = []
    for i, mkt_id in enumerate(market_ids):
        idxs = np.nonzero(within[i])[0]
        counts = pd.Series(ev_q[idxs]).value_counts() if len(idxs) else pd.Series(dtype=int)
        for q in quarters:
            conf_rows.append({"market": mkt_id, "quarter": q,
                              "conf_n": int(counts.get(q, 0))})
    conf_df = pd.DataFrame(conf_rows)

    df = panel.merge(conf_df, on=["market", "quarter"]).sort_values(["market", "t"])
    g = df.groupby("market")
    df["d_act"] = g["act"].diff()
    df["d_conf"] = g["conf_n"].diff()
    est = df.dropna(subset=["d_act", "d_conf"])

    res = fit("d_act ~ d_conf + C(quarter)", est)
    b, se, p = res.params["d_conf"], res.bse["d_conf"], res.pvalues["d_conf"]
    sd = est["d_conf"].std()
    share_exposed = (df.groupby("market")["conf_n"].sum() > 0).mean()
    results.append({
        "radius_km": radius,
        "coef_per_event": b,
        "se": se,
        "p": p,
        "sd_dconf": sd,
        "effect_1sd": b * sd,
        "pct_markets_exposed": 100 * share_exposed,
        "N": int(res.nobs),
    })
    print(f"radius {radius:>2}km: coef = {b:+.4f} (se {se:.4f}, p {p:.2e})   "
          f"SD(d_conf) = {sd:5.2f} -> 1SD effect = {b*sd:+.2f} pts   "
          f"{100*share_exposed:.0f}% of markets exposed   N = {int(res.nobs)}")

out = pd.DataFrame(results)
out.to_csv("radius_sensitivity_results.csv", index=False)
print("\nWrote radius_sensitivity_results.csv")
