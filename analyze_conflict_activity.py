"""
First-differences panel regression of market activity on nearby conflict.

Data: ethiopia_market_activity.json (built by build_ethiopia_data.py), which
holds, per market, a quarterly activity index (2018Q1-2023Q4, market-day
images only) and a quarterly count/fatalities series of ACLED events within
20km of the market.

Model: for market i and quarter t,
    d_act[i,t] = act[i,t] - act[i,t-1]   (only if both quarters observed)
    d_conf[i,t] = conf[i,t] - conf[i,t-1]

  (1) d_act ~ d_conf                      pooled OLS, SEs clustered by market
  (2) d_act ~ d_conf + quarter dummies    absorbs common/seasonal shocks
  (3) same as (2) with d_fatalities instead of event counts
  (4) same as (2) adding lagged d_conf    does conflict predict *next* quarter?

First-differencing removes market fixed effects (baseline size/visibility);
quarter dummies remove nationwide shocks (weather, holidays, sensor mix).
"""
import json

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

with open("ethiopia_market_activity.json", encoding="utf-8") as f:
    data = json.load(f)

quarters = data["quarters"]
qindex = {q: i for i, q in enumerate(quarters)}

rows = []
for mkt_id, m in data["markets"].items():
    conf = m.get("conflict", {})  # missing key -> no nearby events at any point
    for q in quarters:
        act = m["q"].get(q)
        c = conf.get(q, [0, 0])
        rows.append({
            "market": mkt_id,
            "region": m["region"],
            "quarter": q,
            "t": qindex[q],
            "act": act[0] if act else np.nan,
            "n_img": act[1] if act else 0,
            "conf_n": c[0],
            "conf_fat": c[1],
        })

panel = pd.DataFrame(rows).sort_values(["market", "t"]).reset_index(drop=True)
print(f"Panel: {panel['market'].nunique()} markets x {len(quarters)} quarters "
      f"= {len(panel)} cells, {panel['act'].notna().sum()} with activity data")

# first differences, only between consecutive observed quarters
g = panel.groupby("market")
panel["d_act"] = g["act"].diff()
panel["d_conf"] = g["conf_n"].diff()
panel["d_fat"] = g["conf_fat"].diff()
consecutive = g["t"].diff() == 1  # always true given full grid, but explicit
panel.loc[~consecutive, ["d_act", "d_conf", "d_fat"]] = np.nan
panel["d_conf_lag"] = g["d_conf"].shift(1)

est = panel.dropna(subset=["d_act", "d_conf"]).copy()
print(f"Estimation sample: {len(est)} market-quarter differences, "
      f"{est['market'].nunique()} markets\n")

print("Descriptives of the estimation sample:")
print(est[["d_act", "d_conf", "d_fat"]].describe().round(2).to_string(), "\n")

def report(name, res, keys):
    print("=" * 78)
    print(name)
    print("-" * 78)
    for k in keys:
        print(f"  {k:<14} coef = {res.params[k]:>8.4f}   "
              f"se = {res.bse[k]:.4f}   t = {res.tvalues[k]:>6.2f}   "
              f"p = {res.pvalues[k]:.4f}")
    print(f"  N = {int(res.nobs)}   R2 = {res.rsquared:.4f}")
    print()

cluster = {"cov_type": "cluster", "cov_kwds": {"groups": None}}

def fit(formula, df):
    return smf.ols(formula, data=df).fit(
        cov_type="cluster", cov_kwds={"groups": df["market"]})

m1 = fit("d_act ~ d_conf", est)
report("(1) Δactivity ~ Δconflict events (20km), clustered by market", m1, ["d_conf"])

m2 = fit("d_act ~ d_conf + C(quarter)", est)
report("(2) + quarter fixed effects", m2, ["d_conf"])

m3 = fit("d_act ~ d_fat + C(quarter)", est)
report("(3) Δfatalities instead of event counts, + quarter FE", m3, ["d_fat"])

est_lag = est.dropna(subset=["d_conf_lag"])
m4 = fit("d_act ~ d_conf + d_conf_lag + C(quarter)", est_lag)
report("(4) contemporaneous + lagged Δconflict, + quarter FE", m4,
       ["d_conf", "d_conf_lag"])

# effect size framing: what does a large conflict escalation imply?
sd_dconf = est["d_conf"].std()
for name, res in [("(2)", m2)]:
    b = res.params["d_conf"]
    print(f"Effect size: 1 SD increase in Δconflict ({sd_dconf:.1f} events) -> "
          f"{b * sd_dconf:+.2f} points of activity index (spec {name}); "
          f"+10 events -> {b * 10:+.2f} points.")

# save the panel for further work (e.g., Stata, event-study extensions)
est_cols = ["market", "region", "quarter", "act", "conf_n", "conf_fat",
            "d_act", "d_conf", "d_fat", "d_conf_lag", "n_img"]
panel[est_cols].to_csv("conflict_activity_panel.csv", index=False)
print("\nWrote conflict_activity_panel.csv")
