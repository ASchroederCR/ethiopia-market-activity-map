"""
Space x time decay of the conflict -> market activity relationship.

Extends distance_decay.py by adding 1- and 2-quarter time lags to the spatial
model, i.e. a finite distributed-lag (FDL) model in first differences:

    d_act[i,t] = sum_ring sum_{L=0,1,2} beta[ring,L] * d_conf_ring[i,t-L]
                 + quarter FE + e[i,t]

(and the analogous FDL for the continuous exponential-decay exposure). All
specs keep the first-difference + quarter-FE design with market-clustered SEs.

Interpretation of the FDL-in-differences: a *permanent* one-unit step in a
ring's conflict count raises/lowers the change in activity by beta[ring,0] in
the same quarter, beta[ring,1] the next, beta[ring,2] two quarters out; the
cumulative sum beta[ring,0]+[,1]+[,2] is the long-run effect on the *level* of
the activity index of that sustained conflict increase.

Validity checks reported:
  - each ring x lag coefficient (is a given lag significant?);
  - per-ring cumulative effect with its own SE (linear combination);
  - a joint Wald test that ALL lag terms (L=1,2) are zero -> "are lags needed?";
  - AIC of the lagged vs contemporaneous-only model on the SAME sample.

Writes spacetime_decay_rings.csv, spacetime_decay_exp.csv.
"""
import json

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

REPO = "../repo"
RINGS = [(0, 10), (10, 20), (20, 50)]
LAMBDAS_KM = [5, 10, 20, 40]
MAXLAG = 2

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
    out = []
    for q in quarters:
        idx = ev_idx_by_q[q]
        vals = weight_matrix[:, idx].sum(axis=1) if len(idx) else np.zeros(len(market_ids))
        out.append(pd.DataFrame({"market": market_ids, "quarter": q, "v": vals}))
    return pd.concat(out, ignore_index=True)

def add_diff_and_lags(df, base_cols):
    """First-difference activity and each base col, then add L1..MAXLAG lags of
    the conflict differences, respecting market grouping and quarter order."""
    df = df.sort_values(["market", "t"]).copy()
    g = df.groupby("market")
    df["d_act"] = g["act"].diff()
    for c in base_cols:
        df["d_" + c] = g[c].diff()
    for c in base_cols:
        for L in range(1, MAXLAG + 1):
            df[f"d_{c}_L{L}"] = df.groupby("market")["d_" + c].shift(L)
    return df

def fit(formula, df):
    return smf.ols(formula, data=df).fit(
        cov_type="cluster", cov_kwds={"groups": df["market"]})

print(f"Panel: {len(market_ids)} markets x {len(quarters)} quarters; "
      f"{len(conflict)} ACLED events; max lag = {MAXLAG} quarters\n")

# =============================================================================
# (A) Ring FDL:  d_act ~ (ring x {L0,L1,L2}) + quarter FE
# =============================================================================
df = panel.copy()
ring_cols = []
for lo, hi in RINGS:
    col = f"conf_{lo}_{hi}"
    ring_cols.append(col)
    w = ((dist_km > lo) & (dist_km <= hi)).astype(float)
    df = df.merge(quarterly_series(w).rename(columns={"v": col}),
                  on=["market", "quarter"])

df = add_diff_and_lags(df, ring_cols)

# terms
contemp_terms = [f"d_{c}" for c in ring_cols]
lag_terms = [f"d_{c}_L{L}" for c in ring_cols for L in range(1, MAXLAG + 1)]
all_conf_terms = [f"d_{c}_L{L}" if L else f"d_{c}"
                  for c in ring_cols for L in range(0, MAXLAG + 1)]

est = df.dropna(subset=["d_act"] + contemp_terms + lag_terms).copy()

full_formula = "d_act ~ " + " + ".join(all_conf_terms) + " + C(quarter)"
resA = fit(full_formula, est)
# contemporaneous-only baseline on the SAME sample (nested, comparable AIC)
resA0 = fit("d_act ~ " + " + ".join(contemp_terms) + " + C(quarter)", est)

print("=" * 82)
print("(A) Ring finite-distributed-lag model:  d_act ~ ring x {L0,L1,L2} + quarter FE")
print(f"    Estimation sample: N = {int(resA.nobs):,} market-quarter obs "
      f"(t>={MAXLAG+1} within market), {est['market'].nunique()} markets")
print("-" * 82)
print(f"  {'ring':<10}{'lag':<6}{'coef':>10}{'se':>9}{'t':>8}{'p':>10}")
rowsA = []
for c in ring_cols:
    lo, hi = c.split("_")[1:]
    for L in range(0, MAXLAG + 1):
        k = f"d_{c}" if L == 0 else f"d_{c}_L{L}"
        b, se, t, p = resA.params[k], resA.bse[k], resA.tvalues[k], resA.pvalues[k]
        star = "***" if p < .01 else "**" if p < .05 else "*" if p < .1 else ""
        print(f"  {lo+'-'+hi+'km':<10}{('t-'+str(L)) if L else 't':<6}"
              f"{b:>10.4f}{se:>9.4f}{t:>8.2f}{p:>10.2e} {star}")
        rowsA.append({"ring": f"{lo}-{hi}km", "lag": L, "coef": b, "se": se,
                      "t": t, "p": p})
    # cumulative (long-run) effect for this ring: sum of its 3 coefficients
    terms = [f"d_{c}"] + [f"d_{c}_L{L}" for L in range(1, MAXLAG + 1)]
    tt = resA.t_test(" + ".join(terms) + " = 0")
    cb = float(np.ravel(tt.effect)[0])
    cse = float(np.ravel(tt.sd)[0])
    cp = float(np.ravel(tt.pvalue)[0])
    print(f"  {lo+'-'+hi+'km':<10}{'cum':<6}{cb:>10.4f}{cse:>9.4f}"
          f"{cb/cse:>8.2f}{cp:>10.2e}  <- long-run effect of a sustained change")
    rowsA.append({"ring": f"{lo}-{hi}km", "lag": "cumulative", "coef": cb,
                  "se": cse, "t": cb / cse, "p": cp})

# joint validity test: are the lag terms (L1,L2 for all rings) jointly zero?
wald = resA.wald_test(lag_terms, scalar=True)
print("-" * 82)
print(f"  Joint Wald test that ALL {len(lag_terms)} lag terms = 0 "
      f"(i.e. no time lags needed):")
print(f"     F = {float(wald.statistic):.2f},  p = {float(wald.pvalue):.3e}  "
      f"-> {'lags ARE jointly significant' if float(wald.pvalue) < .05 else 'lags NOT significant'}")
print(f"  AIC: contemporaneous-only = {resA0.aic:,.0f}   with lags = {resA.aic:,.0f}   "
      f"(lower is better; delta = {resA.aic - resA0.aic:+,.0f})")
print(f"  R2:  contemporaneous-only = {resA0.rsquared:.4f}   with lags = {resA.rsquared:.4f}\n")

pd.DataFrame(rowsA).to_csv("spacetime_decay_rings.csv", index=False)

# =============================================================================
# (B) Exponential-decay FDL at several decay lengths
# =============================================================================
print("=" * 82)
print("(B) Exponential-decay exposure FDL:  d_act ~ expo(lambda) x {L0,L1,L2} + quarter FE")
print("-" * 82)
print(f"  {'lambda':<8}{'L0':>9}{'L1':>9}{'L2':>9}{'cumul':>10}{'lag p(joint)':>14}{'AIC':>12}")
rowsB = []
for lam in LAMBDAS_KM:
    w = np.exp(-dist_km / lam)
    dfl = panel.merge(quarterly_series(w).rename(columns={"v": "expo"}),
                      on=["market", "quarter"])
    dfl = add_diff_and_lags(dfl, ["expo"])
    lags = [f"d_expo_L{L}" for L in range(1, MAXLAG + 1)]
    terms = ["d_expo"] + lags
    estl = dfl.dropna(subset=["d_act"] + terms)
    res = fit("d_act ~ " + " + ".join(terms) + " + C(quarter)", estl)
    b = [res.params[t] for t in terms]
    tt = res.t_test(" + ".join(terms) + " = 0")
    cum, cump = float(np.ravel(tt.effect)[0]), float(np.ravel(tt.pvalue)[0])
    wj = res.wald_test(lags, scalar=True)
    star = "***" if cump < .01 else "**" if cump < .05 else "*" if cump < .1 else ""
    print(f"  {str(lam)+'km':<8}{b[0]:>9.3f}{b[1]:>9.3f}{b[2]:>9.3f}"
          f"{cum:>10.3f}{float(wj.pvalue):>14.2e}{res.aic:>12,.0f}  {star}")
    rowsB.append({"lambda_km": lam, "coef_L0": b[0], "coef_L1": b[1],
                  "coef_L2": b[2], "cumulative": cum, "cumulative_p": cump,
                  "lag_joint_p": float(wj.pvalue), "aic": res.aic,
                  "N": int(res.nobs)})

pd.DataFrame(rowsB).to_csv("spacetime_decay_exp.csv", index=False)
print("\nWrote spacetime_decay_rings.csv and spacetime_decay_exp.csv")
