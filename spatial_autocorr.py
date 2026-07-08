"""
Are conflict-linked market-activity declines spatially CLUSTERED (nearby markets
decline together = contagion / shared disruption) or do declines coincide with
INCREASES nearby (substitution / trade diversion)?

This tests the spatial autocorrelation of quarterly changes in the activity
index, and whether conflict changes it. Three complementary pieces, all using a
distance-band spatial weight (neighbors within NEIGH_KM, row-standardized):

(1) Global Moran's I of d_act, per quarter (permutation inference).
    I > 0  -> activity changes cluster in space (neighbors move together).
    I < 0  -> checkerboard (a market's change is opposite its neighbors').

(2) Spatial-spillover regressions:
    - own effect:      d_act_i ~ d_conf_i                 (should be negative)
    - SPILLOVER:       Wd_act_i ~ d_conf_i                (KEY: sign answers it)
         Wd_act_i = mean activity change of i's neighbors.
         negative -> i's conflict predicts NEIGHBORS declining  = contagion
         positive -> i's conflict predicts NEIGHBORS rising      = substitution
    - comovement:      d_act_i ~ Wd_act_i * conflict_exposed
         slope on Wd_act = do neighbors' changes move with i's? (contagion if +)

(3) Local Moran (LISA) classification of each market-quarter, to locate hotspots
    and to tabulate, among conflict-linked declines, how often neighbors also
    declined (LL, contagion) vs rose (LH, substitution). The single
    highest-conflict quarter's LISA map is exported for plotting.

Inputs: ethiopia_market_activity.json (activity + per-market nearby-conflict
counts). No raw ACLED needed. Writes several CSVs consumed by
plot_spatial_autocorr.py.
"""
import json

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import sparse

NEIGH_KM = 40          # market-to-market "neighbor" band (trade catchment scale)
N_PERM = 999
RNG = np.random.default_rng(42)

# --- load markets ------------------------------------------------------------
with open("ethiopia_market_activity.json", encoding="utf-8") as f:
    data = json.load(f)

quarters = data["quarters"]
qidx = {q: i for i, q in enumerate(quarters)}
ids = list(data["markets"].keys())
n = len(ids)
lat = np.array([data["markets"][m]["lat"] for m in ids])
lon = np.array([data["markets"][m]["lon"] for m in ids])
region = np.array([data["markets"][m]["region"] for m in ids])

# activity index matrix [markets x quarters], NaN where no data
act = np.full((n, len(quarters)), np.nan)
conf = np.zeros((n, len(quarters)))          # ACLED events within 20km that qtr
for i, m in enumerate(ids):
    md = data["markets"][m]
    for q, rec in md["q"].items():
        act[i, qidx[q]] = rec[0]
    for q, c in (md.get("conflict") or {}).items():
        conf[i, qidx[q]] = c[0]

# quarter-over-quarter change (first difference) of activity and of conflict
d_act = np.full_like(act, np.nan)
d_act[:, 1:] = act[:, 1:] - act[:, :-1]
d_conf = np.zeros_like(conf)
d_conf[:, 1:] = conf[:, 1:] - conf[:, :-1]

# --- spatial weights: distance band, built once -----------------------------
def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0088
    la1, lo1, la2, lo2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat, dlon = la2 - la1, lo2 - lo1
    a = np.sin(dlat/2)**2 + np.cos(la1)*np.cos(la2)*np.sin(dlon/2)**2
    return 2 * R * np.arcsin(np.sqrt(a))

D = haversine_km(lat[:, None], lon[:, None], lat[None, :], lon[None, :])
np.fill_diagonal(D, np.inf)
A = (D <= NEIGH_KM).astype(float)            # binary adjacency
deg = A.sum(1)
print(f"{n} markets; neighbor band <= {NEIGH_KM} km: "
      f"mean {deg.mean():.1f} neighbors, {(deg == 0).sum()} isolates\n")

def row_standardize(Asub):
    d = Asub.sum(1)
    d[d == 0] = 1
    return Asub / d[:, None]

# --- (1) global Moran's I of d_act per quarter -------------------------------
def morans_i(z, W):
    """z: centered values (1D), W: row-standardized dense/sparse weights."""
    lag = W @ z
    num = z @ lag
    den = z @ z
    return num / den if den > 0 else np.nan

moran_rows = []
for q in quarters:
    t = qidx[q]
    avail = np.where(~np.isnan(d_act[:, t]))[0]
    if len(avail) < 30:
        continue
    Asub = A[np.ix_(avail, avail)]
    keep = Asub.sum(1) > 0                    # drop isolates within this quarter
    avail = avail[keep]
    Asub = A[np.ix_(avail, avail)]
    W = row_standardize(Asub)
    z = d_act[avail, t] - d_act[avail, t].mean()
    I = morans_i(z, W)
    # permutation null: shuffle z across locations
    perm = np.empty(N_PERM)
    den = z @ z
    for k in range(N_PERM):
        zp = RNG.permutation(z)
        perm[k] = (zp @ (W @ zp)) / den
    p = (1 + np.sum(perm >= I)) / (N_PERM + 1) if I >= 0 else (1 + np.sum(perm <= I)) / (N_PERM + 1)
    moran_rows.append({"quarter": q, "n": len(avail), "morans_i": I,
                       "p_perm": p, "mean_dconf_abs": np.abs(d_conf[avail, t]).mean(),
                       "total_conf": conf[avail, t].sum()})

moran = pd.DataFrame(moran_rows)
moran.to_csv("spatial_moran_by_quarter.csv", index=False)
sig_pos = ((moran.morans_i > 0) & (moran.p_perm < 0.05)).sum()
sig_neg = ((moran.morans_i < 0) & (moran.p_perm < 0.05)).sum()
print("(1) Global Moran's I of quarterly activity change (d_act)")
print(f"    quarters tested: {len(moran)};  mean I = {moran.morans_i.mean():+.3f}")
print(f"    significant (p<.05): {sig_pos} positive (clustering), {sig_neg} negative (dispersion)")
print(f"    -> {'positive spatial autocorrelation dominates: nearby markets move TOGETHER' if sig_pos > sig_neg else 'negative autocorrelation dominates'}\n")

# --- build long panel with spatial lag Wd_act for regressions ---------------
Wfull = row_standardize(A.copy())
rows = []
for t, q in enumerate(quarters):
    have = ~np.isnan(d_act[:, t])
    # neighbor mean of d_act, over neighbors that have data this quarter
    dz = np.where(have, d_act[:, t], 0.0)
    mask = have.astype(float)
    neigh_sum = A @ dz
    neigh_cnt = A @ mask
    with np.errstate(invalid="ignore", divide="ignore"):
        wd = np.where(neigh_cnt > 0, neigh_sum / neigh_cnt, np.nan)
    for i in range(n):
        if not have[i] or np.isnan(wd[i]):
            continue
        rows.append((ids[i], region[i], q, d_act[i, t], wd[i],
                     conf[i, t], d_conf[i, t]))

panel = pd.DataFrame(rows, columns=["market", "region", "quarter",
                                    "d_act", "Wd_act", "conf", "d_conf"])
panel["conf_exposed"] = (panel["conf"] > 0).astype(int)
print(f"Spillover panel: {len(panel):,} market-quarter obs with a spatial lag\n")

def fit(formula, df):
    return smf.ols(formula, data=df).fit(cov_type="cluster",
                                         cov_kwds={"groups": df["market"]})

reg_rows = []
def report(name, res, keys):
    print(f"  {name}")
    for k in keys:
        b, se, p = res.params[k], res.bse[k], res.pvalues[k]
        print(f"     {k:<26} coef = {b:+.4f}  se = {se:.4f}  p = {p:.2e}")
        reg_rows.append({"model": name, "term": k, "coef": b, "se": se,
                         "t": res.tvalues[k], "p": p, "N": int(res.nobs)})
    print()

print("(2) Spatial-spillover regressions (quarter FE, market-clustered SEs)")
report("own effect: d_act ~ d_conf", fit("d_act ~ d_conf + C(quarter)", panel), ["d_conf"])
report("SPILLOVER: Wd_act ~ d_conf", fit("Wd_act ~ d_conf + C(quarter)", panel), ["d_conf"])
report("SPILLOVER (levels): Wd_act ~ conf", fit("Wd_act ~ conf + C(quarter)", panel), ["conf"])
report("comovement: d_act ~ Wd_act * conf_exposed",
       fit("d_act ~ Wd_act * conf_exposed + C(quarter)", panel),
       ["Wd_act", "Wd_act:conf_exposed", "conf_exposed"])
pd.DataFrame(reg_rows).to_csv("spatial_regressions.csv", index=False)

# --- (3) Local Moran classification + substitution vs contagion tally --------
def lisa_classes(z, W):
    """Return spatial lag and HH/LL/LH/HL label per location (z centered)."""
    lag = W @ z
    lab = np.empty(len(z), dtype=object)
    for i in range(len(z)):
        hi_i, hi_l = z[i] > 0, lag[i] > 0
        lab[i] = ("HH" if hi_i and hi_l else "LL" if not hi_i and not hi_l
                  else "HL" if hi_i and not hi_l else "LH")
    return lag, lab

# pooled: among conflict-exposed own-declines, do neighbors decline or rise?
sub = panel[(panel.conf_exposed == 1) & (panel.d_act < 0)]
non = panel[(panel.conf_exposed == 0) & (panel.d_act < 0)]
def tally(df):
    nb_decline = (df.Wd_act < 0).mean()
    return nb_decline, len(df)
sd_dec, sd_n = tally(sub)
nn_dec, nn_n = tally(non)
print("(3) Among market-quarters where the market itself DECLINED (d_act < 0):")
print(f"    conflict-exposed (event within 20km): neighbors also declined in "
      f"{sd_dec*100:.1f}% of cases  (N={sd_n:,})")
print(f"    NOT conflict-exposed:                 neighbors also declined in "
      f"{nn_dec*100:.1f}% of cases  (N={nn_n:,})")
verdict = ("declines CLUSTER more where there is conflict (contagion / shared "
           "disruption), not substitution" if sd_dec > nn_dec else
           "declines are LESS clustered near conflict (substitution signature)")
print(f"    -> {verdict}\n")

pd.DataFrame([
    {"group": "conflict-exposed declines", "pct_neighbors_also_declined": sd_dec, "N": sd_n},
    {"group": "non-exposed declines", "pct_neighbors_also_declined": nn_dec, "N": nn_n},
]).to_csv("spatial_substitution_summary.csv", index=False)

# LISA map for the single highest-conflict quarter, with per-location
# permutation pseudo-significance so only real local clusters are flagged.
hot_t = int(conf.sum(0).argmax())
hot_q = quarters[hot_t]
avail = np.where(~np.isnan(d_act[:, hot_t]))[0]
Asub = A[np.ix_(avail, avail)]
keep = Asub.sum(1) > 0
avail = avail[keep]
W = row_standardize(A[np.ix_(avail, avail)])
z = d_act[avail, hot_t] - d_act[avail, hot_t].mean()
lag, lab = lisa_classes(z, W)

m2 = (z @ z) / len(z)
local_I = (z / m2) * (W @ z)
ge = np.ones(len(z))                         # permutation counts (|null| >= |obs|)
absobs = np.abs(local_I)
for _ in range(N_PERM):
    zp = RNG.permutation(z)
    li = (zp / m2) * (W @ zp)
    ge += np.abs(li) >= absobs
local_p = ge / (N_PERM + 1)
sig = local_p < 0.05
lab_sig = np.where(sig, lab, "ns")

lisa = pd.DataFrame({
    "market": [ids[i] for i in avail], "lat": lat[avail], "lon": lon[avail],
    "region": region[avail], "d_act": d_act[avail, hot_t],
    "Wd_act_centered": lag, "lisa": lab, "local_p": local_p,
    "lisa_sig": lab_sig, "conf": conf[avail, hot_t],
})
lisa.to_csv("spatial_lisa_hotspot.csv", index=False)
with open("spatial_hotspot_meta.json", "w") as f:
    json.dump({"quarter": hot_q, "neigh_km": NEIGH_KM}, f)
print(f"LISA map exported for highest-conflict quarter {hot_q} "
      f"({len(avail)} markets; significant local clusters at p<.05):")
for cls, desc in [("LL", "decline cluster"), ("HH", "rise cluster"),
                  ("LH", "decline amid rising nbrs"), ("HL", "rise amid declining nbrs")]:
    print(f"    {cls} ({desc}): {int(np.sum(lab_sig == cls))}")
print(f"    ns (not significant): {int(np.sum(lab_sig == 'ns'))}")
print("Wrote spatial_moran_by_quarter.csv, spatial_regressions.csv, "
      "spatial_substitution_summary.csv, spatial_lisa_hotspot.csv")
