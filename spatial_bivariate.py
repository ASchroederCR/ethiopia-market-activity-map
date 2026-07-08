"""
Bivariate spatial association between a market's OWN conflict and its
NEIGHBORS' activity change -- i.e. does conflict at location i drag down the
markets around it? This maps the spillover directly, complementing the
univariate LISA in spatial_autocorr.py.

Bivariate Moran / LISA with:
    x_i = nearby ACLED conflict (events within 20km) at market i, quarter t
    y_i = quarter-over-quarter activity change (d_act) at market i
    (both z-standardized within the quarter; W = markets within 40km, row-std)

  Global:  I_B = z_x' W z_y / n      (negative => conflict areas sit amid
                                       neighbors whose activity is FALLING)
  Local:   I_B,i = z_x,i * (W z_y)_i
     quadrant by (sign z_x,i, sign of neighbor mean (W z_y)_i):
        (+, -)  conflict HIGH, neighbors DECLINING  <- the spillover hotspot
        (+, +)  conflict HIGH, neighbors RISING
        (-, -)  conflict LOW,  neighbors DECLINING
        (-, +)  conflict LOW,  neighbors RISING

Writes spatial_bivariate_by_quarter.csv and spatial_bivariate_hotspot.csv
(the peak-conflict quarter, matching the univariate map).
"""
import json

import numpy as np
import pandas as pd
from scipy import sparse

NEIGH_KM = 40
N_PERM = 999
RNG = np.random.default_rng(42)

with open("ethiopia_market_activity.json", encoding="utf-8") as f:
    data = json.load(f)
quarters = data["quarters"]
qidx = {q: i for i, q in enumerate(quarters)}
ids = list(data["markets"].keys())
n = len(ids)
lat = np.array([data["markets"][m]["lat"] for m in ids])
lon = np.array([data["markets"][m]["lon"] for m in ids])
region = np.array([data["markets"][m]["region"] for m in ids])

act = np.full((n, len(quarters)), np.nan)
conf = np.zeros((n, len(quarters)))
for i, m in enumerate(ids):
    md = data["markets"][m]
    for q, rec in md["q"].items():
        act[i, qidx[q]] = rec[0]
    for q, c in (md.get("conflict") or {}).items():
        conf[i, qidx[q]] = c[0]
d_act = np.full_like(act, np.nan); d_act[:, 1:] = act[:, 1:] - act[:, :-1]

def haversine_km(a1, o1, a2, o2):
    R = 6371.0088
    a1, o1, a2, o2 = map(np.radians, [a1, o1, a2, o2])
    dlat, dlon = a2 - a1, o2 - o1
    h = np.sin(dlat/2)**2 + np.cos(a1)*np.cos(a2)*np.sin(dlon/2)**2
    return 2 * R * np.arcsin(np.sqrt(h))

D = haversine_km(lat[:, None], lon[:, None], lat[None, :], lon[None, :])
np.fill_diagonal(D, np.inf)
A = (D <= NEIGH_KM).astype(float)

def rowstd_sparse(Asub):
    d = Asub.sum(1); d[d == 0] = 1
    return sparse.csr_matrix(Asub / d[:, None])

def zstd(v):
    s = v.std()
    return (v - v.mean()) / s if s > 0 else np.zeros_like(v)

# --- global bivariate Moran's I per quarter ---------------------------------
rows = []
for t in range(1, len(quarters)):
    avail = np.where(~np.isnan(d_act[:, t]))[0]
    if len(avail) < 30:
        continue
    sub = A[np.ix_(avail, avail)]
    keep = sub.sum(1) > 0
    avail = avail[keep]
    W = rowstd_sparse(A[np.ix_(avail, avail)])
    zx = zstd(conf[avail, t])          # own conflict
    zy = zstd(d_act[avail, t])         # own activity change
    m = len(avail)
    Ib = (zx @ (W @ zy)) / m
    perm = np.array([zx @ (W @ RNG.permutation(zy)) / m for _ in range(N_PERM)])
    p = (1 + np.sum(np.abs(perm) >= abs(Ib))) / (N_PERM + 1)
    rows.append({"quarter": quarters[t], "n": m, "biv_moran_i": Ib,
                 "p_perm": p, "total_conf": conf[avail, t].sum()})

biv = pd.DataFrame(rows)
biv.to_csv("spatial_bivariate_by_quarter.csv", index=False)
neg_sig = ((biv.biv_moran_i < 0) & (biv.p_perm < 0.05)).sum()
pos_sig = ((biv.biv_moran_i > 0) & (biv.p_perm < 0.05)).sum()
print("Bivariate Moran's I  (own conflict  x  neighbors' activity change)")
print(f"  quarters tested: {len(biv)};  mean I_B = {biv.biv_moran_i.mean():+.3f}")
print(f"  significant (p<.05): {neg_sig} negative, {pos_sig} positive")
print(f"  -> {'negative dominates: conflict sits amid DECLINING neighbors (spillover)' if neg_sig > pos_sig else 'positive dominates'}\n")

# --- local bivariate LISA for the peak-conflict quarter ----------------------
hot_t = int(conf.sum(0).argmax())
hot_q = quarters[hot_t]
avail = np.where(~np.isnan(d_act[:, hot_t]))[0]
sub = A[np.ix_(avail, avail)]
keep = sub.sum(1) > 0
avail = avail[keep]
W = rowstd_sparse(A[np.ix_(avail, avail)])
zx = zstd(conf[avail, hot_t])
zy = zstd(d_act[avail, hot_t])
lag_y = np.asarray(W @ zy).ravel()
local_Ib = zx * lag_y

# permutation significance (permute the neighbor variable y)
absobs = np.abs(local_Ib)
ge = np.ones(len(avail))
for _ in range(N_PERM):
    li = zx * np.asarray(W @ RNG.permutation(zy)).ravel()
    ge += np.abs(li) >= absobs
local_p = ge / (N_PERM + 1)
sig = local_p < 0.05

def quad(zx_i, lagy_i):
    if zx_i > 0 and lagy_i < 0: return "conf_hi_nbr_down"   # spillover hotspot
    if zx_i > 0 and lagy_i > 0: return "conf_hi_nbr_up"
    if zx_i < 0 and lagy_i < 0: return "conf_lo_nbr_down"
    return "conf_lo_nbr_up"

quads = np.array([quad(zx[i], lag_y[i]) for i in range(len(avail))])
quad_sig = np.where(sig, quads, "ns")

hot = pd.DataFrame({
    "market": [ids[i] for i in avail], "lat": lat[avail], "lon": lon[avail],
    "region": region[avail], "conf": conf[avail, hot_t],
    "d_act": d_act[avail, hot_t], "neighbor_dact_z": lag_y,
    "local_biv_i": local_Ib, "local_p": local_p, "quad_sig": quad_sig,
})
hot.to_csv("spatial_bivariate_hotspot.csv", index=False)
with open("spatial_bivariate_meta.json", "w") as f:
    json.dump({"quarter": hot_q, "neigh_km": NEIGH_KM}, f)

print(f"Local bivariate LISA, {hot_q} (peak conflict); significant cells (p<.05):")
for k, desc in [("conf_hi_nbr_down", "conflict HIGH, neighbors DECLINING (spillover)"),
                ("conf_hi_nbr_up", "conflict HIGH, neighbors rising"),
                ("conf_lo_nbr_down", "conflict low, neighbors declining"),
                ("conf_lo_nbr_up", "conflict low, neighbors rising")]:
    print(f"  {desc:<48} {int(np.sum(quad_sig == k))}")
print(f"  not significant: {int(np.sum(quad_sig == 'ns'))}")
print("\nWrote spatial_bivariate_by_quarter.csv and spatial_bivariate_hotspot.csv")
