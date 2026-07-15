"""
Integrate three layers on a common geography to study whether low measles
vaccination (MCV1) coincides with conflict and with market-activity decline.

Data provenance (all validated / reproduced, see checks printed below):
  - Vaccination: EDHS 2024-25 MCV1, children 12-23 mo. Cluster-level SMOOTHED
    coverage (BYM2/SPDE posterior mean) taken from the report's rendered output;
    RAW cluster coverage independently reproduced here from the DHS microdata
    (KR .DAT via the .DCF dictionary). DHS cluster points are displaced and are
    NOT published; only zone-level aggregates are exported for the public map.
  - Market activity: quarterly MAI index per market, 2018Q1-2023Q4 (this repo).
  - Conflict: ACLED events 2018-2023 aggregated earlier (this repo).

Outputs:
  - ethiopia_zones_vax.geojson   : admin-2 zones with child-weighted MCV1 +
    market + conflict summaries (aggregate; safe to publish on the map).
  - zone_integration_summary.csv : the same zone table, flat.
  - cluster_integration.csv      : cluster-level MCV1 with neighbourhood market
    and conflict covariates (INTERNAL analysis file; not for publication).
"""
import re, json, warnings
import numpy as np, pandas as pd, geopandas as gpd
warnings.filterwarnings("ignore", category=UserWarning)

VAX_DIR = r"C:\Users\Andrew Schroeder\Documents\R_Dir\Ethiopia_Vaccination"
KR = VAX_DIR + r"\MeasureDHS-20260715T125833Z-1-001\MeasureDHS\ET_2024_25\extracted\KR"
GE = VAX_DIR + r"\MeasureDHS-20260715T125833Z-1-001\MeasureDHS\ET_2024_25\extracted\GE\ETGE8AFL.shp"
ADM2 = r"..\repo\datasets\shapefiles\ethiopia_adm2\eth_adm2.shp"
NEIGH_KM = 50   # neighbourhood radius for cluster-level market/conflict covariates

# ---------------------------------------------------------------------------
# 1. Reproduce raw cluster MCV1 from the DHS microdata (validation + weights)
# ---------------------------------------------------------------------------
dcf = open(KR + r"\ETKR8AFL.DCF", encoding="utf-8", errors="replace").read()
items = {}
for blk in dcf.split("[Item]")[1:]:
    n = re.search(r"Name=(\S+)", blk); s = re.search(r"Start=(\d+)", blk); l = re.search(r"Len=(\d+)", blk)
    if n and s and l:
        items[n.group(1).upper()] = (int(s.group(1)), int(l.group(1)))
need = ["V001", "V005", "V024", "V025", "B4", "B5", "B19", "H9"]
def num(x):
    x = x.strip(); return np.nan if x == "" else float(x)
rows = []
with open(KR + r"\ETKR8AFL.DAT", encoding="latin-1") as f:
    for line in f:
        rows.append({v: line[items[v][0]-1:items[v][0]-1+items[v][1]] for v in need})
kr = pd.DataFrame(rows)
for v in need:
    kr[v] = kr[v].map(num)
kr["wt"] = kr["V005"] / 1e6
kr["mcv1"] = np.where(kr["H9"].isin([1, 2, 3]), 1.0, np.where(kr["H9"] == 0, 0.0, np.nan))
samp = kr[(kr["B5"] == 1) & (kr["B19"].between(12, 23)) & (kr["mcv1"].notna())].copy()
nat = np.average(samp["mcv1"], weights=samp["wt"])
print(f"[validate] sample={len(samp)} children, clusters={samp['V001'].nunique()}, "
      f"national weighted MCV1={nat*100:.1f}% (report ~52%, official ~51%)")

raw_cl = samp.groupby("V001").apply(lambda g: pd.Series({
    "n_children": len(g), "n_vacc": int((g["mcv1"] == 1).sum()),
    "raw_wt": np.average(g["mcv1"], weights=g["wt"])})).reset_index()

# ---------------------------------------------------------------------------
# 2. Cluster coordinates (GPS shapefile) + smoothed coverage (report output)
# ---------------------------------------------------------------------------
gps = gpd.read_file(GE)[["DHSCLUST", "LATNUM", "LONGNUM", "URBAN_RURA"]]
gps = gps.rename(columns={"DHSCLUST": "V001", "LATNUM": "lat", "LONGNUM": "lon"})
gps = gps[(gps["lat"] != 0) | (gps["lon"] != 0)]

html = open(VAX_DIR + r"\ethiopia_measles_vaccination_report.html", encoding="utf-8").read()
blk = re.findall(r'<script type="application/json"[^>]*>(.*?)</script>', html, re.S)[1]
cm = [c for c in json.loads(blk)["x"]["calls"] if c["method"] == "addCircleMarkers"][0]
sm = pd.DataFrame([{
    "V001": int(re.search(r"Cluster:</b> (\d+)", ph).group(1)),
    "region": re.search(r"Region:</b> ([^<]+)", ph).group(1),
    "smooth": float(re.search(r"Smoothed MCV1 coverage:</b> (\d+)%", ph).group(1)) / 100,
} for ph in cm["args"][8]])

cl = (raw_cl.merge(gps, on="V001").merge(sm, on="V001"))
print(f"[clusters] {len(cl)} with coords + raw + smoothed coverage")
# cross-check raw vs smoothed agreement
print(f"[validate] corr(raw cluster, smoothed cluster) = "
      f"{np.corrcoef(cl['raw_wt'], cl['smooth'])[0,1]:.2f}")

# ---------------------------------------------------------------------------
# 3. Market activity metrics (level + trend), and conflict events
# ---------------------------------------------------------------------------
act = json.load(open("ethiopia_market_activity.json", encoding="utf-8"))
quarters = act["quarters"]; qi = {q: i for i, q in enumerate(quarters)}
mk = []
for mid, m in act["markets"].items():
    qs = [(qi[q], v[0]) for q, v in m["q"].items()]
    if not qs:
        continue
    xs = np.array([a for a, _ in qs]); ys = np.array([b for _, b in qs])
    trend = np.polyfit(xs, ys, 1)[0] if len(qs) >= 6 else np.nan
    mk.append({"lat": m["lat"], "lon": m["lon"], "mean_idx": ys.mean(),
               "trend": trend, "n_q": len(qs)})
mk = pd.DataFrame(mk)
print(f"[markets] {len(mk)} markets with activity; mean index={mk['mean_idx'].mean():.1f}")

conf = json.load(open("ethiopia_conflict_events.json", encoding="utf-8"))
ev = pd.DataFrame(conf["events"])
print(f"[conflict] {len(ev)} ACLED events 2018-2023")

# ---------------------------------------------------------------------------
# 4. Cluster-level neighbourhood covariates (haversine within NEIGH_KM)
# ---------------------------------------------------------------------------
def hav(lat1, lon1, lat2, lon2):
    R = 6371.0088
    a1, o1, a2, o2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat, dlon = a2 - a1, o2 - o1
    h = np.sin(dlat/2)**2 + np.cos(a1)*np.cos(a2)*np.sin(dlon/2)**2
    return 2 * R * np.arcsin(np.sqrt(h))

clat = cl["lat"].to_numpy(); clon = cl["lon"].to_numpy()
Dm = hav(clat[:, None], clon[:, None], mk["lat"].to_numpy()[None, :], mk["lon"].to_numpy()[None, :])
De = hav(clat[:, None], clon[:, None], ev["lat"].to_numpy()[None, :], ev["lon"].to_numpy()[None, :])
mk_idx = mk["mean_idx"].to_numpy(); mk_tr = mk["trend"].to_numpy()
near_mk = Dm <= NEIGH_KM
near_ev = De <= NEIGH_KM
cl["n_markets_near"] = near_mk.sum(1)
cl["mkt_idx_near"] = [mk_idx[near_mk[i]].mean() if near_mk[i].any() else np.nan for i in range(len(cl))]
cl["mkt_trend_near"] = [np.nanmean(mk_tr[near_mk[i]]) if near_mk[i].any() else np.nan for i in range(len(cl))]
cl["conflict_near"] = near_ev.sum(1)
cl["fatal_near"] = [ev["fat"].to_numpy()[near_ev[i]].sum() if near_ev[i].any() else 0 for i in range(len(cl))]
cl.to_csv("cluster_integration.csv", index=False)
print(f"[clusters] neighbourhood covariates written (radius {NEIGH_KM} km)")

# ---------------------------------------------------------------------------
# 5. Zone-level aggregation (admin-2), safe-to-publish aggregate
# ---------------------------------------------------------------------------
zones = gpd.read_file(ADM2).to_crs(4326)[["shapeName", "geometry"]].rename(columns={"shapeName": "zone"})

def to_gdf(df):
    return gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df["lon"], df["lat"]), crs=4326)

cl_z = gpd.sjoin(to_gdf(cl), zones, predicate="within").drop(columns="index_right")
mk_z = gpd.sjoin(to_gdf(mk), zones, predicate="within").drop(columns="index_right")
ev_z = gpd.sjoin(to_gdf(ev), zones, predicate="within").drop(columns="index_right")

def wmean(v, w):
    v = np.asarray(v, float); w = np.asarray(w, float)
    ok = ~np.isnan(v)
    return np.average(v[ok], weights=w[ok]) if ok.any() and w[ok].sum() > 0 else np.nan

zrows = []
for zname, g in cl_z.groupby("zone"):
    zrows.append({"zone": zname, "n_clusters": len(g),
                  "n_children": int(g["n_children"].sum()),
                  "mcv1_smooth": wmean(g["smooth"], g["n_children"]),
                  "mcv1_raw": wmean(g["raw_wt"], g["n_children"])})
zc = pd.DataFrame(zrows)
zm = mk_z.groupby("zone").agg(n_markets=("mean_idx", "size"),
                              mkt_idx=("mean_idx", "mean"),
                              mkt_trend=("trend", "mean")).reset_index()
ze = ev_z.groupby("zone").agg(conflict_events=("fat", "size"),
                              fatalities=("fat", "sum")).reset_index()

zt = zones.merge(zc, on="zone", how="left").merge(zm, on="zone", how="left").merge(ze, on="zone", how="left")
for c in ["conflict_events", "fatalities", "n_markets"]:
    zt[c] = zt[c].fillna(0)
# suppress MCV1 for sparsely-sampled zones (aggregate-disclosure + noise)
zt["mcv1_reliable"] = zt["n_clusters"] >= 3
zt.loc[~zt["mcv1_reliable"].fillna(False), ["mcv1_smooth", "mcv1_raw"]] = np.nan

zt.drop(columns="geometry").to_csv("zone_integration_summary.csv", index=False)
zt["geometry"] = zt.geometry.simplify(0.01, preserve_topology=True)
zt.to_file("ethiopia_zones_vax.geojson", driver="GeoJSON")
nz = int(zt["mcv1_smooth"].notna().sum())
print(f"[zones] {len(zt)} admin-2 zones; {nz} with reliable MCV1 (>=3 clusters)")
print(f"[zones] {int((zt['n_markets']>0).sum())} have markets, "
      f"{int((zt['conflict_events']>0).sum())} have conflict events")

# quick validation vs report's stated zone extremes
low = zt[zt["mcv1_reliable"]].nsmallest(6, "mcv1_smooth")[["zone", "mcv1_smooth", "n_clusters"]]
print("\n[validate] lowest-MCV1 zones (report said Somali zones lowest, <10%):")
print(low.assign(mcv1=(low["mcv1_smooth"]*100).round(0)).to_string(index=False))
