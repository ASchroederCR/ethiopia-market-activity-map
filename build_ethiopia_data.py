"""
Aggregate Ethiopia market-activity raw data (per-image observations, 2016-2025)
into a quarterly per-market activity index, and export supporting GeoJSON layers,
for use in a time-enabled Leaflet map.
"""
import glob
import json
import geopandas as gpd
import numpy as np
import pandas as pd

REPO = "../repo"
CONFLICT_RADIUS_KM = 20  # "nearby" catchment radius for market <-> ACLED event linkage

# ---------------------------------------------------------------------------
# 1. Load and aggregate the raw per-image activity observations to quarterly
# ---------------------------------------------------------------------------
NORM_YEARS = range(2018, 2024)  # only years 2018-2023 have a norm_ column

cols_needed = [
    "Location", "date", "mktDay", "weekday", "activity_measure",
    "clear_percent", "marketlat", "marketlon", "admLvl1",
] + [f"activity_measure_norm_{y}" for y in NORM_YEARS]

frames = []
files = sorted(glob.glob(f"{REPO}/datasets/activity_raw/df_ETH_20251216_batch*.csv"))
print(f"Found {len(files)} batch files")

for f in files:
    df = pd.read_csv(f, usecols=lambda c: c in cols_needed)
    df["date"] = pd.to_datetime(df["date"])
    df["year"] = df["date"].dt.year

    # pick the normalized activity value that matches each row's own year
    df["act_idx"] = pd.NA
    for y in NORM_YEARS:
        col = f"activity_measure_norm_{y}"
        if col in df.columns:
            mask = df["year"] == y
            df.loc[mask, "act_idx"] = df.loc[mask, col]
    df["act_idx"] = pd.to_numeric(df["act_idx"], errors="coerce")
    # a handful of locations have a near-zero market-day/non-market-day baseline
    # gap, which blows up the normalized ratio to +/-tens of thousands for a
    # single image; winsorize per-image values before averaging so a couple of
    # extreme images don't dominate a market's whole quarterly mean.
    df["act_idx"] = df["act_idx"].clip(-150, 400)

    df["quarter"] = df["date"].dt.to_period("Q").astype(str)  # e.g. '2019Q3'

    frames.append(df[["Location", "quarter", "mktDay", "weekday", "act_idx",
                       "clear_percent", "marketlat", "marketlon", "admLvl1"]])

raw = pd.concat(frames, ignore_index=True)
print("Total raw observation rows:", len(raw))
print("Unique locations:", raw["Location"].nunique())

# static per-location metadata (lat/lon/region) — take the first non-null value
meta = (raw.dropna(subset=["marketlat", "marketlon"])
           .groupby("Location")
           .agg(lat=("marketlat", "first"),
                lon=("marketlon", "first"),
                admLvl1=("admLvl1", "first"))
           .reset_index())

# market day(s) of week actually observed as active (mktDay == 1)
mkt_days = (raw[raw["mktDay"] == 1]
            .groupby("Location")["weekday"]
            .apply(lambda s: sorted(set(int(x) for x in s.dropna().unique())))
            .reset_index()
            .rename(columns={"weekday": "market_weekdays"}))

# quarterly aggregation of the activity index -- restricted to images captured
# on the market's own designated market day(s). Non-market-day images are the
# ~0 baseline used to *construct* act_idx, not a signal of market activity, so
# mixing them in would dilute/distort the index.
quarterly = (raw[raw["mktDay"] == 1].dropna(subset=["act_idx"])
                .groupby(["Location", "quarter"])
                .agg(act_idx=("act_idx", "mean"), n_obs=("act_idx", "size"))
                .reset_index())

print("Quarterly rows:", len(quarterly))
print("Quarter range:", quarterly["quarter"].min(), "-", quarterly["quarter"].max())

# ---------------------------------------------------------------------------
# 2. Build compact per-market JSON: {loc: {lat, lon, region, market_weekdays, q: {quarter: [idx, n]}}}
# ---------------------------------------------------------------------------
meta = meta.merge(mkt_days, on="Location", how="left")
meta["market_weekdays"] = meta["market_weekdays"].apply(
    lambda v: v if isinstance(v, list) else [])

markets = {}
for row in meta.itertuples():
    markets[row.Location] = {
        "lat": round(row.lat, 5),
        "lon": round(row.lon, 5),
        "region": row.admLvl1 if isinstance(row.admLvl1, str) else "Unknown",
        "mktDays": row.market_weekdays,
        "q": {},
    }

for row in quarterly.itertuples():
    if row.Location in markets:
        markets[row.Location]["q"][row.quarter] = [
            round(float(row.act_idx), 1), int(row.n_obs)
        ]

all_quarters = sorted(quarterly["quarter"].unique().tolist())

# ---------------------------------------------------------------------------
# 2b. ACLED conflict events: (a) a standalone point layer, quarter-filtered,
# and (b) a per-market "nearby conflict" quarterly series (spatial join on a
# CONFLICT_RADIUS_KM catchment) merged into each market's own record so it can
# be plotted alongside that market's activity-index sparkline.
# ---------------------------------------------------------------------------
conflict = pd.read_csv(f"{REPO}/datasets/conflict/2012-07-01-2025-07-01-Ethiopia.csv")
conflict["event_date"] = pd.to_datetime(conflict["event_date"], format="mixed")
conflict["quarter"] = conflict["event_date"].dt.to_period("Q").astype(str)
# align to the same quarter range as the market activity data, so both layers
# are driven by the same time slider
conflict = conflict[conflict["quarter"].isin(all_quarters)].reset_index(drop=True)
print("Conflict events in range:", len(conflict))

event_types = sorted(conflict["event_type"].dropna().unique().tolist())

def trim(s, n=220):
    if not isinstance(s, str):
        return ""
    return s if len(s) <= n else s[:n].rsplit(" ", 1)[0] + "…"

events_out = []
for row in conflict.itertuples():
    events_out.append({
        "lat": round(row.latitude, 4),
        "lon": round(row.longitude, 4),
        "q": row.quarter,
        "type": row.event_type,
        "sub": row.sub_event_type,
        "fat": int(row.fatalities) if pd.notna(row.fatalities) else 0,
        "date": row.event_date.strftime("%Y-%m-%d"),
        "loc": row.location,
        "a1": row.actor1 if isinstance(row.actor1, str) else "",
        "notes": trim(row.notes, 140),
    })

with open("ethiopia_conflict_events.json", "w", encoding="utf-8") as f:
    json.dump({"event_types": event_types, "events": events_out}, f,
               separators=(",", ":"), ensure_ascii=False)
print("Wrote ethiopia_conflict_events.json with", len(events_out), "events across",
      len(event_types), "event types")

# spatial join: haversine distance between every market and every event
def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0088
    lat1r, lon1r, lat2r, lon2r = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat, dlon = lat2r - lat1r, lon2r - lon1r
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1r) * np.cos(lat2r) * np.sin(dlon / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))

mkt_lat = meta["lat"].to_numpy()
mkt_lon = meta["lon"].to_numpy()
ev_lat = conflict["latitude"].to_numpy()
ev_lon = conflict["longitude"].to_numpy()
ev_quarter = conflict["quarter"].to_numpy()
ev_fatal = conflict["fatalities"].fillna(0).to_numpy()

dist_km = haversine_km(mkt_lat[:, None], mkt_lon[:, None], ev_lat[None, :], ev_lon[None, :])
within = dist_km <= CONFLICT_RADIUS_KM

n_with_nearby_conflict = 0
for i, loc in enumerate(meta["Location"]):
    idxs = np.nonzero(within[i])[0]
    if len(idxs) == 0 or loc not in markets:
        continue
    n_with_nearby_conflict += 1
    sub_q = ev_quarter[idxs]
    sub_f = ev_fatal[idxs]
    qdict = {}
    for q in all_quarters:
        m = sub_q == q
        cnt = int(m.sum())
        if cnt:
            qdict[q] = [cnt, int(sub_f[m].sum())]
    markets[loc]["conflict"] = qdict

print(f"{n_with_nearby_conflict}/{len(markets)} markets have >=1 ACLED event within "
      f"{CONFLICT_RADIUS_KM}km at some point in the period")

# ---------------------------------------------------------------------------
# 2c. Write the combined per-market activity+conflict JSON
# ---------------------------------------------------------------------------
output = {
    "quarters": all_quarters,
    "conflict_radius_km": CONFLICT_RADIUS_KM,
    "markets": markets,
}

with open("ethiopia_market_activity.json", "w") as f:
    json.dump(output, f, separators=(",", ":"))

print("Wrote ethiopia_market_activity.json with", len(markets), "markets and",
      len(all_quarters), "quarters")

# ---------------------------------------------------------------------------
# 3. Export admin-1 boundaries as GeoJSON (simplified) for basemap context
# ---------------------------------------------------------------------------
adm1 = gpd.read_file(f"{REPO}/datasets/shapefiles/Eth_Adm1.shp")
adm1 = adm1.to_crs(4326)
adm1["geometry"] = adm1.geometry.simplify(0.01, preserve_topology=True)
adm1.to_file("ethiopia_adm1.geojson", driver="GeoJSON")
print("Wrote ethiopia_adm1.geojson with", len(adm1), "features")

# ---------------------------------------------------------------------------
# 4. Export detected-market polygons (for reference / non-time layer) as GeoJSON
# ---------------------------------------------------------------------------
mkts = gpd.read_file(f"{REPO}/datasets/shapefiles/eth_detectedMarkets.shp")
mkts = mkts.to_crs(4326)
mkts["marketDays"] = mkts["marketDays"].astype(str)
mkts.to_file("ethiopia_detected_markets.geojson", driver="GeoJSON")
print("Wrote ethiopia_detected_markets.geojson with", len(mkts), "features")
