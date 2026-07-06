"""
Aggregate Ethiopia market-activity raw data (per-image observations, 2016-2025)
into a quarterly per-market activity index, and export supporting GeoJSON layers,
for use in a time-enabled Leaflet map.
"""
import glob
import json
import geopandas as gpd
import pandas as pd

REPO = "../repo"

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

output = {
    "quarters": all_quarters,
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
