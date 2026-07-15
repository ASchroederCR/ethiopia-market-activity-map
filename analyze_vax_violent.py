"""
Re-run the vaccination-conflict relationships restricting ACLED to ARMED
VIOLENCE only -- excluding 'Protests', 'Riots', and 'Strategic developments',
which are heavily urban and mostly non-lethal and are a prime suspect for the
population/urbanicity confound. Keeps: Battles, Violence against civilians,
Explosions/Remote violence. Recomputes cluster- and zone-level conflict
covariates from scratch and compares to the all-events results.
"""
import json
import numpy as np, pandas as pd, geopandas as gpd
from scipy import stats
import statsmodels.formula.api as smf

VAX = r"C:\Users\Andrew Schroeder\Documents\R_Dir\Ethiopia_Vaccination"
GC = VAX + r"\MeasureDHS-20260715T125833Z-1-001\MeasureDHS\ET_2024_25\extracted\GC\ETGC8AFL.csv"
KEEP = {"Battles", "Violence against civilians", "Explosions/Remote violence"}
DROP = {"Protests", "Riots", "Strategic developments"}
NEIGH_KM = 50

cl = pd.read_csv("cluster_integration.csv")   # V001, lat, lon, smooth, n_children, region, conflict_near(all), fatal_near(all)
conf = json.load(open("ethiopia_conflict_events.json", encoding="utf-8"))
ev = pd.DataFrame(conf["events"])
print("ACLED event mix (2018-2023):")
print(ev["type"].value_counts().to_string())
evv = ev[ev["type"].isin(KEEP)].reset_index(drop=True)
print(f"\nKept {len(evv)} of {len(ev)} events ({len(evv)/len(ev)*100:.0f}%) as armed violence; "
      f"dropped {len(ev)-len(evv)} protests/riots/strategic.\n")

def hav(lat1, lon1, lat2, lon2):
    R = 6371.0088
    a1, o1, a2, o2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat, dlon = a2 - a1, o2 - o1
    h = np.sin(dlat/2)**2 + np.cos(a1)*np.cos(a2)*np.sin(dlon/2)**2
    return 2 * R * np.arcsin(np.sqrt(h))

# recompute cluster neighbourhood conflict using violent events only
clat, clon = cl["lat"].to_numpy(), cl["lon"].to_numpy()
De = hav(clat[:, None], clon[:, None], evv["lat"].to_numpy()[None, :], evv["lon"].to_numpy()[None, :])
near = De <= NEIGH_KM
cl["viol_near"] = near.sum(1)
cl["viol_fatal_near"] = [evv["fat"].to_numpy()[near[i]].sum() if near[i].any() else 0 for i in range(len(cl))]

gc = pd.read_csv(GC)[["DHSCLUST", "UN_Population_Density_2020"]].rename(
    columns={"DHSCLUST": "V001", "UN_Population_Density_2020": "pop_dens"})
cl = cl.merge(gc, on="V001", how="left")
cl["log_dens"] = np.log1p(cl["pop_dens"])

def pcorr(y, x, z=None):
    d = pd.DataFrame({"y": y, "x": x})
    if z is not None: d["z"] = z
    d = d.dropna()
    if z is None:
        r, p = stats.pearsonr(d["y"], d["x"]); return r, p, len(d)
    ry = smf.ols("y ~ z", d).fit().resid; rx = smf.ols("x ~ z", d).fit().resid
    r, p = stats.pearsonr(ry, rx); return r, p, len(d)

print("=" * 76)
print("CLUSTER LEVEL (n=722): all-events vs armed-violence-only")
print("=" * 76)
rows = [("events: ALL types", "conflict_near"), ("events: violence only", "viol_near"),
        ("fatalities: ALL types", "fatal_near"), ("fatalities: violence only", "viol_fatal_near")]
for label, col in rows:
    r, p, n = pcorr(cl["smooth"], cl[col])
    rp, pp, _ = pcorr(cl["smooth"], cl[col], cl["log_dens"])
    print(f"  MCV1 vs {label:<26} raw r={r:+.2f} (p={p:.1e})   partial|density r={rp:+.2f} (p={pp:.1e})")

print("\nWeighted OLS (children-weighted), armed-violence event count:")
d2 = cl.dropna(subset=["log_dens"])
for label, formula, term in [
    ("MCV1 ~ viol_near", "smooth ~ viol_near", "viol_near"),
    ("MCV1 ~ viol_near + log_dens", "smooth ~ viol_near + log_dens", "viol_near"),
    ("MCV1 ~ viol_near + log_dens + C(region)", "smooth ~ viol_near + log_dens + C(region)", "viol_near"),
]:
    m = smf.wls(formula, data=d2, weights=d2["n_children"]).fit()
    print(f"  {label:<46} beta={m.params[term]:+.5f} (p={m.pvalues[term]:.3f})")

# ---- zone level: recompute violent conflict per zone ----------------------
print("\n" + "=" * 76)
print("ZONE LEVEL (n=67 zones with reliable MCV1)")
print("=" * 76)
zt = gpd.read_file("ethiopia_zones_vax.geojson")[["zone", "mcv1_smooth", "mcv1_reliable", "geometry"]]
evg = gpd.GeoDataFrame(evv, geometry=gpd.points_from_xy(evv["lon"], evv["lat"]), crs=4326)
evz = gpd.sjoin(evg, zt, predicate="within")
viol_by_zone = evz.groupby("zone").agg(viol_events=("fat", "size"), viol_fatal=("fat", "sum")).reset_index()
z = zt.drop(columns="geometry").merge(viol_by_zone, on="zone", how="left")
z[["viol_events", "viol_fatal"]] = z[["viol_events", "viol_fatal"]].fillna(0)
z = z[z["mcv1_smooth"].notna()]
for label, col in [("armed-violence events", "viol_events"), ("armed-violence fatalities", "viol_fatal")]:
    r, p, n = pcorr(z["mcv1_smooth"], z[col])
    rho, _ = stats.spearmanr(z["mcv1_smooth"], z[col])
    print(f"  MCV1 vs {label:<26} r={r:+.2f} (p={p:.1e})  rho={rho:+.2f}  n={n}")

print("\nRegional means: all events vs armed-violence-only near clusters")
reg = cl.groupby("region").agg(
    mcv1=("smooth", lambda s: np.average(s, weights=cl.loc[s.index, "n_children"])),
    all_events=("conflict_near", "mean"),
    violence=("viol_near", "mean")).sort_values("mcv1")
reg["pct_violent"] = (reg["violence"] / reg["all_events"] * 100).round(0)
print((reg.assign(mcv1=(reg["mcv1"]*100).round(0), all_events=reg["all_events"].round(0),
                  violence=reg["violence"].round(0))).to_string())
reg.to_csv("vax_region_violent.csv")
print("\nWrote vax_region_violent.csv")
