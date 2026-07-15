"""
Robustness: repeat the vaccination-conflict relationships using conflict
FATALITIES instead of event COUNTS. Fatalities measure severity and are (a bit)
less mechanically tied to population than raw event counts, so this is a useful
check on the confounding story. Prints events-vs-fatalities side by side.
"""
import re
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.formula.api as smf

VAX = r"C:\Users\Andrew Schroeder\Documents\R_Dir\Ethiopia_Vaccination"
GC = VAX + r"\MeasureDHS-20260715T125833Z-1-001\MeasureDHS\ET_2024_25\extracted\GC\ETGC8AFL.csv"

cl = pd.read_csv("cluster_integration.csv")
zt = pd.read_csv("zone_integration_summary.csv")
gc = pd.read_csv(GC)[["DHSCLUST", "UN_Population_Density_2020"]].rename(
    columns={"DHSCLUST": "V001", "UN_Population_Density_2020": "pop_dens"})
cl = cl.merge(gc, on="V001", how="left")
cl["log_dens"] = np.log1p(cl["pop_dens"])
cl["log_fatal"] = np.log1p(cl["fatal_near"])

def pcorr(y, x, z=None):
    d = pd.DataFrame({"y": y, "x": x});
    if z is not None: d["z"] = z
    d = d.dropna()
    if z is None:
        r, p = stats.pearsonr(d["y"], d["x"]); return r, p, len(d)
    ry = smf.ols("y ~ z", d).fit().resid; rx = smf.ols("x ~ z", d).fit().resid
    r, p = stats.pearsonr(ry, rx); return r, p, len(d)

print("Distribution of nearby fatalities (cluster, 50km):")
print(cl["fatal_near"].describe().round(1).to_string(), "\n")

print("=" * 74)
print("CLUSTER LEVEL (n=722): MCV1 vs conflict, EVENTS vs FATALITIES")
print("=" * 74)
for label, col in [("events", "conflict_near"), ("fatalities", "fatal_near"),
                   ("log(1+fatalities)", "log_fatal")]:
    r, p, n = pcorr(cl["smooth"], cl[col])
    rp, pp, _ = pcorr(cl["smooth"], cl[col], cl["log_dens"])
    print(f"  MCV1 vs {label:<18} raw r={r:+.2f} (p={p:.1e})   "
          f"partial|density r={rp:+.2f} (p={pp:.1e})")

print("\nWeighted OLS of MCV1 (children-weighted), fatalities measure:")
d2 = cl.dropna(subset=["log_dens"]).copy()
for label, formula in [
    ("MCV1 ~ fatal_near", "smooth ~ fatal_near"),
    ("MCV1 ~ fatal_near + log_dens", "smooth ~ fatal_near + log_dens"),
    ("MCV1 ~ fatal_near + log_dens + C(region)", "smooth ~ fatal_near + log_dens + C(region)"),
    ("MCV1 ~ log_fatal + log_dens + C(region)", "smooth ~ log_fatal + log_dens + C(region)"),
]:
    m = smf.wls(formula, data=d2, weights=d2["n_children"]).fit()
    term = "fatal_near" if "fatal_near" in formula and "log_fatal" not in formula else "log_fatal"
    b, p = m.params[term], m.pvalues[term]
    print(f"  {label:<46} {term} beta={b:+.5f} (p={p:.3f})")

print("\n" + "=" * 74)
print("ZONE LEVEL (n=67 zones with reliable MCV1)")
print("=" * 74)
z = zt[zt["mcv1_smooth"].notna()].copy()
for label, col in [("conflict events", "conflict_events"), ("fatalities", "fatalities")]:
    r, p, n = pcorr(z["mcv1_smooth"], z[col])
    rho, _ = stats.spearmanr(z["mcv1_smooth"], z[col])
    print(f"  MCV1 vs {label:<16} r={r:+.2f} (p={p:.1e})  rho={rho:+.2f}  n={n}")

# case-fatality-style: fatalities per event (severity), zone level
z["fatal_per_event"] = np.where(z["conflict_events"] > 0, z["fatalities"] / z["conflict_events"], np.nan)
r, p, n = pcorr(z["mcv1_smooth"], z["fatal_per_event"])
print(f"  MCV1 vs fatalities/event   r={r:+.2f} (p={p:.1e})  n={n}")

print("\nRegional means: events vs fatalities near clusters (are they redundant?)")
reg = cl.groupby("region").agg(
    mcv1=("smooth", lambda s: np.average(s, weights=cl.loc[s.index, "n_children"])),
    events=("conflict_near", "mean"),
    fatalities=("fatal_near", "mean")).sort_values("mcv1")
print((reg.assign(mcv1=(reg["mcv1"]*100).round(0), events=reg["events"].round(0),
                  fatalities=reg["fatalities"].round(0))).to_string())
print(f"\ncorr(mean events, mean fatalities) across regions = "
      f"{reg['events'].corr(reg['fatalities']):.2f}")
reg.to_csv("vax_region_fatalities.csv")
