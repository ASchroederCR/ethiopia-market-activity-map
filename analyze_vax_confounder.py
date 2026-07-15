"""
Is the (surprising) POSITIVE conflict-vaccination correlation just confounding
by population density / urbanicity? ACLED counts events, not events per capita,
so populous, urbanised areas rack up more recorded conflict AND have better
health-system access (higher MCV1). This tests whether the conflict-MCV1 link
survives adjustment for a population-density / urbanicity proxy from the DHS
geographic covariates extract.
"""
import re, json
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.formula.api as smf

VAX = r"C:\Users\Andrew Schroeder\Documents\R_Dir\Ethiopia_Vaccination"
GC = VAX + r"\MeasureDHS-20260715T125833Z-1-001\MeasureDHS\ET_2024_25\extracted\GC\ETGC8AFL.csv"

cl = pd.read_csv("cluster_integration.csv")
gc = pd.read_csv(GC)[["DHSCLUST", "UN_Population_Density_2020", "All_Population_Count_2020",
                      "Nightlights_Composite", "Global_Human_Footprint"]]
gc = gc.rename(columns={"DHSCLUST": "V001", "UN_Population_Density_2020": "pop_dens",
                        "Nightlights_Composite": "nightlights",
                        "Global_Human_Footprint": "footprint"})
d = cl.merge(gc, on="V001", how="left")
d["log_dens"] = np.log1p(d["pop_dens"])

def pr(name, x, y):
    m = (~pd.isna(x)) & (~pd.isna(y))
    r, p = stats.pearsonr(np.asarray(x)[m], np.asarray(y)[m])
    print(f"  {name:<46} r={r:+.2f} (p={p:.1e})")

print("Urbanicity proxies vs conflict and vs MCV1 (cluster level):")
pr("conflict_near   vs  population density", d["conflict_near"], d["pop_dens"])
pr("conflict_near   vs  nightlights", d["conflict_near"], d["nightlights"])
pr("MCV1 (smooth)   vs  population density", d["smooth"], d["pop_dens"])
pr("MCV1 (smooth)   vs  nightlights", d["smooth"], d["nightlights"])

print("\nWeighted OLS of MCV1 (raw idea: does conflict survive density control?)")
d2 = d.dropna(subset=["log_dens"]).copy()
for label, formula in [
    ("MCV1 ~ conflict_near", "smooth ~ conflict_near"),
    ("MCV1 ~ conflict_near + log_dens", "smooth ~ conflict_near + log_dens"),
    ("MCV1 ~ conflict_near + log_dens + nightlights", "smooth ~ conflict_near + log_dens + nightlights"),
    ("MCV1 ~ conflict_near + log_dens + C(region)", "smooth ~ conflict_near + log_dens + C(region)"),
]:
    m = smf.wls(formula, data=d2, weights=d2["n_children"]).fit()
    b, p = m.params["conflict_near"], m.pvalues["conflict_near"]
    print(f"  {label:<52} conflict beta={b:+.5f} (p={p:.3f})")

# partial correlation of MCV1 and conflict controlling for log density
def partial_r(y, x, z):
    dd = pd.DataFrame({"y": y, "x": x, "z": z}).dropna()
    ry = smf.ols("y ~ z", dd).fit().resid
    rx = smf.ols("x ~ z", dd).fit().resid
    r, p = stats.pearsonr(ry, rx)
    return r, p, len(dd)

r, p, n = partial_r(d["smooth"], d["conflict_near"], d["log_dens"])
print(f"\nPartial corr(MCV1, conflict | log population density): r={r:+.2f} (p={p:.1e}, n={n})")
print("  (compare with the +0.57 raw correlation)")
