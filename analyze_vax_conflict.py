"""
Do low vaccination, conflict, and market decline coincide spatially?

Tests the association at two levels, and — critically — checks whether any
apparent link survives adjustment for the dominant confounder: Ethiopia's
highland/lowland (agrarian/pastoralist) development gradient, which drives both
where markets exist and where vaccination is low, largely independently of the
2018-23 conflict (which concentrated in the northern highlands).
"""
import json
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.formula.api as smf

cl = pd.read_csv("cluster_integration.csv")
zt = pd.read_csv("zone_integration_summary.csv")

def corr_line(name, x, y, w=None):
    m = (~pd.isna(x)) & (~pd.isna(y))
    x, y = np.asarray(x)[m], np.asarray(y)[m]
    r, p = stats.pearsonr(x, y)
    rho, _ = stats.spearmanr(x, y)
    print(f"  {name:<44} n={m.sum():4d}  r={r:+.2f} (p={p:.1e})  rho={rho:+.2f}")

print("=" * 78)
print("CLUSTER LEVEL (n=722 DHS clusters; MCV1 = smoothed posterior mean)")
print("=" * 78)
print("Neighbourhood covariates within 50 km of each cluster:")
corr_line("MCV1  vs  nearby conflict events", cl["smooth"], cl["conflict_near"])
corr_line("MCV1  vs  nearby conflict fatalities", cl["smooth"], cl["fatal_near"])
corr_line("MCV1  vs  nearby market activity (level)", cl["smooth"], cl["mkt_idx_near"])
corr_line("MCV1  vs  nearby market activity (trend)", cl["smooth"], cl["mkt_trend_near"])
corr_line("MCV1  vs  has any market within 50km", cl["smooth"], (cl["n_markets_near"] > 0).astype(float))

print("\nAmong clusters WITH >=1 nearby market (n={}):".format((cl["n_markets_near"] > 0).sum()))
sub = cl[cl["n_markets_near"] > 0]
corr_line("MCV1  vs  nearby market activity (level)", sub["smooth"], sub["mkt_idx_near"])
corr_line("MCV1  vs  nearby market activity (trend)", sub["smooth"], sub["mkt_trend_near"])
corr_line("MCV1  vs  nearby conflict events", sub["smooth"], sub["conflict_near"])

print("\n" + "=" * 78)
print("ZONE LEVEL (admin-2; MCV1 child-weighted, zones with >=3 clusters)")
print("=" * 78)
z = zt[zt["mcv1_smooth"].notna()].copy()
corr_line("MCV1  vs  conflict events", z["mcv1_smooth"], z["conflict_events"])
corr_line("MCV1  vs  fatalities", z["mcv1_smooth"], z["fatalities"])
corr_line("MCV1  vs  market activity (mean index)", z["mcv1_smooth"], z["mkt_idx"])
corr_line("MCV1  vs  market activity (trend)", z["mcv1_smooth"], z["mkt_trend"])
corr_line("MCV1  vs  number of markets", z["mcv1_smooth"], z["n_markets"])

print("\nCross-tab: mean MCV1 by conflict exposure and market presence (zones)")
z["hi_conflict"] = z["conflict_events"] > z["conflict_events"].median()
z["has_market"] = z["n_markets"] > 0
tab = z.groupby(["has_market", "hi_conflict"])["mcv1_smooth"].agg(["mean", "size"])
print((tab.assign(mean=(tab["mean"]*100).round(0))).to_string())

print("\n" + "=" * 78)
print("CONFOUNDING CHECK: does conflict predict low MCV1 within regions?")
print("=" * 78)
# region FE at cluster level, weighted by children; use raw coverage (honest)
cl2 = cl.dropna(subset=["mkt_trend_near"]).copy()
m_uncond = smf.wls("smooth ~ conflict_near", data=cl, weights=cl["n_children"]).fit()
m_region = smf.wls("smooth ~ conflict_near + C(region)", data=cl, weights=cl["n_children"]).fit()
print(f"  MCV1 ~ conflict_near               : beta={m_uncond.params['conflict_near']:+.5f} "
      f"(p={m_uncond.pvalues['conflict_near']:.3f})")
print(f"  MCV1 ~ conflict_near + region FE    : beta={m_region.params['conflict_near']:+.5f} "
      f"(p={m_region.pvalues['conflict_near']:.3f})")
mt_uncond = smf.wls("smooth ~ mkt_trend_near", data=cl2, weights=cl2["n_children"]).fit()
mt_region = smf.wls("smooth ~ mkt_trend_near + C(region)", data=cl2, weights=cl2["n_children"]).fit()
print(f"  MCV1 ~ market_trend                 : beta={mt_uncond.params['mkt_trend_near']:+.4f} "
      f"(p={mt_uncond.pvalues['mkt_trend_near']:.3f})")
print(f"  MCV1 ~ market_trend + region FE      : beta={mt_region.params['mkt_trend_near']:+.4f} "
      f"(p={mt_region.pvalues['mkt_trend_near']:.3f})")

print("\nRegional means (are conflict & low-vax in the same regions?)")
reg = cl.groupby("region").agg(
    mcv1=("smooth", lambda s: np.average(s, weights=cl.loc[s.index, "n_children"])),
    conflict_near=("conflict_near", "mean"),
    mkt_trend=("mkt_trend_near", "mean"),
    n_clusters=("V001", "size")).sort_values("mcv1")
print((reg.assign(mcv1=(reg["mcv1"]*100).round(0),
                  conflict_near=reg["conflict_near"].round(0),
                  mkt_trend=reg["mkt_trend"].round(2))).to_string())

reg.to_csv("vax_region_summary.csv")
z.to_csv("vax_zone_analysis.csv", index=False)
print("\nWrote vax_region_summary.csv, vax_zone_analysis.csv")
