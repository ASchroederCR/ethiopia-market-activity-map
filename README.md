# Ethiopia Rural Market Activity Index — Interactive Map

**Live map:** https://aschroedercr.github.io/ethiopia-market-activity-map/

A time-enabled Leaflet map of 1,770 satellite-detected rural marketplaces across Ethiopia, built from the
data underlying:

> Doran, Boehnke & Krusell, *"Using satellite imagery to map rural marketplaces and monitor their activity
> at high frequency"* ([arXiv:2407.12953](https://arxiv.org/abs/2407.12953)),
> replication package: https://github.com/pauldingus/MAI-replication-package

Each market is plotted at its detected location and colored/sized by a quarterly **activity index**
(2018Q1–2023Q4), aggregated from ~2.8M per-image satellite observations down to one value per
market per quarter, restricted to images captured on the market's own trading day. An index of ~0
corresponds to a typical non-market-day baseline; ~100 corresponds to typical market-day activity for
that market in that year.

## Features

- Quarter slider with play/pause animation
- Region (admin-1) filter and boundary overlay
- Zoom-dependent hex-bin summary: zoomed out, markets are aggregated into ~25 km geographic hexagons
  colored by the mean activity index of the markets inside; zoom in past level 8 to reveal the
  individual market points. Toggleable, and respects the quarter slider and region filter. Click a hex
  for a popup with a sparkline of the bin's quarterly mean activity index over the full period
- [ACLED](https://acleddata.com/) conflict events plotted as color-coded squares for the selected
  quarter, toggleable, drawn on top of both the hexagons and the market points to visually compare
  conflict incidence against market activity
- Click any market for a popup with a sparkline of its full quarterly activity-index time series
  (highlighting the currently selected quarter) with a second series overlaid showing nearby ACLED
  conflict events per quarter (count within 20km, and fatalities), so activity dips can be checked
  against conflict spikes market by market
- Self-contained `index.html` — the aggregated data is embedded directly in the page, no server-side
  component required

## Files

- `index.html` — the map (also the GitHub Pages entry point)
- `map_template.html` — the HTML/JS template before data is embedded
- `build_ethiopia_data.py` — aggregates the raw per-image CSVs and ACLED conflict data (from the
  replication package) into `ethiopia_market_activity.json` (which also carries each market's nearby-conflict
  series), `ethiopia_conflict_events.json`, and `ethiopia_adm1.geojson`
- `build_html.py` — embeds those JSON/GeoJSON files into `map_template.html` to produce `index.html`
- `ethiopia_market_activity.json`, `ethiopia_conflict_events.json`, `ethiopia_adm1.geojson` — the
  aggregated data, also embedded in `index.html`

## Analysis: does conflict depress market activity?

A set of first-differenced panel regressions relate quarter-over-quarter changes in each market's
activity index to changes in nearby ACLED conflict. First-differencing removes each market's fixed
characteristics (baseline size/visibility); quarter fixed effects absorb nationwide shocks (weather,
seasons, sensor mix); standard errors are clustered by market. Panel: 1,764 markets × 2018Q1–2023Q4,
N = 31,772 market-quarter changes.

![Spatial decay of the conflict effect on market activity](distance_decay.png)

**Headline result:** conflict events are followed by measurable drops in nearby market activity, and the
per-event effect decays steeply with distance — an event within 10 km is associated with roughly a
−0.86-point change in the activity index (where 100 ≈ normal market-day activity) versus −0.29 at
20–50 km. The relationship is significant at every distance and the effect persists (weakly) into the
following quarter.

Scripts and outputs, in order of dependency:

- `analyze_conflict_activity.py` — the main model at a fixed 20 km catchment (pooled, + quarter FE,
  fatalities instead of counts, and a lagged-conflict spec). Writes `conflict_activity_panel.csv`, the
  assembled market-quarter panel (levels + first differences) for reuse in Stata/R.
- `radius_sensitivity.py` → `radius_sensitivity_results.csv` — re-estimates the main spec at 10 / 20 /
  50 km catchments.
- `distance_decay.py` → `distance_decay_rings.csv`, `distance_decay_results.csv` — estimates the spatial
  gradient directly: (A) concentric-ring event counts (0–10 / 10–20 / 20–50 km) in a single joint model,
  and (B) a continuous exponential-decay exposure Σ exp(−distance/λ) at λ = 5/10/20/40 km. Model fit
  (AIC) favors a decay length of ~20–40 km, so a ~20 km catchment captures most of the signal while
  staying attributable to a specific market.
- `plot_distance_decay.py` → `distance_decay.png` — the figure above.
- `plot_distance_decay_by_region.py` → `distance_decay_by_region.png`, `distance_decay_by_region.csv` —
  the same ring model re-estimated separately per admin-1 region, to compare the gradient across areas
  (see below).
- `spacetime_decay.py` → `spacetime_decay_rings.csv`, `spacetime_decay_exp.csv` — adds 1- and 2-quarter
  time lags to the ring and exponential-decay models (distributed-lag), testing whether the effect
  persists over time (see "Time lags" below).
- `plot_spacetime_decay.py` → `spacetime_decay.png` — the space × time figure.

### Regional variation

![Conflict–activity distance decay by region](distance_decay_by_region.png)

Re-estimating the ring model region by region (each panel repeats the national gradient in grey for
reference) shows the national average masks real heterogeneity:

- **Amhara** tracks the national gradient almost exactly, with tight confidence intervals — it is the
  main driver of the countrywide result.
- **Oromia** is noisy at close range (the 0–10 km estimate is even positive but not significant) and only
  becomes precise in the 20–50 km ring.
- **SNNPR** points to a much larger near-market effect (≈ −2 to −3 points per event) but with wide CIs, so
  the magnitude is uncertain.
- **Tigray** is the clear outlier: little or no per-event effect near markets, and a *positive* 20–50 km
  coefficient. The 2018–2023 window spans the Tigray war, when conflict was near-ubiquitous and markets
  were disrupted region-wide, so a market-level event count within a small radius is a poor proxy for
  exposure there — a caution against reading the national coefficient as uniform.

Only regions with ≥150 markets are shown; Beneshangul Gumu (39), Afar (15), Gambela (3) and Somali (3)
have too few markets for reliable market-clustered inference and are omitted.

### Time lags: does the effect persist? (space × time decay)

The distance-decay model above is contemporaneous — conflict in quarter *t* vs the activity change in
*t*. To check whether conflict effects also carry forward in time, `spacetime_decay.py` extends it to a
**finite distributed-lag model in first differences**: each distance ring (and each exponential-decay
exposure) enters at lags 0, 1 and 2 quarters, keeping the quarter-FE, market-clustered design. Writes
`spacetime_decay_rings.csv` and `spacetime_decay_exp.csv`; `plot_spacetime_decay.py` →
`spacetime_decay.png`.

![Space × time decay of the conflict effect](spacetime_decay.png)

Findings:

- **Time lags are jointly significant** (Wald test that all six lag terms = 0: F = 26.7, p = 1.6e-4), so
  the effect is not purely contemporaneous — this validates adding the lag structure.
- **Persistence is concentrated near the market.** At 0–10 km the effect is still significant one quarter
  later (t−1 ≈ −0.70, p < 0.05); the 10–20 km ring is essentially same-quarter only. So near-market
  conflict has both the largest *and* the longest-lived footprint — decay in space and in time together.
- **The cumulative (long-run) effect exceeds the same-quarter effect** for every ring — e.g. 20–50 km:
  −0.31 contemporaneous vs −0.49 summed over three quarters — because the lag terms mostly reinforce
  rather than reverse. The long-run effect is the impact on the *level* of activity of a sustained
  conflict increase.
- **But the lags add little explanatory power** (ΔAIC = −7, R² 0.0347 → 0.0353 on N = 28,945), and the
  best-fitting spatial scale is unchanged (λ ≈ 20 km). Practically: the contemporaneous model is a good
  approximation; lags refine it (and reveal near-market persistence) without overturning the spatial
  conclusion.

These are associations, not causal estimates: conflict may co-move with other local disruptions
(displacement, road closures), and cloud cover can correlate with season/region in ways quarter dummies
only partly absorb.

Note: the analysis scripts read `ethiopia_market_activity.json` from this repo, but the ones that rebuild
conflict exposure by radius (`radius_sensitivity.py`, `distance_decay.py`, `spacetime_decay.py`) also need
the raw ACLED event CSV from the [MAI replication package](https://github.com/pauldingus/MAI-replication-package)
(`datasets/conflict/2012-07-01-2025-07-01-Ethiopia.csv`), expected at `../repo/` relative to the script.
ACLED data is subject to [ACLED's terms of use](https://acleddata.com/terms-of-use/) and is therefore not
redistributed here beyond the derived quarterly aggregates.

This is an independent, derived visualization — not part of the official replication package, and not
affiliated with the paper's authors.
