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
- [ACLED](https://acleddata.com/) conflict events plotted as color-coded triangles for the selected
  quarter, toggleable, to visually compare conflict incidence against market activity
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

## Analysis: distance decay of the conflict–activity relationship

- `distance_decay.py` — first-differenced panel regressions (quarter fixed effects, SEs clustered by
  market) of quarter-over-quarter changes in each market's activity index on changes in nearby ACLED
  conflict, in two forms: (A) concentric-ring event counts (0–10 / 10–20 / 20–50 km) estimated jointly,
  and (B) a continuous exponential-decay exposure Σ exp(−distance/λ) at λ = 5/10/20/40 km.
- `distance_decay_rings.csv`, `distance_decay_results.csv` — the estimates. Headline findings: the
  per-event effect on market activity is negative and significant at all distances and decays steeply
  (~−0.86 points per event within 10 km vs ~−0.29 at 20–50 km); model fit favors a decay length of
  roughly 20–40 km, so a ~20 km catchment captures most of the signal while staying attributable to a
  specific market.

Note: `distance_decay.py` reads `ethiopia_market_activity.json` from this repo, but also needs the raw
ACLED event CSV from the [MAI replication package](https://github.com/pauldingus/MAI-replication-package)
(`datasets/conflict/2012-07-01-2025-07-01-Ethiopia.csv`), expected at `../repo/` relative to the script.
ACLED data is subject to [ACLED's terms of use](https://acleddata.com/terms-of-use/) and is therefore not
redistributed here beyond the derived quarterly aggregates.

This is an independent, derived visualization — not part of the official replication package, and not
affiliated with the paper's authors.
