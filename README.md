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
- Click any market for a popup with a sparkline of its full quarterly time series, highlighting the
  currently selected quarter
- Self-contained `index.html` — the aggregated data is embedded directly in the page, no server-side
  component required

## Files

- `index.html` — the map (also the GitHub Pages entry point)
- `map_template.html` — the HTML/JS template before data is embedded
- `build_ethiopia_data.py` — aggregates the raw per-image CSVs (from the replication package) into
  `ethiopia_market_activity.json` and `ethiopia_adm1.geojson`
- `build_html.py` — embeds those JSON/GeoJSON files into `map_template.html` to produce `index.html`
- `ethiopia_market_activity.json`, `ethiopia_adm1.geojson` — the aggregated data, also embedded in `index.html`

This is an independent, derived visualization — not part of the official replication package, and not
affiliated with the paper's authors.
