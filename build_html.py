import json

with open("ethiopia_market_activity.json", encoding="utf-8") as f:
    activity_json = f.read()
with open("ethiopia_adm1.geojson", encoding="utf-8") as f:
    adm1_json = f.read()
with open("ethiopia_conflict_events.json", encoding="utf-8") as f:
    conflict_json = f.read()
with open("ethiopia_zones_vax.geojson", encoding="utf-8") as f:
    vax_json = f.read()

html_template = open("map_template.html", encoding="utf-8").read()
html = (html_template
        .replace("__ACTIVITY_DATA__", activity_json)
        .replace("__ADM1_DATA__", adm1_json)
        .replace("__CONFLICT_DATA__", conflict_json)
        .replace("__VAX_DATA__", vax_json))

with open("ethiopia_market_activity_map.html", "w", encoding="utf-8") as f:
    f.write(html)

print("Wrote ethiopia_market_activity_map.html", len(html) / 1e6, "MB")
