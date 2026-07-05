#!/usr/bin/env python3
"""Build the static dashboard data bundle: compact UTLA geometry + the latest
committed risk map + national summaries, in one small JSON the widget inlines.

The dexetera pattern is static-baked: posteriors are precomputed server-side
(cmd/forecast, cmd/skill-ceiling) and the reader gets a self-contained page. This
script assembles that page's data.

  geometry : ONS BUC (ultra-generalised) England UTLAs, coords rounded to 3dp
  risk     : latest data/predictions/<date>.json (per-UTLA + national blocks)
  skill    : data/skill_ceiling.json

Output: app/measlesdash/data/dashboard.json

Usage: python3 app/measlesdash/build_data.py
"""

import glob
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(HERE, "data", "dashboard.json")

BASE = "https://services1.arcgis.com/ESMARspQHYMw9BZ9/arcgis/rest/services"
LAYER = "Counties_and_Unitary_Authorities_December_2024_Boundaries_UK_BUC"


def fetch_geometry():
    url = (f"{BASE}/{LAYER}/FeatureServer/0/query"
           "?where=CTYUA24CD+LIKE+'E%25'&outFields=CTYUA24CD,CTYUA24NM"
           "&returnGeometry=true&outSR=4326&f=geojson")
    p = subprocess.run(["curl", "-sf", "--max-time", "120", url],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"geometry fetch failed: {p.stderr}")
    return json.loads(p.stdout)


def round_ring(ring):
    return [[round(x, 3), round(y, 3)] for x, y in ring]


def round_geom(g):
    if g["type"] == "Polygon":
        return {"t": "P", "c": [round_ring(r) for r in g["coordinates"]]}
    if g["type"] == "MultiPolygon":
        return {"t": "M", "c": [[round_ring(r) for r in poly] for poly in g["coordinates"]]}
    return None


def main():
    gj = fetch_geometry()
    geom = {f["properties"]["CTYUA24CD"]: round_geom(f["geometry"])
            for f in gj["features"]}

    risk_path = sorted(glob.glob(os.path.join(ROOT, "data/predictions/*.json")))[-1]
    rm = json.load(open(risk_path))
    utlas = {e["code"]: {
        "n": e["name"], "s": e["susceptibility"], "pR": e["prob_r_local_gt_1"],
        "pL": e["prob_large_outbreak"], "imp": e["importation_adjusted_pressure"],
        "c50": e["cluster_p50"],
    } for e in rm["utlas"]}

    bundle = {
        "generated": rm["generated_date"],
        "geom": geom, "utlas": utlas,
        "importation": rm.get("national_importation_pressure"),
        "calibration": rm.get("cluster_calibration"),
        "joint": rm.get("joint_cooccurrence"),
        "caseLink": rm.get("susceptibility_case_link"),
        "skill": None,
    }
    sk = os.path.join(ROOT, "data/skill_ceiling.json")
    if os.path.exists(sk):
        bundle["skill"] = json.load(open(sk))

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(bundle, open(OUT, "w"), separators=(",", ":"))
    print(f"wrote {OUT}: {os.path.getsize(OUT)} bytes "
          f"({len(geom)} geoms, {len(utlas)} UTLAs, risk {os.path.basename(risk_path)})")


if __name__ == "__main__":
    sys.exit(main())
