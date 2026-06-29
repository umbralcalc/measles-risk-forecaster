#!/usr/bin/env python3
"""Ingest COVER MMR coverage by UTLA from the UKHSA dashboard API (the
susceptibility layer; gating check #2 confirmed it is exposed at UTLA).

For both doses (MMR1, MMR2) and both age milestones (stratum 24m, 5y), this
pulls the full annual coverage history for every England Upper Tier Local
Authority and writes a tidy long CSV. Each data row carries the ONS geography
code (e.g. E09000002), so downstream joins to the boundary/adjacency set (#4)
are by code, not by fragile name matching.

Output: dat/cover_mmr_coverage.csv
  columns: geography_code, geography, dose, stratum, year, date, coverage_pct

Licence: UK Health Security Agency, Open Government Licence v3.0. (SOURCES.md §3.)

Usage: python3 dat/fetch_coverage.py
"""

import csv
import json
import os
import subprocess
import sys
import time
import urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "cover_mmr_coverage.csv")
NODES = os.path.join(HERE, "utla_nodes.csv")

API = "https://api.ukhsa-dashboard.data.gov.uk"
GEO_TYPE = "Upper Tier Local Authority"
DOSES = ["MMR1", "MMR2"]


def get_json(url, retries=3):
    # Fetch via curl: the python.org Python on macOS ships without a usable CA
    # bundle, while curl uses the system trust store. Keeps JSON logic in Python.
    for attempt in range(retries):
        proc = subprocess.run(
            ["curl", "-sf", "--max-time", "60", url],
            capture_output=True, text=True,
        )
        if proc.returncode == 0:
            return json.loads(proc.stdout)
        if attempt == retries - 1:
            raise RuntimeError(
                f"curl failed ({proc.returncode}) for {url}: {proc.stderr}")
        time.sleep(1.5 * (attempt + 1))
    raise RuntimeError("unreachable")


def topic_base(dose):
    gt = urllib.parse.quote(GEO_TYPE)
    return (f"{API}/themes/immunisation/sub_themes/childhood-vaccines/"
            f"topics/{dose}/geography_types/{gt}")


def geographies(dose):
    return [g["name"] for g in get_json(topic_base(dose) + "/geographies")]


def coverage_rows(dose, geography):
    gq = urllib.parse.quote(geography)
    metric = f"{dose}_coverage_coverageByYear"
    url = (f"{topic_base(dose)}/geographies/{gq}/metrics/{metric}"
           f"?page_size=365")
    rows = []
    while url:
        page = get_json(url)
        for r in page["results"]:
            rows.append({
                "geography_code": r["geography_code"],
                "geography": r["geography"],
                "dose": dose,
                "stratum": r["stratum"],
                "year": r["year"],
                "date": r["date"],
                "coverage_pct": r["metric_value"],
            })
        url = page.get("next")
    return rows


def main():
    all_rows = []
    for dose in DOSES:
        geos = geographies(dose)
        print(f"{dose}: {len(geos)} UTLAs ...", flush=True)
        for i, g in enumerate(geos, 1):
            all_rows.extend(coverage_rows(dose, g))
            if i % 40 == 0:
                print(f"  {dose} {i}/{len(geos)}", flush=True)

    fields = ["geography_code", "geography", "dose", "stratum",
              "year", "date", "coverage_pct"]
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(all_rows)

    # --- summary ---
    codes = {r["geography_code"] for r in all_rows}
    years = sorted({r["year"] for r in all_rows})
    strata = sorted({r["stratum"] for r in all_rows})
    print(f"\nwrote {len(all_rows)} rows -> {os.path.basename(OUT)}")
    print(f"  UTLAs: {len(codes)}   years: {years[0]}-{years[-1]}   "
          f"strata: {strata}   doses: {DOSES}")

    # --- join check against the boundary/adjacency set (#4) ---
    if os.path.exists(NODES):
        with open(NODES) as f:
            boundary = {row["code"]: row["name"]
                        for row in csv.DictReader(f)}
        cover_names = {r["geography_code"]: r["geography"] for r in all_rows}
        only_cover = set(cover_names) - set(boundary)
        only_boundary = set(boundary) - set(cover_names)
        print(f"\njoin to boundaries (#4): "
              f"{len(set(cover_names) & set(boundary))} matched by code")
        if only_cover:
            print(f"  in COVER not boundaries ({len(only_cover)}): "
                  f"{[(c, cover_names[c]) for c in sorted(only_cover)]}")
        if only_boundary:
            print(f"  in boundaries not COVER ({len(only_boundary)}): "
                  f"{[(c, boundary[c]) for c in sorted(only_boundary)]}")
    else:
        print(f"\n(skip join check: {NODES} not found — run build_adjacency.py)")


if __name__ == "__main__":
    sys.exit(main())
