#!/usr/bin/env python3
"""Ingest annual WHO European Region measles incidence — the importation-pressure
covariate (coarse importation index; see PLAN.md importation-feasibility verdict).

England's measles importations track upstream European activity (the 2024 surge
followed Europe's). This is the *coarse* covariate the feasibility study endorsed:
annual resolution is what the data reliably supports — monthly European series are
not cleanly machine-readable for the current period, and the UK import-count
calibration is prose-only and tiny. So the index is deliberately low-resolution
and wide-uncertainty, not a precise weekly forecast.

Source: Our World in Data "Reported cases of measles" (WHO-sourced), entity
"Europe (WHO)" = WHO European Region. CC-BY; underlying data WHO/UNICEF.

Output: dat/european_measles_annual.csv  (year, cases)

Usage: python3 dat/fetch_european_measles.py
"""

import csv
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "european_measles_annual.csv")
URL = ("https://ourworldindata.org/grapher/reported-cases-of-measles.csv"
       "?v=1&csvType=full")
ENTITY = "Europe (WHO)"
VALUE_COL = "Measles - number of reported cases"


def main():
    p = subprocess.run(["curl", "-sfL", "--max-time", "90", URL],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"curl failed: {p.stderr}")

    rows = []
    for r in csv.DictReader(p.stdout.splitlines()):
        if r["Entity"] == ENTITY and r[VALUE_COL].strip():
            rows.append({"year": int(r["Year"]), "cases": int(r[VALUE_COL])})
    rows.sort(key=lambda x: x["year"])

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["year", "cases"])
        w.writeheader()
        w.writerows(rows)

    recent = {r["year"]: r["cases"] for r in rows if r["year"] >= 2018}
    print(f"wrote {len(rows)} years -> {os.path.basename(OUT)}")
    print(f"  recent: {recent}")


if __name__ == "__main__":
    sys.exit(main())
