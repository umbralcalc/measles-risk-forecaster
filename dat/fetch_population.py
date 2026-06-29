#!/usr/bin/env python3
"""Ingest ONS mid-year population estimates by UTLA (the denominator layer).

Population is the offset in the susceptibility->cases observation model (seam 2):
bigger areas see more cases for the same per-capita risk, so case counts are
modelled as log mu_i = log(pop_i) + alpha + beta * z_i. Source: Nomis dataset
NM_2002_1 (ONS population estimates), queried by the GSS codes from the boundary
set so it joins to the adjacency graph (#4) by code.

Output: dat/population_utla.csv  (code, name, year, population)

Licence: Office for National Statistics via Nomis, Open Government Licence v3.0.

Usage: python3 dat/fetch_population.py
"""

import csv
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
NODES = os.path.join(HERE, "utla_nodes.csv")
OUT = os.path.join(HERE, "population_utla.csv")

NOMIS = "https://www.nomisweb.co.uk/api/v01/dataset/NM_2002_1.data.csv"


def curl(url):
    p = subprocess.run(["curl", "-sf", "--max-time", "90", url],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"curl failed: {p.stderr}")
    return p.stdout


def main():
    with open(NODES) as f:
        codes = [row["code"] for row in csv.DictReader(f)]

    rows = []
    # Batch the GSS codes to keep the URL well within limits.
    BATCH = 50
    for i in range(0, len(codes), BATCH):
        chunk = codes[i:i + BATCH]
        url = (f"{NOMIS}?geography={','.join(chunk)}"
               f"&date=latest&gender=0&c_age=200&measures=20100"
               f"&select=date_name,geography_name,geography_code,obs_value")
        reader = csv.DictReader(curl(url).splitlines())
        for r in reader:
            rows.append({
                "code": r["GEOGRAPHY_CODE"],
                "name": r["GEOGRAPHY_NAME"],
                "year": r["DATE_NAME"],
                "population": r["OBS_VALUE"],
            })

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["code", "name", "year", "population"])
        w.writeheader()
        w.writerows(rows)

    got = {r["code"] for r in rows}
    missing = [c for c in codes if c not in got]
    years = sorted({r["year"] for r in rows})
    print(f"wrote {len(rows)} rows -> {os.path.basename(OUT)}  (years {years})")
    print(f"  matched {len(got)}/{len(codes)} boundary UTLAs")
    if missing:
        print(f"  MISSING population for: {missing}")


if __name__ == "__main__":
    sys.exit(main())
