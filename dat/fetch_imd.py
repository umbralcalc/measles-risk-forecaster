#!/usr/bin/env python3
"""Ingest deprivation per UTLA — a community-structure covariate NOT derived from
MMR coverage, to test whether it adds spatial cluster-prediction skill beyond
coverage (avenue #2). Measles under-vaccination often concentrates in specific
communities that area-average coverage smooths away; deprivation is one open
proxy for that structure.

Source: MHCLG English Indices of Deprivation 2019, File 11 (upper-tier local
authority summaries), IMD sheet. Average score (higher = more deprived) and the
proportion of neighbourhoods in the most-deprived 10% nationally.

2019 UTLA codes are crosswalked to the 2024 boundary set (the 2020-23 unitary
reorganisations); split counties inherit the old county's score.

Output: dat/imd_utla.csv  (code, name, imd_avg_score, imd_prop_deprived10)
Licence: Open Government Licence v3.0.

Usage: python3 dat/fetch_imd.py
"""

import csv
import os
import subprocess
import sys
import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "imd_utla.csv")
XLSX = os.path.join(HERE, "_imd_utla.xlsx")
URL = ("https://assets.publishing.service.gov.uk/media/5d8b3d7aed915d0369518030/"
       "File_11_-_IoD2019_Local_Authority_District_Summaries__upper-tier__.xlsx")

# 2019 upper-tier code -> 2024 boundary code(s) for reorganised areas.
OLD_TO_NEW = {
    "E10000002": ["E06000060"],               # Buckinghamshire (2020)
    "E10000021": ["E06000061", "E06000062"],  # Northamptonshire -> N/W Northants
    "E10000006": ["E06000063", "E06000064"],  # Cumbria -> Cumberland, Westmorland
    "E10000023": ["E06000065"],               # North Yorkshire (unitary 2023)
    "E10000027": ["E06000066"],               # Somerset (unitary 2023)
}


def main():
    p = subprocess.run(["curl", "-sfL", "--max-time", "90", URL, "-o", XLSX])
    if p.returncode != 0:
        raise RuntimeError("IMD download failed")

    wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
    ws = wb["IMD"]
    rows = list(ws.iter_rows(values_only=True))
    # columns: 0 code, 1 name, 4 avg score, 6 proportion in most deprived 10%
    out = {}
    for r in rows[1:]:
        code, name = r[0], r[1]
        if not code or not str(code).startswith("E"):
            continue
        rec = {"name": name, "imd_avg_score": r[4], "imd_prop_deprived10": r[6]}
        out[code] = rec
        for new in OLD_TO_NEW.get(code, []):
            out[new] = dict(rec)  # split successor inherits the old county's score

    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["code", "name", "imd_avg_score", "imd_prop_deprived10"])
        for code in sorted(out):
            rec = out[code]
            w.writerow([code, rec["name"], rec["imd_avg_score"], rec["imd_prop_deprived10"]])
    os.remove(XLSX)
    print(f"wrote {OUT}: {len(out)} UTLAs")


if __name__ == "__main__":
    sys.exit(main())
