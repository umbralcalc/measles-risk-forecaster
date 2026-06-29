#!/usr/bin/env python3
"""Ingest measles case data (the dynamics + fine-spatial layers; gating check #1).

Two sources at mismatched granularity (see SOURCES.md §1-2):

  * Region x week  -- UKHSA dashboard API, metric measles_cases_casesByOnsetWeek,
    9 UKHSA regions, weekly, with the in_reporting_delay_period flag (feeds the
    sub-model C nowcast). -> dat/measles_cases_region_week.csv

  * UTLA cumulative -- the GOV.UK epidemiology report HTML table ("Number of
    laboratory-confirmed measles cases by region and UTLA"). Cumulative per
    report window; suppression is ROW OMISSION (UTLAs with <10 cases are absent,
    not shown as 0 or "<10"). -> dat/measles_cases_utla.csv. The censored set is
    reconstructed downstream as (boundary UTLAs) - (listed UTLAs).

Both OGL v3.0, UK Health Security Agency.

Usage: python3 dat/fetch_cases.py
"""

import csv
import json
import os
import subprocess
import sys
import time
from html.parser import HTMLParser

HERE = os.path.dirname(os.path.abspath(__file__))
REGION_OUT = os.path.join(HERE, "measles_cases_region_week.csv")
UTLA_OUT = os.path.join(HERE, "measles_cases_utla.csv")

API = "https://api.ukhsa-dashboard.data.gov.uk"
MEASLES = (f"{API}/themes/infectious_disease/sub_themes/vaccine_preventable/"
           f"topics/Measles")
REGION_GT = "UKHSA%20Region"

REPORT = ("https://www.gov.uk/government/publications/"
          "measles-epidemiology-2023-to-2026/confirmed-cases-of-measles-in-"
          "england-by-month-age-region-and-upper-tier-local-authority-{year}")
REPORT_YEARS = [2024, 2025, 2026]


def curl(url, retries=3):
    for attempt in range(retries):
        p = subprocess.run(["curl", "-sfL", "--max-time", "60", url],
                           capture_output=True, text=True)
        if p.returncode == 0:
            return p.stdout
        if attempt == retries - 1:
            raise RuntimeError(f"curl {url} failed: {p.stderr}")
        time.sleep(1.5 * (attempt + 1))


def get_json(url):
    return json.loads(curl(url))


# --- Region x week from the API -------------------------------------------

def fetch_region_week():
    geos = [g["name"] for g in
            get_json(f"{MEASLES}/geography_types/{REGION_GT}/geographies")]
    rows = []
    for name in geos:
        from urllib.parse import quote
        url = (f"{MEASLES}/geography_types/{REGION_GT}/geographies/"
               f"{quote(name)}/metrics/measles_cases_casesByOnsetWeek"
               f"?page_size=500")
        while url:
            page = get_json(url)
            for r in page["results"]:
                rows.append({
                    "geography_code": r["geography_code"],
                    "region": r["geography"],
                    "year": r["year"],
                    "epiweek": r["epiweek"],
                    "date": r["date"],
                    "cases": r["metric_value"],
                    "in_reporting_delay_period": r["in_reporting_delay_period"],
                })
            url = page.get("next")
    with open(REGION_OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "geography_code", "region", "year", "epiweek", "date",
            "cases", "in_reporting_delay_period"])
        w.writeheader()
        w.writerows(rows)
    weeks = sorted({(r["year"], r["epiweek"]) for r in rows})
    print(f"region x week: {len(rows)} rows, {len(geos)} regions, "
          f"{len(weeks)} weeks {weeks[0]}..{weeks[-1]}")
    return rows


# --- UTLA cumulative from the report HTML ---------------------------------

class TableParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.tables, self.cur, self.row, self.cell, self.intd = [], None, None, "", False

    def handle_starttag(self, t, a):
        if t == "table":
            self.cur = []
        elif t == "tr" and self.cur is not None:
            self.row = []
        elif t in ("td", "th") and self.row is not None:
            self.cell, self.intd = "", True

    def handle_endtag(self, t):
        if t == "table" and self.cur is not None:
            self.tables.append(self.cur); self.cur = None
        elif t == "tr" and self.row is not None:
            self.cur.append(self.row); self.row = None
        elif t in ("td", "th") and self.intd:
            self.row.append(self.cell.strip()); self.intd = False

    def handle_data(self, d):
        if self.intd:
            self.cell += d


def fetch_utla():
    rows = []
    for year in REPORT_YEARS:
        try:
            html = curl(REPORT.format(year=year))
        except RuntimeError as e:
            print(f"  {year}: report not available ({e})")
            continue
        p = TableParser(); p.feed(html)
        utla_table = None
        for tb in p.tables:
            if tb and tb[0] and tb[0][0].lower().startswith("upper tier"):
                utla_table = tb
                break
        if not utla_table:
            print(f"  {year}: no UTLA table found")
            continue
        has_region = len(utla_table[0]) == 3
        n = 0
        for r in utla_table[1:]:
            if not r or not r[0]:
                continue
            utla = " ".join(r[0].split())  # normalise whitespace
            region = " ".join(r[1].split()) if has_region else ""
            total = r[-1].replace(",", "")
            try:
                cases = int(total)
            except ValueError:
                continue
            rows.append({"report_year": year, "utla": utla,
                         "region": region, "total_cases": cases})
            n += 1
        print(f"  {year}: {n} UTLAs listed (>=10 cases); the rest are "
              f"row-omitted / censored [0,9]")
    with open(UTLA_OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["report_year", "utla", "region",
                                          "total_cases"])
        w.writeheader()
        w.writerows(rows)
    print(f"UTLA cumulative: {len(rows)} listed rows across "
          f"{len({r['report_year'] for r in rows})} report years")
    return rows


def main():
    print("Fetching measles cases...")
    fetch_region_week()
    fetch_utla()
    print("Done.")


if __name__ == "__main__":
    sys.exit(main())
