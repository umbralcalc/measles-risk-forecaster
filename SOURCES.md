# measles-risk-forecaster — SOURCES

Data sources, their granularity, licences, and access constraints. All sources
are **Open Government Licence v3.0 (OGL v3.0)**, open, England. Licences verified
2026-06-25 (gating check #1 + licence de-risking).

## Licence summary

| # | Source | Used for | Granularity | Licence | Verified |
|---|--------|----------|-------------|---------|----------|
| 1 | UKHSA dashboard API | measles cases | Nation, UKHSA Region × week | OGL v3.0 | own footer |
| 2 | GOV.UK measles epidemiology report | measles cases | UTLA (cumulative window) | OGL v3.0 (© Crown 2026) | own footer |
| 3 | COVER / vaccine uptake (UKHSA) | MMR coverage | UTLA × year | OGL v3.0 | own footer |
| 4 | ONS Open Geography Portal | UTLA boundaries + adjacency | UTLA | OGL v3.0 (+ OS) | ONS licences page |
| 5 | ONS population estimates | denominators | UTLA, mid-year | OGL v3.0 | ONS T&Cs |

All five carry the standard OGL *"except where otherwise stated"* carve-out for
third-party/personal data. Immaterial here: every source is already aggregated
and disclosure-controlled, so no personal data sits under the exception.

## 1. Measles cases — Region × week (the fast/coarse dynamics feed)

- **Endpoint:** `api.ukhsa-dashboard.data.gov.uk`, path
  `themes/infectious_disease/sub_themes/vaccine_preventable/topics/Measles`.
- **Geography types exposed:** **Nation** and **UKHSA Region only** — the 9
  regions (East Midlands, East of England, London, North East, North West,
  South East, South West, West Midlands, Yorkshire and Humber).
- **No UTLA/LTLA in the API.** Confirmed by walking `…/Measles/geography_types`
  (gating check #1). Clean JSON, fast cadence; fine in time, coarse in space.
- **Licence:** OGL v3.0 — *"All content is available under the Open Government
  Licence v3.0, except where otherwise stated"*, © Crown copyright.

## 2. Measles cases — UTLA (the fine/slow spatial layer)

- **Source:** GOV.UK "Confirmed cases of measles in England by month, age, region
  and upper tier local authority" report. Updated **every 2 weeks**.
- **Format:** **inline HTML tables — no CSV/ODS/Excel download.** Ingestion needs
  an HTML scraper; markup has varied year-to-year (2023 ≠ 2026 structure), so the
  parser must be version-aware and brittle by nature.
- **Granularity — CONFIRMED (2024, 2025, 2026 reports):** UTLA is **only ever a
  single cumulative total** over the report window. There is **no UTLA×month panel
  anywhere** — the page title "by month, age, region and UTLA" is misleading. Each
  report has just two tables: Table 1 = age×region (cumulative); Table 2 = UTLA×total
  (cumulative). Month counts appear only in prose or at region/national level, never
  crossed with UTLA. Consequences:
  - **The fine spatial layer carries spatial, not temporal, signal.** UTLA tells you
    *where within England* cases concentrate (a cross-section refreshed fortnightly),
    not *when*. Temporal dynamics live only at Region×week (source #1).
  - **Architecture follows from this:** Region×week drives the time course; UTLA
    cumulative shares spatially disaggregate it within region (areal-disaggregation
    framing — fits the CAR/GP spatial prior). This *is* the mismatched-granularity
    reconciliation the project centres on.
  - **A UTLA temporal panel only exists prospectively, by self-harvesting.** The live
    report is year-to-date cumulative; differencing successive fortnightly snapshots
    yields per-UTLA increments (for UTLAs above suppression in *both* snapshots — a
    UTLA crossing the <10 threshold shows up as absent→jump, which the censoring logic
    must handle). Proof-of-commit harvesting builds this panel going forward.
  - **Backtest depth at UTLA resolution is thin:** finalized history is one cumulative
    annual UTLA total per year (2024, 2025) plus 2026-to-date, unless reconstructed
    from Wayback snapshots of the then-live report. Honest UTLA scoring is
    cross-sectional risk-ordering against the next snapshot's increments.
  - **Sub-model C (reporting-lag nowcast) is a Region-level model only** — onset-week
    × delay data exists at Region, not UTLA.
- **Suppression is row-omission, not value-censoring.** Counts <10 do not appear
  as `<10` or a blank cell — **the UTLA is dropped from the table entirely.**
  Footnote: *"Case counts have been suppressed to not present any UTLAs with fewer
  than 10 cases."* → The censored set = (full ONS UTLA list) − (UTLAs listed).
  Ingestion must reconstruct it by joining against canonical UTLA geography (#4);
  the likelihood treats those UTLAs as interval-censored `[0,9]`, never zero.
- **Licence:** OGL v3.0, © Crown copyright 2026.

## 3. MMR coverage — UTLA × year (the susceptibility layer)

- **Source:** COVER (Cover of Vaccination Evaluated Rapidly) programme via the
  UKHSA dashboard API. **Gating check #2 resolved:** unlike measles cases, coverage
  IS exposed at **Upper Tier Local Authority via the API** — clean JSON, no scraping.
  Path: `immunisation/childhood-vaccines/{MMR1,MMR2}/geography_types/Upper Tier Local
  Authority`. Geography types: Nation, Region, United Kingdom, UTLA.
- **Coverage / structure (confirmed from live API 2026-06-25):**
  - **152 UTLAs** (full England set; matches ONS UTLA count).
  - **Annual cadence**, financial-year end (`date` = 31 March), **~12 years history**
    per UTLA (2014–2025; count=24 = 12 years × 2 strata).
  - Age milestone is the **`stratum`** field — **both `24m` and `5y` available for
    MMR1 *and* MMR2** (richer than "MMR1@24m, MMR2@5y"; supports real cohort
    accounting in sub-model A).
  - Metrics: `MMR{1,2}_coverage_coverageByYear` (coverage %, e.g. 88.09) and
    `MMR{1,2}_coverage_oneYearChange`.
- **Documented biases (to encode as explicit flags, not ignore):** London coverage
  may currently be *underestimated* (system/methodology change); small-population
  LAs are combined for disclosure; coverage is by current age-milestone, not
  historical cohort.
- **Licence:** OGL v3.0, © Crown copyright.

## 4. UTLA boundaries + adjacency (spatial prior)

- **Source:** ONS Open Geography Portal (`geoportal.statistics.gov.uk`),
  ArcGIS FeatureServer org `ESMARspQHYMw9BZ9`. **Gating check #4 resolved.**
- **Layer chosen:** `Counties_and_Unitary_Authorities_December_2024_Boundaries_UK_BGC`,
  filtered to England (`CTYUA24CD LIKE 'E%'`) → **153 areas**. Fields used:
  `CTYUA24CD`, `CTYUA24NM`. Fetched via `dat/fetch_geography.sh` (BGC = generalised
  clipped; topologically consistent, so neighbours share exact boundary vertices).
  - The older `Upper_Tier_Local_Authorities_December_2022` layer was **rejected**:
    it groups metropolitan districts into metropolitan counties (122 England areas)
    — wrong granularity — and predates the 2023 unitary reorganisations.
- **Reconciliation to COVER (153 vs 152):** the boundary set has 153 England UTLAs;
  COVER reports 152. The difference is the documented small-LA combination (e.g.
  City of London, Isles of Scilly). Joined by name at ingest; mismatches flagged,
  not dropped. (Same small-LA-combination caveat noted for source #3.)
- **Adjacency:** built by `dat/build_adjacency.py` (pure Python, no geometry lib —
  exact shared-segment detection on the topological boundaries). Result:
  **153 nodes, 360 edges, mean degree 4.77**; Isle of Wight and Isles of Scilly are
  isolated (no land border — handled by the CAR prior, which requires each
  connected component to carry at least one anchoring observation). Committed
  derived artifacts: `dat/utla_nodes.csv`, `dat/utla_adjacency.csv` (the Go spatial
  tests depend on them); the 7.4 MB source GeoJSON is gitignored (regenerable).
- **Licence:** OGL v3.0, with a **two-part** mandatory attribution:
  - *"Source: Office for National Statistics licensed under the Open Government
    Licence v3.0"*
  - *"Contains OS data © Crown copyright and database right [year]"*
  - The second (Ordnance Survey) line is easy to miss — it is OS data carried
    inside ONS boundaries and **must** appear wherever boundaries are reproduced.

## 5. Population denominators (per-capita susceptibility / importation pressure)

- **Source:** ONS mid-year population estimates.
- **Licence:** OGL v3.0 — ONS T&Cs: *"Most content on this website is subject to
  Crown copyright protection and is published under the Open Government Licence."*

## Attribution string (for the dashboard footer)

> Contains public sector information from the UK Health Security Agency (measles
> epidemiology and COVER vaccination coverage) and the Office for National
> Statistics (geographies and population), licensed under the Open Government
> Licence v3.0. Contains OS data © Crown copyright and database right [year].
> This project is a non-commercial public-interest methodological exercise, is not
> affiliated with or endorsed by UKHSA or ONS, and is not an official source of
> measles risk assessment — for that, consult UKHSA directly.

(Adds the OS Crown-copyright-and-database-right line missing from the PLAN.md draft.)
