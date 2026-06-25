# measles-risk-forecaster — PLAN

A sub-national **measles spatial-risk map** for England: a calibrated, regularly-
refreshed estimate of which areas are primed for measles transmission given local
susceptibility, published as a frozen interactive dashboard with a proof-of-commit
calibration loop. One Go + stochadex repo, the public-health flagship companion to
the (completed) AMR work and the bathing-water companion, built to the same
honesty discipline as the other forecasters.

This file is the design and the road to here. Methodology lives in `README.md`;
data sources and their licences in `SOURCES.md`.

## Why this project, and the competitive position (stated up front)

Measles is a deliberately chosen *thin-market* public-health target. Unlike
respiratory viruses — where formal forecasting hubs (US FluSight, European
Respicast, the COVID-19 Forecast Hub) run standing, publicly-scored, multi-team
ensembles — there is no equivalent standing public measles *forecasting* operation
with a live calibration record. The measles literature is rich but mostly
retrospective outbreak-modelling papers, not a continuously-scored public product.
That is the gap we occupy.

The headline deliverable is a **spatial risk map**, not an outbreak-timing
forecast — and that choice is principled, not evasive. The definitive England &
Wales analysis established that as vaccination rose, the signature of spatially
*predictable* spread diminished and infection increasingly arrived from
unidentifiable, effectively random importation sources (Eilenberg/Grenfell-lineage
competing-risks work, 2020). So in a high-coverage setting:

- **Where** susceptibility pockets sit is partly predictable (coverage gaps
  accumulate slowly and spatially cluster) — this we can map and forecast.
- **When and where** an importation spark lands is largely irreducible.

We therefore forecast *spatial risk* (a susceptibility-and-exposure surface) and
are explicit that outbreak *timing* skill is bounded by importation randomness.
That decomposition — predictable pocket structure vs. irreducible spark — is the
distinctive, literature-backed honest finding, not a limitation we hide. It is the
measles analogue of the bathing-water skill-ceiling figure.

Why it clears the bar: open and live (UKHSA + COVER, OGL v3.0), unimpeachably
civic (nobody bets on measles), no incumbent calibrated public forecaster, and a
genuinely hard multi-resolution problem that is hard to replicate.

## What we forecast

For each English upper-tier local authority (UTLA), refreshed on each data update:

- **Headline — a calibrated transmission-risk score**: an estimate of local
  effective susceptibility translated into `R_local` (the expected secondary
  cases per importation given local immunity), and from it `P(self-sustaining
  transmission | importation)` i.e. P(R_local > 1).
- **Outbreak-potential**: conditional on a seeded case, the predictive
  distribution of resulting cluster size (a branching-process readout), surfacing
  the "large outbreak" tail rather than a point estimate.
- **Supporting nowcast**: a de-biased estimate of *current* activity by area,
  correcting the reporting-lag right-truncation (a sub-model, not the headline —
  see architecture).

We forecast distributions, simulating the local transmission process many times
and reading risk quantities off the ensemble.

**Named precisely.** The target the risk score is validated against is the
subsequent occurrence and size of laboratory-confirmed measles clusters by UTLA,
as published in the UKHSA epidemiology tables. We forecast *risk of* transmission,
scored against *observed* clusters — and we are explicit that a low-risk area can
still see an imported case (that is the irreducible part), so scoring rewards
calibrated risk ordering, not deterministic outbreak calling.

## Data

Full detail and licences in `SOURCES.md`. All OGL v3.0, open, England.

Three layers at **mismatched granularities** — reconciling them is the core
modelling work, not an afterthought:

- **Dynamics layer — case counts, Region × week.** UKHSA dashboard API
  (machine-readable Swagger API at `api.ukhsa-dashboard.data.gov.uk`, plus bulk
  downloads), measles cases by onset week and age group by UKHSA Region from
  Oct 2023. Fast cadence, ~9 regions — coarse in space, fine in time.
- **Fine spatial layer — case counts, UTLA × month.** The GOV.UK measles
  epidemiology reports (updated ~fortnightly), confirmed cases by month, age,
  region and UTLA. **Known constraint:** counts below 10 per UTLA are *suppressed*
  for disclosure control — so this layer is itself **interval-censored at the low
  end**, exactly where most UTLAs sit. Treated as censored data (see honesty
  rules), not dropped or zero-filled. *Gating check: confirm whether any UTLA/LTLA
  case granularity is exposed via the API or only via the report tables, and
  whether report tables are machine-parseable CSV/ODS.*
- **Susceptibility layer — MMR coverage, UTLA × year.** COVER programme via the
  dashboard, annual childhood vaccination coverage (MMR1 at 24 months, MMR2 at 5
  years) by UTLA, region and country. **Known caveats, stated in SOURCES.md:**
  London coverage may currently be *underestimated* due to system/methodology
  changes; small-population LAs are combined for disclosure; coverage is by current
  age-milestone, not historical cohort. These bias the susceptibility surface and
  are modelled/flagged, not ignored.
- **Spatial structure.** UTLA boundaries and adjacency (ONS open geographies) for
  the spatial smoothing prior; population denominators (ONS mid-year estimates) for
  per-capita susceptibility and importation pressure.

## The model

A **hierarchical spatial susceptibility-and-transmission model**, built from
three linked sub-models in the conditional-probability style that the measles
decision-modelling literature validates (the Bayesian-hierarchical integration of
dynamics + partial observability + coverage uncertainty; cf. the retrospective
outbreak decision-making work).

**Sub-model A — Susceptibility surface.** For UTLA *i*, age cohort *a*, time *t*,
estimate the effective susceptible fraction `s_iat` by accumulating cohort MMR
coverage (one- and two-dose, with assumed efficacy) net of natural immunity and
adjusting for the documented coverage biases. Because coverage is UTLA-annual and
clusters spatially, a **CAR / Gaussian-process spatial prior** over UTLAs pools
strength across neighbours — the disaggregation-by-spatial-regression approach
(INLA-SPDE-style) that the coverage-mapping literature uses to turn coarse areal
coverage into a smooth susceptibility field. This is the lever that lets sparse,
censored local data still produce a usable surface.

**Sub-model B — Local transmission.** Convert the susceptibility surface into a
local reproduction number `R_local,i = R0 · s_i` (R0 in the measles 12–18 range,
treated as an uncertain parameter, surface-pooled), and run a **stochastic
branching process** per UTLA seeded by importation pressure to get the
cluster-size predictive distribution and P(R_local > 1). The stochadex latent
factor carries a shared national epidemic anomaly (importation waves coupling
regions).

**Sub-model C — Reporting-lag nowcast.** A right-truncation correction on the
reporting triangle (onset week × reporting delay), so current activity is de-biased
before it informs the risk surface. This is well-trodden ground — the NobBS /
Bayesian-smoothing / generative-renewal nowcasting lineage, one of which used
measles as a test case and noted a slightly time-varying delay — so we adopt an
established method rather than invent, and cite it.

Sub-models linked via conditional probability so uncertainty propagates end to
end. Fit by **simulation-based inference**. Engine in `internal/risk`.

## The modelling unit (settled empirically)

A sweep (`cmd/unit-sweep`) over spatial unit (UTLA vs LTLA where available vs
region) and the spatial-prior strength / age-cohort resolution. Hypothesis to
test, not assume: UTLA with strong spatial pooling wins (fine enough to locate
pockets, coarse enough that the censored counts still inform), while LTLA is too
sparse to identify without heavy borrowing. Pick on out-of-sample cluster-
occurrence calibration, not intuition.

## Validation & honesty rules

- **Low-count suppression is first-class censoring.** UTLA counts <10 are
  interval-censored `[0,9]`, not zero. The likelihood treats them as censored,
  never as zero or as the cap. `cmd/suppression-ablation` quantifies the bias from
  naively zero-filling — an honest figure and a methodological contribution in its
  own right (mirrors the bathing-water and AMR censoring discipline).
- **The predictable-vs-irreducible decomposition — the load-bearing figure.**
  `cmd/skill-ceiling` partitions outcome variance into the spatially predictable
  component (susceptibility pocket structure) and the irreducible component
  (importation timing/location), explicitly reproducing and quantifying the
  high-coverage-unpredictability result the England & Wales literature established.
  We state this ceiling up front: we forecast *where is primed*, not *when it
  sparks*.
- **Backtest.** Expanding-window, no-leakage (`cmd/risk-backtest`), scoring the
  risk ordering against subsequent observed clusters, vs baselines: (a) raw
  most-recent case counts ("where it's already happening"), and (b) raw MMR
  coverage alone ("low coverage = high risk", the intuition the public already
  has). Beating *coverage-alone* — by adage of spatial pooling, age structure, and
  importation pressure — is the bar. If we can't, we say so.
- **Coverage-bias sensitivity (`cmd/coverage-sensitivity`).** Because London
  coverage may be underestimated and small LAs are combined, re-run the surface
  under stated alternative coverage assumptions and publish how much the risk map
  moves. Honesty about input bias propagated to output.
- **Proof of commit.** `cmd/forecast` commits each refresh's UTLA risk map *before*
  the next epidemiology-table release, to `data/predictions/`. `cmd/resolve`
  settles the risk ordering against the subsequently published clusters, refusing
  un-committed periods and asserting no leakage. Cadence follows the ~fortnightly
  report; the calibration curve gains points over a season.

## Scoring

- **Calibration of P(R_local > 1)** and of cluster-occurrence, as reliability
  curves — the running plot that gains points each release.
- **Ranked / ordinal skill** (does the risk ordering put subsequently-affected
  UTLAs high?) via AUC-style and rank-correlation metrics, since absolute
  occurrence is rare and ordering is the decision-relevant output.
- **CRPS** on the cluster-size predictive distribution where clusters occur.
- **Rare-event aware**: Brier skill score vs the coverage-alone baseline; we report
  skill *relative to the obvious heuristic*, never raw accuracy (which a
  "no outbreak anywhere" null scores well on).

## Repo structure

```
cmd/
  ingest-cases/        UKHSA dashboard API (region×week) + epidemiology tables (UTLA×month)
  ingest-coverage/     COVER MMR coverage (UTLA×year), with documented-bias flags
  ingest-geo/          ONS UTLA boundaries/adjacency + population denominators
  unit-sweep/          empirical spatial-unit / prior-strength selection
  forecast/            commit the UTLA risk map before the next report release
  resolve/             settle risk ordering vs published clusters, assert no leakage
  risk-backtest/       expanding-window scoring vs case-count + coverage-alone baselines
  suppression-ablation/ cost of treating <10 suppression as zero (censoring honesty)
  skill-ceiling/       predictable pocket structure vs irreducible importation
  coverage-sensitivity/ how far the map moves under alternative coverage assumptions
internal/
  risk/                susceptibility surface + branching-process transmission + ensemble
  nowcast/             reporting-triangle right-truncation correction (sub-model C)
  spatial/             CAR/GP spatial prior, UTLA adjacency, disaggregation
  data/                UKHSA/COVER/ONS loaders, suppression-censoring handling
data/
  raw/                 cached pulls (gitignored beyond samples)
  predictions/         committed risk maps (the proof-of-commit ledger)
SOURCES.md
README.md
PLAN.md
```

## Dashboard & distribution

- Static-baked via the dexetera pattern: UKHSA/COVER/ONS pulled and posteriors
  precomputed server-side, the branching-process simulator compiled to WASM so the
  reader runs per-UTLA outbreak ensembles client-side (inline driver, R2 binary).
- Hero widgets: (1) a **choropleth risk map** of England by UTLA, coloured by
  P(R_local > 1), click a UTLA for its susceptibility breakdown and cluster-size
  distribution; (2) a **susceptibility explorer** — adjust assumed coverage or
  efficacy and watch the risk surface respond, with the predictable/irreducible
  split made visible; (3) a **nowcast panel** showing recent weeks de-biased for
  reporting lag, with the truncation uncertainty shaded honestly.
- Landing-page hook: a recurring "measles risk update" with the map thumbnail as an
  SVG drop; each fortnightly resolve feeds the public calibration curve.
- Audience: a topical, high-public-interest, unimpeachably civic subject (the live
  MMR-coverage and resurgence debate) that anchors the public-health spine of the
  blog alongside AMR — non-bettable, no incumbent forecaster, no data owner to
  antagonise.

## Status

Draft / not started. Gating checks before build:
1. Confirm the finest **machine-readable** case granularity: does the UKHSA API
   expose sub-region (UTLA/LTLA) measles counts, or is UTLA only in the GOV.UK
   report tables — and are those tables parseable (CSV/ODS) at acceptable cadence?
   The whole spatial layer depends on this; resolve it first.
2. Pull COVER MMR coverage by UTLA and build the age-cohort susceptibility
   accounting, with the London-underestimate and small-LA-combination caveats
   encoded as explicit flags.
3. Implement the censored likelihood for <10 suppression and unit-test it
   (recover known rates from synthetic suppressed data) before any real inference.
4. Build the spatial prior on real UTLA adjacency and sanity-check that pooling
   produces a sensible susceptibility surface on a known period (e.g. the 2024
   London/West Midlands concentration) before adding transmission and nowcast.
5. One end-to-end committed risk map for one report cycle before scaling.

## Attribution (draft)

Contains public sector information from the UK Health Security Agency (measles
epidemiology and COVER vaccination coverage) and the Office for National
Statistics (geographies and population), licensed under the Open Government
Licence v3.0. This project is a non-commercial public-interest methodological
exercise, is not affiliated with or endorsed by UKHSA or ONS, and is not an
official source of measles risk assessment — for that, consult UKHSA directly.