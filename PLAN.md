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
- **Supporting nowcast**: a de-biased estimate of *current* activity,
  correcting the reporting-lag right-truncation (a sub-model, not the headline —
  see architecture). **Region-level only:** onset-week × reporting-delay data
  exists at Region, not UTLA (gating check #1), so the nowcast runs on the 9
  regions and informs the surface, rather than nowcasting each UTLA directly.

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
- **Fine spatial layer — case counts, UTLA cumulative (refreshed ~fortnightly).**
  The GOV.UK measles epidemiology reports. **Gating check #1 resolved (see
  SOURCES.md §2):** UTLA is exposed *only* in these reports (not the API), *only* as
  **inline HTML tables** (no CSV/ODS), and *only* as a **single cumulative total**
  per window — there is **no UTLA×month panel**. So this layer carries spatial
  signal (*where* within England cases concentrate), not temporal; the time course
  lives at Region×week. A UTLA temporal panel only accrues prospectively by
  differencing successive fortnightly snapshots (proof-of-commit harvesting).
  **Suppression is row-omission:** UTLAs with <10 cases are dropped from the table
  entirely (not shown as `<10`), so the censored set = (full ONS UTLA list) −
  (listed UTLAs), treated as interval-censored `[0,9]` (see honesty rules), never
  zero. This is exactly where most UTLAs sit.
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
  importation pressure — is the bar. If we can't, we say so. **Depth caveat
  (gating check #1):** at UTLA resolution, finalized history is one cumulative total
  per year (2024, 2025) plus the live 2026-to-date, so the at-launch backtest is
  thin; a proper UTLA temporal panel accrues prospectively via proof-of-commit
  harvesting (and, where available, reconstructed from archived report snapshots).
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

In build. Repo follows the stochadex house style (`pkg/`, `cfg/`, `dat/`, `app/`,
`nbs/`), matching the AMR sibling, rather than the `cmd/`+`internal/` sketch this
file originally proposed. Gating checks before build:
1. ✅ **Resolved (2026-06-25).** Finest machine-readable case granularity is
   **UKHSA Region** via the API; UTLA is report-only, **HTML tables (no CSV/ODS)**,
   **cumulative-only (no UTLA×month)**, with **row-omission** suppression. See
   SOURCES.md §1–§2; data-architecture consequences folded into the Data section.
2. ✅ **Resolved (2026-06-25).** COVER MMR coverage **is exposed at UTLA via the
   API** (152 LAs, annual, ~12yr history, 24m+5y strata for MMR1 and MMR2) — clean
   JSON, no scraping. See SOURCES.md §3. Age-cohort accounting + bias flags still
   to build on top.
3. ✅ **Done (2026-06-25).** Censored likelihood for <10 suppression implemented
   and unit-tested in `pkg/measles/censored_poisson.go`. Synthetic-recovery test:
   at 91% suppression, censored MLE recovers true rate 6.00→5.996 while naive
   zero-fill collapses to 0.93 (85% bias). Implements the stochadex
   `inference.LikelihoodDistribution` interface, SBI-ready.
4. ✅ **Done (2026-06-25).** Spatial prior built on **real** ONS UTLA adjacency
   (CUA Dec 2024, 153 England areas, 360 edges; `dat/fetch_geography.sh` +
   `dat/build_adjacency.py`). Intrinsic CAR smoother in
   `pkg/measles/spatial_car.go`, tested: 64% RMSE reduction recovering a smooth
   surface from noisy obs; fully-suppressed (no-obs) UTLAs imputed from neighbours
   (exact ICAR invariant) and beating zero-fill — the #3↔#4 link. Remaining for a
   full check: rerun on the *real* 2024 London/West Midlands data once cases +
   coverage are ingested and joined (overlaps with #5).
5. ✅ **Done (2026-06-29).** End-to-end committed risk map. `ingest-coverage`
   (`dat/fetch_coverage.py`) + `ingest-cases` (`dat/fetch_cases.py`: Region×week
   API + UTLA HTML scrape with row-omission censoring). Sub-models, all in
   `pkg/measles` and unit-tested (20 tests):
   - **A susceptibility** (`coverage.go`): COVER → effective susceptibility → CAR
     surface, COVER↔boundary mismatches reconciled. London ≈0.21 vs national ≈0.10.
   - **B transmission** (`transmission.go`): R_local=R0·s, negative-binomial
     branching process → P(R_local>1), cluster-size tail. Validated vs analytic
     mean total progeny 1/(1−m).
   - **C nowcast** (`nowcast.go`): reporting-lag right-truncation, Gamma-posterior
     interval; 95% coverage = 0.95 on synthetic truncated data.
   - **`cmd/forecast`** writes `data/predictions/<date>.json`. First committed map
     (2026-06-29): top risk = Hackney + London boroughs (the real 2025 hotspots);
     high-coverage areas (e.g. South Tyneside) below the herd-immunity threshold.

   **Next (post-gating, productionisation):** real-snapshot reporting-delay (replace
   the assumed curve), `resolve`/`risk-backtest`/ablation commands, dashboard.

## Stochadex integration (post-gating)

The gating-check kernels were built standalone (correctness-first) and are now
being wired into the stochadex engine across three seams:

1. ✅ **Done (2026-06-29).** Sub-model B as a stochadex `Iteration`:
   `BranchingProcessIteration` (`pkg/measles/transmission_iteration.go` +
   `transmission_settings.yaml`) advances the outbreak one generation per step;
   an ensemble of seeds is the cluster-size distribution. Runs under the
   stochadex partition coordinator and test harness; shares `nextGeneration()`
   with the standalone kernel so the two cannot drift. This is the form that
   compiles to WASM for the dashboard's client-side per-UTLA ensembles.

   **Single-partition is correct for the per-UTLA marginals.** The *joint*
   co-occurrence question needed multi-partition — built as a genuine 2-partition
   stochadex simulation (`pkg/measles/joint_simulation.go`): a shared
   `NationalImportationIteration` (partition 0) draws a national seed total from
   the importation-index band and feeds a `JointOutbreakIteration` (partition 1,
   all UTLAs in one state vector) via `params_from_upstream`. The shared latent is
   what correlates UTLAs — validated: mean outbreak correlation **0.37 vs 0.00**
   for the independent baseline, and **6.7× the national-total variance**. Committed
   in the map as `joint_cooccurrence`: co-occurrence *breadth* (≈39 UTLAs mean, 61
   p90) and *relative* tail. Only breadth and correlation structure are honest
   outputs.

   ✅ **Susceptible depletion (2026-07-05).** The transmission model is now a
   stochastic SIR, not a pure branching process: R_eff = R_local,0 · S/S_0 declines
   as the local susceptible pool burns down (`simulateCluster` in transmission.go;
   the single-UTLA and joint Iterations deplete identically). This replaces the
   arbitrary cap with a mechanistic final-size distribution. **Honest finding:**
   depletion alone does NOT make absolute sizes realistic — homogeneous mixing
   across a UTLA's full susceptible pool (s·N) gives a takeoff that infects ~all of
   it (Hackney → ~67k vs observed 133). So the pool is bounded by a
   `ReachableClusterSize` (a school/neighbourhood-scale well-mixed community).

   ✅ **Calibrated (2026-07-05).** `pkg/measles/calibration.go` fits
   `ReachableClusterSize` (and an importation scale) to the observed censored UTLA
   totals via SBI. The observed total is a compound of importations × depleting
   chains, so `mu_i = importScale · receptivity_i · meanChainSize(R_local,i, RC)`;
   `meanChainSize` (precomputed grid) saturates near RC for supercritical UTLAs but
   stays small for subcritical ones, so the high-vs-low-susceptibility *contrast*
   identifies RC. Validated by synthetic recovery (true 400 → CI [223,556]).
   **Real 2025 fit: RC ≈ 400–510 susceptibles, CI ~[210, 1120]** (wide — a
   structural parameter weakly pinned by one year's cross-section). The forecast
   now calibrates on the latest *complete* year and uses the fitted RC, so absolute
   outbreak sizes are **community-scaled and data-driven** (Hackney p50 ≈ 500, vs
   133 observed — right order) rather than 67k or an arbitrary guess. A taken-off
   outbreak fills the reachable community, so conditional sizes saturate near RC;
   **P(R_local>1)** and **P(large outbreak)** remain the robust discriminators.
   Committed as `cluster_calibration`. Future: multi-year / spatial-heterogeneity
   refinement to tighten the wide CI.
2. ✅ **Done (2026-06-29).** Fitting via stochadex SBI. `pkg/measles/caselink.go`
   fits the susceptibility→cases observation model
   `log μ_i = α + log(pop_i/popBar) + β·z_i` to the **real censored UTLA counts**
   using `inference.ComputePosterior` + `inference.UniformPrior` and the censored
   Poisson likelihood (#3) — iterated importance sampling (population Monte Carlo,
   ESS ≈ 12.8k/20k). Validated by synthetic recovery (α,β within ~1 sd). Real
   finding: **β = 0.80, 95% CI [0.73, 0.87]** — susceptibility significantly
   predicts where cases concentrated. Committed into the risk map as
   `susceptibility_case_link`. Needed `ingest-population` (`dat/fetch_population.py`,
   ONS via Nomis). *Note: programmatic SBI, not the `cfg/*.yaml` codegen path,
   which is heavier and built for time series — this cross-sectional fit is cleaner
   via `ComputePosterior` directly.*
3. 🔨 **Dashboard started (2026-07-05).** `app/measlesdash/` — a self-contained
   static-baked widget (dexetera pattern, house style like AMR's `app/`):
   `build_data.py` bakes the committed risk map + ONS BUC geometry into a compact
   `data/dashboard.json`; `build_widget.py` inlines it into `widget.html` → a
   dependency-free `index.html`. **Hero choropleth done:** England UTLAs coloured by
   a metric toggle (susceptibility / P(Rₗ>1) / P(large outbreak) / importation
   pressure), hand-rolled SVG (no CDN deps), hover + click-to-detail, sequential
   single-hue risk ramp with the herd-immunity threshold marked, both light/dark
   themes, plus a stat strip surfacing the skill ceiling, importation index,
   case-link β, calibrated cluster size, and joint co-occurrence.
   **Next:** the interactive WASM branching-ensemble (compile
   `BranchingProcessIteration` per seam 1 to WASM for client-side per-UTLA outbreak
   simulation on click) + the susceptibility explorer; optionally `cmd/forecast` →
   stochadex CLI/YAML.

## Skill ceiling — the load-bearing honest figure

✅ **Done (2026-07-03).** `cmd/skill-ceiling` (`pkg/measles/skill_ceiling.go`)
computes the predictable-vs-irreducible decomposition on **real 2025 data**,
writing `data/skill_ceiling.json`:

- **Spatial: McFadden R² = 0.51, irreducible = 0.49** (LR = 759, hugely
  significant). Susceptibility explains ~half of *where* cases concentrated beyond
  population; **the other half is irreducible spark placement** — exactly the
  high-coverage-unpredictability result the England & Wales literature established,
  now quantified. This is the measles analogue of the bathing-water skill ceiling.
- **Temporal: momentum R² = 0.81, lag-1 autocorr = 0.88** — but explicitly flagged
  (in code and output) as **transmission persistence, not importation skill**. It
  is an upper bound; conflating it with importation predictability is the trap.
- Validated: the decomposition reports McFadden R²≈0.85 when susceptibility truly
  drives synthetic cases, ≈0.00 when cases follow population only.

## Importation predictability — feasibility verdict (2026-07-03)

Can a shared national importation latent process add temporal/co-occurrence
forecast skill? **Possible in principle, mechanism confirmed, but data-limited.**
- **Identifiability exists but coarse:** UKHSA quarterly "Lab-confirmed measles,
  rubella and mumps" reports classify cases imported / import-related / indigenous
  (Q1 2025: 14.6% imported, 3.5% import-related) — national, quarterly, small-count.
- **Leading covariate available:** ECDC monthly measles cases by country (EU/EEA);
  WHO European Region 127,350 in 2024 (2× 2023). England's surge followed Europe's.
  Flight-passenger dynamics (Nature 2024) usable for receptivity weights w_i.
- **Verdict:** enough to establish the *mechanism*, not (yet) to fit a
  well-identified weekly λ_nat(t). Build as a **coarse, wide-uncertainty national
  importation-pressure index** feeding seed counts — not a precise weekly forecast.
  Improves as quarters accrue. Do NOT sell the R²=0.81 momentum as importation skill.

✅ **Built (2026-07-03).** `pkg/measles/importation.go` + `dat/fetch_european_measles.py`.
The index is WHO European Region measles incidence (OWID, annual) relative to a
baseline geomean, with a deliberately **wide band** from the [0.5, 1.5] elasticity
uncertainty on the Europe→UK link. 2024: **~4–5× baseline, band [2×, 8×]**; 2021
≈0 (pandemic collapse). Distributed to UTLAs by population-share receptivity and
combined with local amplification P(self-sustaining) into a per-UTLA
`importation_adjusted_pressure` score (the co-occurrence mechanism), committed in
the risk map. Honest resolution is **annual** — monthly European series are not
reliably machine-readable for 2023–26. Result is a distinct lens: pure P(large)
tops with London boroughs (per-importation risk); importation-adjusted pressure
tops with large counties (expected outbreak count). Future refinement:
flight-connectivity receptivity weights; UKHSA quarterly import-count calibration.

All five sources confirmed OGL v3.0 (SOURCES.md); ONS geography needs the extra
OS Crown-copyright-and-database-right attribution line.

## Attribution (draft)

Contains public sector information from the UK Health Security Agency (measles
epidemiology and COVER vaccination coverage) and the Office for National
Statistics (geographies and population), licensed under the Open Government
Licence v3.0. This project is a non-commercial public-interest methodological
exercise, is not affiliated with or endorsed by UKHSA or ONS, and is not an
official source of measles risk assessment — for that, consult UKHSA directly.