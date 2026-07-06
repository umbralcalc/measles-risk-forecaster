# measles-risk-forecaster

A calibrated, sub-national **measles transmission-risk surface** for England's 153
upper-tier local authorities (UTLAs), built on UKHSA, COVER and ONS open data and
the [stochadex](https://github.com/umbralcalc/stochadex) simulation engine. For
each area it estimates local effective susceptibility, converts it to the chance
an introduction leads to self-sustaining transmission, and simulates the resulting
outbreak-size distribution — then scores that risk ordering openly against the
baseline the public already has: *"low coverage = high risk."*

The point of the project is **calibration and an honest decomposition, not a claim
to beat the obvious heuristic** — and, as it turns out, the headline result is a
negative one, stated up front.

## The honest headline

1. **Coverage is a near-sufficient spatial statistic.** In a no-leakage,
   out-of-sample backtest of *which areas get clusters* (≥10 confirmed cases),
   the naive coverage heuristic scores **AUC 0.791**. Our full spatial surface
   reaches **0.799** — a +0.008 edge from spatial pooling that is within noise over
   the two complete years available. **No modelling we tried robustly beats plain
   coverage.** We say so.

2. **The predictable-vs-irreducible decomposition is the real finding.** A deviance
   decomposition (`cmd/skill-ceiling`) shows susceptibility explains **~51%** of
   *where* 2025 cases concentrated beyond population; the other **~49% is
   irreducible spark placement** — reproducing, with a live number, the
   high-coverage-unpredictability result the England & Wales literature
   established. We forecast *where is primed*, not *when it sparks*.

3. **The two most promising improvement avenues were tested and failed, honestly.**
   *Age-cohort accumulation* (using the school-age coverage window, not a single
   milestone) helps 2024 but hurts 2025 — no robust gain (0.794). *Deprivation*
   (MHCLG IMD 2019) is near chance (**AUC 0.545**): measles under-vaccination is
   community-specific (Orthodox Jewish, Traveller, vaccine-hesitant pockets), not
   deprivation-driven, so area-average deprivation is a poor proxy.

4. **What the machinery *is* good for.** Even where it doesn't beat the heuristic,
   the pipeline delivers things the heuristic cannot: a calibrated per-area risk
   with propagated uncertainty, a **data-driven** outbreak scale, and inference
   that is honest about the disclosure censoring in the source data (below).

## The model

Three linked sub-models, fit by simulation-based inference; engine in
`pkg/measles`.

**A — Susceptibility surface.** For UTLA *i*, effective susceptibility from
two-dose MMR accounting,

```
s_i = 1 − ( c2·e2 + (c1−c2)·e1 )        c1,c2 = MMR1,MMR2 coverage;  e1,e2 = efficacy
```

then **CAR-smoothed** over the real ONS adjacency graph so sparse, disclosure-
suppressed areas borrow strength from neighbours (the smoother reduces error 64%
on a synthetic surface and imputes unobserved areas by the exact neighbour-mean).

**B — Local transmission.** `R_local = R0·s` (R0 ~ 12–18), run as a **stochastic
SIR / susceptible-depleting branching process**: the effective reproduction number
declines as the local susceptible pool burns down, so outbreaks self-limit rather
than growing to an arbitrary cap. It reads off P(R_local > 1) and the cluster-size
distribution, validated against the analytic mean total progeny 1/(1−m). The same
kernel is a stochadex `Iteration`, and a two-partition version with a shared
national importation latent reproduces genuine outbreak **co-occurrence**
(correlation 0.37 vs 0.00 independent).

**C — Reporting-lag nowcast.** A right-truncation correction on the Region×week
series (Gamma-posterior interval), 95% coverage on synthetic truncated data.
Region-level, because onset-week × delay data exists only at Region.

**Calibration.** The observed UTLA total is a compound of importations × depleting
chains, so `μ_i = importScale · receptivity_i · meanChainSize(R_local,i,
reachableCluster)`. Because the mean chain size saturates near the reachable
cluster for supercritical areas but stays small for subcritical ones, the
high-vs-low-susceptibility *contrast* identifies the reachable cluster size —
fitted to observed outbreak sizes at **~400 susceptibles** (95% CI [280, 520]).

**Censoring is first-class.** UTLA counts below 10 are *suppressed* for disclosure
(the area is dropped from the table, not shown as 0), so the fine spatial layer is
interval-censored `[0,9]` exactly where most areas sit. The likelihood treats it
as censored, proven by synthetic recovery: at **91% suppression** it recovers a
true rate of 6.00 as **5.996**, while naive zero-filling collapses to 0.93 (an 85%
bias). Data sources and the verified API surface are in [SOURCES.md](SOURCES.md).

## Results

All ranking scores are out-of-sample (coverage as-of year *t−1*; the surface never
sees case data), target = "UTLA had ≥10 confirmed cases", pooled over 2024–25
(~90 clusters).

| finding | result |
|---|---|
| **Does anything beat coverage-alone?** | coverage **0.791**; spatial surface 0.799 (+0.008, noise); cohort accumulation 0.794; **no robust gain** (AUC) |
| **Deprivation as a covariate** | AUC **0.545** — near chance; wrong proxy for community-specific hesitancy |
| **Predictable vs irreducible** | susceptibility explains **51%** of where cases concentrate beyond population; **49% irreducible** (McFadden R², LR 759) |
| **Does susceptibility predict cases at all?** | yes — fitted link **β ≈ 0.8, 95% CI excludes 0** (censored-likelihood SBI) |
| **Censoring, done right** | recovers 6.00 → **5.996** at 91% suppression; zero-fill → 0.93 |
| **Outbreak scale, calibrated** | reachable cluster **~400** susceptibles, CI [280, 520] |

## Why we can't (yet) beat the heuristic — and what would

The skill ceiling already says ~half the spatial variance is irreducible
importation randomness, and coverage captures most of the predictable half. The
one place genuinely-new predictive information plausibly lives is the
**teen/young-adult susceptible cohorts** — those who missed MMR in the 1998–2004
Wakefield-scare dip and drove a real share of 2024–25 cases — but their coverage is
only knowable at UTLA level from pre-2014 data, which COVER does not reach. With
the available open data, coverage is near the attainable spatial ceiling. The
remaining honest lever is simply **more seasons** (two years / ~90 clusters cannot
resolve a small edge).

Because of this, the project ships **no interactive dashboard**: we will not
present a "better forecast" until a backtest can quote a pooled AUC that beats
coverage-alone.

## Reproduce

```bash
go test ./...                       # 39 tests: kernels, SBI recovery, spatial, calibration
go run ./cmd/forecast               # commit a per-UTLA risk map to data/predictions/
go run ./cmd/skill-ceiling          # the predictable-vs-irreducible decomposition
go run ./cmd/risk-backtest          # the out-of-sample AUC verdict (the honest headline)
```

Data ingestion (Python + curl, all OGL v3.0 / CC-BY):

```bash
dat/fetch_geography.sh              # ONS UTLA boundaries -> adjacency graph
python3 dat/build_adjacency.py
python3 dat/fetch_coverage.py       # COVER MMR coverage by UTLA
python3 dat/fetch_cases.py          # UKHSA Region×week + UTLA report scrape (censored)
python3 dat/fetch_population.py     # ONS mid-year population (Nomis)
python3 dat/fetch_european_measles.py  # WHO European incidence (importation covariate)
python3 dat/fetch_imd.py            # MHCLG deprivation (tested, no signal)
```

## Attribution

Contains public sector information from the UK Health Security Agency (measles
epidemiology and COVER vaccination coverage), the Office for National Statistics
(geographies and population), and the Ministry of Housing, Communities & Local
Government (indices of deprivation), licensed under the
[Open Government Licence v3.0](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/);
contains OS data © Crown copyright and database right 2024. European measles
incidence: WHO via [Our World in Data](https://ourworldindata.org/measles)
(CC-BY 4.0). A non-commercial public-interest methodological exercise, not
affiliated with or endorsed by any of these bodies, and **not an official source
of measles risk assessment** — for that, consult UKHSA directly.
