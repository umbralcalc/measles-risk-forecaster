// Command forecast produces one committed UTLA measles risk map (gating check
// #5, the proof-of-commit ledger). It wires the three sub-models together:
//
//	A (susceptibility)  COVER coverage -> effective susceptibility -> CAR surface
//	B (transmission)    R_local = R0*s -> branching process -> P(R_local>1),
//	                    cluster-size distribution, large-outbreak tail
//	C (nowcast)         Region x week, reporting-lag de-biased (national summary)
//
// and writes data/predictions/<date>.json — the map committed before the next
// epidemiology-table release, so it can later be scored against observed clusters.
//
// Usage:
//
//	go run ./cmd/forecast [-datadir dat] [-outdir data/predictions] [-seed 1]
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"math/rand/v2"
	"os"
	"path/filepath"
	"sort"
	"time"

	"github.com/umbralcalc/measles-risk-forecaster/pkg/measles"
	"github.com/umbralcalc/stochadex/pkg/inference"
)

type riskEntry struct {
	Code              string  `json:"code"`
	Name              string  `json:"name"`
	Susceptibility    float64 `json:"susceptibility"`
	Imputed           bool    `json:"susceptibility_imputed"`
	ProbRLocalGt1     float64 `json:"prob_r_local_gt_1"`
	MedianRLocal      float64 `json:"median_r_local"`
	ClusterP50        int     `json:"cluster_p50"`
	ClusterP90        int     `json:"cluster_p90"`
	ClusterP99        int     `json:"cluster_p99"`
	ProbLargeOutbreak float64 `json:"prob_large_outbreak"`
	// ImportationPressure = national index x receptivity x P(self-sustaining):
	// a coarse relative expected-outbreak-pressure score (co-occurrence-aware).
	ImportationPressure float64 `json:"importation_adjusted_pressure"`
}

type nowcastEntry struct {
	AgeWeeks     int     `json:"age_weeks"`
	Observed     float64 `json:"observed"`
	Completeness float64 `json:"completeness"`
	Point        float64 `json:"point"`
	Lower        float64 `json:"lower"`
	Upper        float64 `json:"upper"`
}

// caseLinkBlock is the seam-2 SBI result: the posterior on how strongly the
// susceptibility surface predicts where cases actually concentrated, fitted to
// the censored UTLA counts. Committed alongside the map so the calibration claim
// travels with the forecast.
type caseLinkBlock struct {
	ReportYear     int       `json:"report_year"`
	Alpha          float64   `json:"alpha"`
	Beta           float64   `json:"beta"`
	BetaStd        float64   `json:"beta_std"`
	BetaCI95       []float64 `json:"beta_ci95"`
	BetaPredictsCases bool   `json:"beta_predicts_cases"`
	LogMarginalLik float64   `json:"log_marginal_lik"`
	EffectiveSampleSize float64 `json:"effective_sample_size"`
	NumParticles   int       `json:"num_particles"`
}

// importationBlock is the coarse national importation-pressure index committed
// with the map: European measles activity relative to baseline, wide-band.
type importationBlock struct {
	Year            int     `json:"european_data_year"`
	EuropeanCases   float64 `json:"european_cases"`
	BaselineGeomean float64 `json:"baseline_geomean"`
	RelativeCentral float64 `json:"relative_central"`
	RelativeLow     float64 `json:"relative_low"`
	RelativeHigh    float64 `json:"relative_high"`
	Note            string  `json:"note"`
}

// jointBlock summarises the joint co-occurrence simulation: a shared national
// importation latent feeding all UTLAs, so the national-total tail and the number
// of UTLAs lighting up together are captured with their correct correlation.
type jointBlock struct {
	NScenarios            int       `json:"n_scenarios"`
	SeedBand              []float64 `json:"national_seed_band"`
	ExpectedUTLAsOutbreak float64   `json:"expected_utlas_with_outbreak"`
	P90UTLAsOutbreak      float64   `json:"p90_utlas_with_outbreak"`
	// Relative tail fatness (p90/p50, p99/p50 of national activity), the
	// co-occurrence signature. Reported as ratios, NOT absolute case counts,
	// because the pure branching process has no susceptible depletion and is
	// therefore uncalibrated in absolute magnitude.
	NationalTailRatioP90 float64 `json:"national_tail_ratio_p90"`
	NationalTailRatioP99 float64 `json:"national_tail_ratio_p99"`
	Note                 string  `json:"note"`
}

// calibrationBlock is the reachable-cluster-size calibration: the well-mixed
// susceptible community per introduction, fitted to the observed UTLA totals so
// absolute outbreak sizes are data-driven rather than assumed.
type calibrationBlock struct {
	ReachableCluster    float64   `json:"reachable_cluster"`
	ReachableClusterCI  []float64 `json:"reachable_cluster_ci95"`
	ImportScale         float64   `json:"import_scale"`
	R0Used              float64   `json:"r0_used"`
	EffectiveSampleSize float64   `json:"effective_sample_size"`
	ReportYear          int       `json:"report_year"`
}

type riskMap struct {
	GeneratedDate string         `json:"generated_date"`
	DataVintage   map[string]any `json:"data_vintage"`
	Model         map[string]any `json:"model"`
	ClusterCalibration *calibrationBlock `json:"cluster_calibration"`
	Attribution   string         `json:"attribution"`
	SusceptibilityCaseLink *caseLinkBlock `json:"susceptibility_case_link"`
	NationalImportationPressure *importationBlock `json:"national_importation_pressure"`
	JointCoOccurrence *jointBlock `json:"joint_cooccurrence"`
	UTLAs         []riskEntry    `json:"utlas"`
	NationalNowcast []nowcastEntry `json:"national_nowcast"`
}

const attribution = "Contains public sector information from UKHSA (measles " +
	"epidemiology and COVER vaccination coverage) and the Office for National " +
	"Statistics (geographies and population), licensed under the Open Government " +
	"Licence v3.0. Contains OS data © Crown copyright and database right 2024. " +
	"European measles incidence: WHO via Our World in Data (CC-BY 4.0). " +
	"Not an official source of measles risk assessment."

func main() {
	datadir := flag.String("datadir", "dat", "directory of ingested CSVs")
	outdir := flag.String("outdir", "data/predictions", "committed risk-map directory")
	seed := flag.Uint64("seed", 1, "RNG seed (committed maps are reproducible)")
	stratum := flag.String("stratum", "5y", "coverage age milestone (24m|5y)")
	flag.Parse()

	rng := rand.New(rand.NewPCG(*seed, *seed))

	// --- load layers ---
	graph, err := measles.LoadAdjacency(
		filepath.Join(*datadir, "utla_nodes.csv"),
		filepath.Join(*datadir, "utla_adjacency.csv"))
	must(err, "load adjacency")

	coverage, err := measles.LoadCoverage(filepath.Join(*datadir, "cover_mmr_coverage.csv"))
	must(err, "load coverage")

	regionWeek, err := measles.LoadRegionWeekCases(
		filepath.Join(*datadir, "measles_cases_region_week.csv"))
	must(err, "load region-week cases")

	// --- sub-model A: susceptibility surface ---
	surf, err := measles.BuildSusceptibilitySurface(
		graph, coverage, *stratum,
		measles.DefaultMMR1Efficacy, measles.DefaultMMR2Efficacy,
		3.0 /*tau*/, 50.0 /*obsPrecision*/)
	must(err, "build susceptibility surface")

	// --- population (denominator for susceptible-depletion pools) ---
	population := loadPopulationVector(graph, *datadir)

	// --- seam 2: fit the susceptibility->cases link (SBI) ---
	caseLink := fitCaseLink(graph, surf, *datadir, *seed)

	// --- coarse national importation-pressure index + per-UTLA receptivity ---
	imp, receptivity := importationPressure(graph, *datadir, time.Now().Year())

	// --- calibrate the reachable-cluster size to observed outbreak sizes ---
	calib := calibrateCluster(graph, surf, receptivity, *datadir)

	// --- sub-model B: per-UTLA transmission risk (with susceptible depletion) ---
	tp := measles.DefaultTransmissionParams()
	if calib != nil {
		tp.ReachableClusterSize = calib.ReachableCluster
	}
	entries := make([]riskEntry, graph.NumNodes())
	for i := 0; i < graph.NumNodes(); i++ {
		popI := 200000.0 // fallback if population unavailable
		if population != nil {
			popI = population[i]
		}
		r := measles.ForecastUTLARisk(surf.Smoothed[i], popI, tp, rng)
		e := riskEntry{
			Code: graph.Codes[i], Name: graph.Names[i],
			Susceptibility: round3(surf.Smoothed[i]), Imputed: surf.Imputed[i],
			ProbRLocalGt1: round3(r.ProbRLocalGt1), MedianRLocal: round2(r.MedianRLocal),
			ClusterP50: r.ClusterP50, ClusterP90: r.ClusterP90, ClusterP99: r.ClusterP99,
			ProbLargeOutbreak: round3(r.ProbLargeOutbreak),
		}
		if imp != nil && receptivity != nil {
			e.ImportationPressure = round4(
				imp.RelativeCentral * receptivity[i] * r.ProbRLocalGt1)
		}
		entries[i] = e
	}
	sort.Slice(entries, func(a, b int) bool {
		if entries[a].ProbLargeOutbreak != entries[b].ProbLargeOutbreak {
			return entries[a].ProbLargeOutbreak > entries[b].ProbLargeOutbreak
		}
		return entries[a].Susceptibility > entries[b].Susceptibility
	})

	// --- joint co-occurrence simulation (shared importation latent) ---
	joint := jointCoOccurrence(surf, receptivity, population, imp, tp.ReachableClusterSize, *seed)

	// --- sub-model C: national reporting-lag nowcast ---
	nowcast, latestWeek := nationalNowcast(regionWeek)

	// --- assemble & commit ---
	rm := riskMap{
		GeneratedDate: time.Now().Format("2006-01-02"),
		DataVintage: map[string]any{
			"coverage_stratum": *stratum,
			"latest_case_week": latestWeek,
		},
		Model: map[string]any{
			"r0_min": tp.R0Min, "r0_max": tp.R0Max, "dispersion": tp.Dispersion,
			"num_sims": tp.NumSims, "large_outbreak_threshold": tp.LargeOutbreakThreshold,
			"effective_pool_fraction": tp.EffectivePoolFraction,
			"reachable_cluster_size": tp.ReachableClusterSize,
			"car_tau": 3.0, "mmr1_efficacy": measles.DefaultMMR1Efficacy,
			"mmr2_efficacy": measles.DefaultMMR2Efficacy, "seed": *seed,
			"cluster_size_caveat": "Reachable-cluster size is calibrated to observed " +
				"outbreak sizes (see cluster_calibration) but with wide uncertainty; a " +
				"taken-off outbreak fills the reachable community, so absolute sizes " +
				"saturate near it. P(R_local>1) and P(large_outbreak) are the robust, " +
				"calibration-free discriminators.",
		},
		ClusterCalibration: calib,
		Attribution: attribution, SusceptibilityCaseLink: caseLink,
		NationalImportationPressure: imp, JointCoOccurrence: joint,
		UTLAs: entries, NationalNowcast: nowcast,
	}

	must(os.MkdirAll(*outdir, 0o755), "mkdir outdir")
	outpath := filepath.Join(*outdir, rm.GeneratedDate+".json")
	f, err := os.Create(outpath)
	must(err, "create output")
	enc := json.NewEncoder(f)
	enc.SetIndent("", "  ")
	must(enc.Encode(rm), "encode")
	must(f.Close(), "close output")

	// --- report ---
	fmt.Printf("Committed risk map -> %s\n", outpath)
	fmt.Printf("  %d UTLAs | coverage stratum %s | latest case week %s\n",
		len(entries), *stratum, latestWeek)
	fmt.Printf("  Top 8 by P(large outbreak):\n")
	for _, e := range entries[:8] {
		fmt.Printf("    %-26s s=%.3f  P(R>1)=%.2f  P(large)=%.2f  cluster p50/p90/p99=%d/%d/%d\n",
			e.Name, e.Susceptibility, e.ProbRLocalGt1, e.ProbLargeOutbreak,
			e.ClusterP50, e.ClusterP90, e.ClusterP99)
	}
	if caseLink != nil {
		fmt.Printf("  Susceptibility->cases link (%d, SBI): beta=%.2f CI[%.2f,%.2f] "+
			"-> predicts cases=%v (ESS=%.0f/%d)\n",
			caseLink.ReportYear, caseLink.Beta, caseLink.BetaCI95[0], caseLink.BetaCI95[1],
			caseLink.BetaPredictsCases, caseLink.EffectiveSampleSize, caseLink.NumParticles)
	}
	if calib != nil {
		fmt.Printf("  Reachable-cluster calibration (%d): %.0f susceptibles "+
			"CI[%.0f,%.0f]  importScale=%.1f (ESS=%.0f)\n",
			calib.ReportYear, calib.ReachableCluster, calib.ReachableClusterCI[0],
			calib.ReachableClusterCI[1], calib.ImportScale, calib.EffectiveSampleSize)
	}
	if joint != nil {
		fmt.Printf("  Joint co-occurrence (%d scenarios): UTLAs outbreaking mean=%.1f "+
			"p90=%.0f  national tail p90/p50=%.1fx p99/p50=%.1fx (relative, uncalibrated)\n",
			joint.NScenarios, joint.ExpectedUTLAsOutbreak, joint.P90UTLAsOutbreak,
			joint.NationalTailRatioP90, joint.NationalTailRatioP99)
	}
	if imp != nil {
		fmt.Printf("  National importation pressure (%d European data): %.1fx baseline "+
			"[%.1f, %.1f] (coarse, wide-band)\n",
			imp.Year, imp.RelativeCentral, imp.RelativeLow, imp.RelativeHigh)
		byImp := make([]riskEntry, len(entries))
		copy(byImp, entries)
		sort.Slice(byImp, func(a, b int) bool {
			return byImp[a].ImportationPressure > byImp[b].ImportationPressure
		})
		fmt.Printf("  Top 5 by importation-adjusted pressure:")
		for _, e := range byImp[:5] {
			fmt.Printf(" %s(%.4f)", e.Name, e.ImportationPressure)
		}
		fmt.Println()
	}
}

// calibrateCluster fits the reachable-cluster size to the latest report year's
// observed censored UTLA totals. Returns nil (with a warning) if inputs are
// unavailable, so the forecast falls back to the default.
func calibrateCluster(
	graph *measles.AdjacencyGraph, surf *measles.SusceptibilitySurface,
	receptivity []float64, datadir string,
) *calibrationBlock {
	if receptivity == nil {
		return nil
	}
	cases, err := measles.LoadUTLACases(filepath.Join(datadir, "measles_cases_utla.csv"))
	if err != nil {
		fmt.Printf("  (skip calibration: %v)\n", err)
		return nil
	}
	// Calibrate on the latest COMPLETE calendar year (the current year is partial,
	// so its outbreaks are half-developed and would bias the cluster size low).
	currentYear := time.Now().Year()
	year, maxYear := 0, 0
	for _, c := range cases {
		if c.ReportYear > maxYear {
			maxYear = c.ReportYear
		}
		if c.ReportYear < currentYear && c.ReportYear > year {
			year = c.ReportYear
		}
	}
	if year == 0 {
		year = maxYear
	}
	obs := measles.BuildCensoredCaseObservations(graph, cases, year)
	post, err := measles.CalibrateReachableCluster(
		surf.Smoothed, receptivity, obs.Counts, measles.DefaultSuppressionThreshold,
		15.0 /*r0*/, 0.5 /*dispersion*/,
		inference.UniformPrior{Lo: 1, Hi: 300},
		inference.UniformPrior{Lo: 20, Hi: 3000}, 20000, 1)
	if err != nil {
		fmt.Printf("  (skip calibration: %v)\n", err)
		return nil
	}
	return &calibrationBlock{
		ReachableCluster:   round1(post.ReachableCluster),
		ReachableClusterCI: []float64{round1(post.ReachableClusterCI95[0]), round1(post.ReachableClusterCI95[1])},
		ImportScale:        round2(post.ImportScale),
		R0Used:             post.R0Used,
		EffectiveSampleSize: round1(post.EffectiveSampleSize),
		ReportYear:         year,
	}
}

// loadPopulationVector loads the per-UTLA population aligned to the graph, or nil
// (with a warning) if unavailable — used for susceptible-depletion pools.
func loadPopulationVector(graph *measles.AdjacencyGraph, datadir string) []float64 {
	popMap, err := measles.LoadPopulation(filepath.Join(datadir, "population_utla.csv"))
	if err != nil {
		fmt.Printf("  (population unavailable, using fallback: %v)\n", err)
		return nil
	}
	pop, err := measles.PopulationVector(graph, popMap)
	if err != nil {
		fmt.Printf("  (population unavailable, using fallback: %v)\n", err)
		return nil
	}
	return pop
}

// jointSeedBase is the order-of-magnitude expected number of nationally-seeded
// measles chains per season at baseline European activity. Coarse; scaled by the
// importation index. Documented as an assumption, not a fitted value.
const jointSeedBase = 30.0

// jointCoOccurrence runs the shared-latent joint simulation on the real surface
// and receptivity, using the importation index band as the national seed band,
// and summarises the national-total tail and co-occurring-outbreak count.
func jointCoOccurrence(
	surf *measles.SusceptibilitySurface, receptivity, population []float64,
	imp *importationBlock, reachableCluster float64, seed uint64,
) *jointBlock {
	if imp == nil || receptivity == nil || population == nil {
		return nil
	}
	// Per-UTLA susceptible pool = min(s_i * N_i, reachable cluster), matching the
	// single-UTLA transmission model's calibrated reachable-community bound.
	pool := make([]float64, len(surf.Smoothed))
	for i := range pool {
		p := surf.Smoothed[i] * population[i]
		if p > reachableCluster {
			p = reachableCluster
		}
		pool[i] = p
	}
	const scenarios = 1000
	p := measles.JointParams{
		R0: 15, Dispersion: 0.5,
		SeedLow:  jointSeedBase * imp.RelativeLow,
		SeedHigh: jointSeedBase * imp.RelativeHigh,
		MaxGenerations: 30,
	}
	sc := measles.RunJointEnsemble(surf.Smoothed, receptivity, pool, p, scenarios, seed)

	totals := make([]float64, len(sc))
	outbreakCounts := make([]float64, len(sc))
	for k, s := range sc {
		for _, c := range s.Cumulative {
			totals[k] += c
			if c >= 10 {
				outbreakCounts[k]++
			}
		}
	}
	sort.Float64s(totals)
	sort.Float64s(outbreakCounts)
	mean := 0.0
	for _, v := range outbreakCounts {
		mean += v
	}
	mean /= float64(len(outbreakCounts))

	p50 := quantile(totals, 0.50)
	tailP90, tailP99 := 0.0, 0.0
	if p50 > 0 {
		tailP90 = quantile(totals, 0.90) / p50
		tailP99 = quantile(totals, 0.99) / p50
	}

	return &jointBlock{
		NScenarios:            scenarios,
		SeedBand:              []float64{round1(p.SeedLow), round1(p.SeedHigh)},
		ExpectedUTLAsOutbreak: round2(mean),
		P90UTLAsOutbreak:      quantile(outbreakCounts, 0.90),
		NationalTailRatioP90:  round2(tailP90),
		NationalTailRatioP99:  round2(tailP99),
		Note: "Shared national importation latent -> correlated per-UTLA outbreaks. " +
			"Meaningful outputs are the co-occurrence BREADTH (UTLAs lighting up " +
			"together) and RELATIVE tail fatness. Absolute case counts are NOT " +
			"reported: the pure branching process has no susceptible depletion, so " +
			"it is uncalibrated in magnitude (see README.md).",
	}
}

func quantile(sorted []float64, q float64) float64 {
	if len(sorted) == 0 {
		return 0
	}
	i := int(q * float64(len(sorted)-1))
	return sorted[i]
}

// importationPressure builds the coarse national index for forecastYear and the
// per-UTLA population-share receptivity weights. Returns (nil, nil) if inputs are
// unavailable, so the map still commits without the importation layer.
func importationPressure(
	graph *measles.AdjacencyGraph, datadir string, forecastYear int,
) (*importationBlock, []float64) {
	series, err := measles.LoadEuropeanIncidence(
		filepath.Join(datadir, "european_measles_annual.csv"))
	if err != nil {
		fmt.Printf("  (skip importation: %v)\n", err)
		return nil, nil
	}
	popMap, err := measles.LoadPopulation(filepath.Join(datadir, "population_utla.csv"))
	if err != nil {
		fmt.Printf("  (skip importation: %v)\n", err)
		return nil, nil
	}
	pop, err := measles.PopulationVector(graph, popMap)
	if err != nil {
		fmt.Printf("  (skip importation: %v)\n", err)
		return nil, nil
	}
	base := measles.RecentBaselineYears(series, forecastYear, 8)
	p, err := measles.ComputeImportationPressure(series, forecastYear, base)
	if err != nil {
		fmt.Printf("  (skip importation: %v)\n", err)
		return nil, nil
	}
	block := &importationBlock{
		Year: p.Year, EuropeanCases: p.EuropeanCases,
		BaselineGeomean: round1(p.BaselineGeomean),
		RelativeCentral: round2(p.RelativeCentral),
		RelativeLow:     round2(p.RelativeLow), RelativeHigh: round2(p.RelativeHigh),
		Note: "Coarse annual importation-pressure proxy: WHO European Region measles " +
			"activity relative to baseline. Wide band reflects weak Europe->UK " +
			"identification. Not a weekly forecast (see README.md).",
	}
	return block, measles.ReceptivityWeights(pop)
}

// fitCaseLink runs the seam-2 SBI: load population + the latest report year's
// censored UTLA case counts and fit how strongly susceptibility predicts case
// concentration. Returns nil (with a warning) if the inputs are unavailable, so
// the map still commits.
func fitCaseLink(
	graph *measles.AdjacencyGraph,
	surf *measles.SusceptibilitySurface,
	datadir string,
	seed uint64,
) *caseLinkBlock {
	popMap, err := measles.LoadPopulation(filepath.Join(datadir, "population_utla.csv"))
	if err != nil {
		fmt.Printf("  (skip case-link: %v)\n", err)
		return nil
	}
	pop, err := measles.PopulationVector(graph, popMap)
	if err != nil {
		fmt.Printf("  (skip case-link: %v)\n", err)
		return nil
	}
	cases, err := measles.LoadUTLACases(filepath.Join(datadir, "measles_cases_utla.csv"))
	if err != nil {
		fmt.Printf("  (skip case-link: %v)\n", err)
		return nil
	}
	year := 0
	for _, c := range cases {
		if c.ReportYear > year {
			year = c.ReportYear
		}
	}
	obs := measles.BuildCensoredCaseObservations(graph, cases, year)
	post, err := measles.FitCaseLink(surf.Smoothed, pop, obs.Counts,
		measles.DefaultSuppressionThreshold, 20000, measles.DefaultCaseLinkPriors(), seed)
	if err != nil {
		fmt.Printf("  (skip case-link: %v)\n", err)
		return nil
	}
	return &caseLinkBlock{
		ReportYear: year, Alpha: round2(post.Alpha), Beta: round2(post.Beta),
		BetaStd:  round3(post.BetaStd),
		BetaCI95: []float64{round2(post.BetaCI95[0]), round2(post.BetaCI95[1])},
		BetaPredictsCases:   post.BetaExcludesZero(),
		LogMarginalLik:      round1(post.LogMarginalLik),
		EffectiveSampleSize: round1(post.EffectiveSampleSize), NumParticles: post.NumParticles,
	}
}

// nationalNowcast aggregates Region x week to a national weekly series, orders
// it, and de-biases the trailing reporting-delay-flagged weeks (sub-model C).
func nationalNowcast(recs []measles.RegionWeekCase) ([]nowcastEntry, string) {
	type wk struct {
		year, week int
		date       string
		cases      float64
		inDelay    bool
	}
	byWeek := map[[2]int]*wk{}
	for _, r := range recs {
		k := [2]int{r.Year, r.Epiweek}
		w := byWeek[k]
		if w == nil {
			w = &wk{year: r.Year, week: r.Epiweek, date: r.Date}
			byWeek[k] = w
		}
		w.cases += r.Cases
		if r.InReportingDelay {
			w.inDelay = true
		}
	}
	weeks := make([]*wk, 0, len(byWeek))
	for _, w := range byWeek {
		weeks = append(weeks, w)
	}
	sort.Slice(weeks, func(a, b int) bool {
		if weeks[a].year != weeks[b].year {
			return weeks[a].year < weeks[b].year
		}
		return weeks[a].week < weeks[b].week
	})

	counts := make([]float64, len(weeks))
	tail := 0
	for i, w := range weeks {
		counts[i] = w.cases
		if w.inDelay {
			tail++
		}
	}
	if tail == 0 {
		tail = 1
	}
	// Assumed reporting-delay shape until empirical delays accrue from harvested
	// report snapshots (documented in SOURCES.md).
	delay := measles.ExponentialReportingDelay(2.0, tail+2)
	res := measles.NowcastTail(counts, delay, tail)

	out := make([]nowcastEntry, len(res))
	for i, r := range res {
		out[i] = nowcastEntry{
			AgeWeeks: r.AgeWeeks, Observed: r.Observed,
			Completeness: round3(r.Completeness), Point: round1(r.Point),
			Lower: round1(r.Lower), Upper: round1(r.Upper),
		}
	}
	latest := ""
	if n := len(weeks); n > 0 {
		latest = weeks[n-1].date
	}
	return out, latest
}

func must(err error, ctx string) {
	if err != nil {
		log.Fatalf("%s: %v", ctx, err)
	}
}

func round1(x float64) float64 { return float64(int(x*10+0.5)) / 10 }
func round2(x float64) float64 { return float64(int(x*100+0.5)) / 100 }
func round3(x float64) float64 { return float64(int(x*1000+0.5)) / 1000 }
func round4(x float64) float64 { return float64(int(x*10000+0.5)) / 10000 }
