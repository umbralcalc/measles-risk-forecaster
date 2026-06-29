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

type riskMap struct {
	GeneratedDate string         `json:"generated_date"`
	DataVintage   map[string]any `json:"data_vintage"`
	Model         map[string]any `json:"model"`
	Attribution   string         `json:"attribution"`
	SusceptibilityCaseLink *caseLinkBlock `json:"susceptibility_case_link"`
	UTLAs         []riskEntry    `json:"utlas"`
	NationalNowcast []nowcastEntry `json:"national_nowcast"`
}

const attribution = "Contains public sector information from UKHSA (measles " +
	"epidemiology and COVER vaccination coverage) and the Office for National " +
	"Statistics (geographies and population), licensed under the Open Government " +
	"Licence v3.0. Contains OS data © Crown copyright and database right 2024. " +
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

	// --- seam 2: fit the susceptibility->cases link (SBI) ---
	caseLink := fitCaseLink(graph, surf, *datadir, *seed)

	// --- sub-model B: per-UTLA transmission risk ---
	tp := measles.DefaultTransmissionParams()
	entries := make([]riskEntry, graph.NumNodes())
	for i := 0; i < graph.NumNodes(); i++ {
		r := measles.ForecastUTLARisk(surf.Smoothed[i], tp, rng)
		entries[i] = riskEntry{
			Code: graph.Codes[i], Name: graph.Names[i],
			Susceptibility: round3(surf.Smoothed[i]), Imputed: surf.Imputed[i],
			ProbRLocalGt1: round3(r.ProbRLocalGt1), MedianRLocal: round2(r.MedianRLocal),
			ClusterP50: r.ClusterP50, ClusterP90: r.ClusterP90, ClusterP99: r.ClusterP99,
			ProbLargeOutbreak: round3(r.ProbLargeOutbreak),
		}
	}
	sort.Slice(entries, func(a, b int) bool {
		if entries[a].ProbLargeOutbreak != entries[b].ProbLargeOutbreak {
			return entries[a].ProbLargeOutbreak > entries[b].ProbLargeOutbreak
		}
		return entries[a].Susceptibility > entries[b].Susceptibility
	})

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
			"num_sims": tp.NumSims, "outbreak_cap": tp.OutbreakCap,
			"car_tau": 3.0, "mmr1_efficacy": measles.DefaultMMR1Efficacy,
			"mmr2_efficacy": measles.DefaultMMR2Efficacy, "seed": *seed,
		},
		Attribution: attribution, SusceptibilityCaseLink: caseLink,
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
