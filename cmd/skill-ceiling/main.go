// Command skill-ceiling computes and commits the predictable-vs-irreducible
// decomposition (PLAN.md's load-bearing honest figure): how much of where and
// when measles happens is predictable from the susceptibility surface, versus
// irreducible. Writes data/skill_ceiling.json.
//
// Usage: go run ./cmd/skill-ceiling [-datadir dat] [-outdir data] [-year 2025]
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"sort"

	"github.com/umbralcalc/measles-risk-forecaster/pkg/measles"
)

func main() {
	datadir := flag.String("datadir", "dat", "directory of ingested CSVs")
	outdir := flag.String("outdir", "data", "output directory")
	year := flag.Int("year", 2025, "report year of UTLA case counts to score against")
	flag.Parse()

	graph, err := measles.LoadAdjacency(
		filepath.Join(*datadir, "utla_nodes.csv"),
		filepath.Join(*datadir, "utla_adjacency.csv"))
	must(err, "adjacency")
	cov, err := measles.LoadCoverage(filepath.Join(*datadir, "cover_mmr_coverage.csv"))
	must(err, "coverage")
	surf, err := measles.BuildSusceptibilitySurface(
		graph, cov, "5y", measles.DefaultMMR1Efficacy, measles.DefaultMMR2Efficacy, 3.0, 50.0)
	must(err, "surface")
	popMap, err := measles.LoadPopulation(filepath.Join(*datadir, "population_utla.csv"))
	must(err, "population")
	pop, err := measles.PopulationVector(graph, popMap)
	must(err, "population vector")
	cases, err := measles.LoadUTLACases(filepath.Join(*datadir, "measles_cases_utla.csv"))
	must(err, "utla cases")
	regionWeek, err := measles.LoadRegionWeekCases(
		filepath.Join(*datadir, "measles_cases_region_week.csv"))
	must(err, "region-week cases")

	// --- Spatial: fit susceptibility->cases, then decompose skill vs irreducible.
	obs := measles.BuildCensoredCaseObservations(graph, cases, *year)
	post, err := measles.FitCaseLink(surf.Smoothed, pop, obs.Counts,
		measles.DefaultSuppressionThreshold, 40000, measles.DefaultCaseLinkPriors(), 1)
	must(err, "fit case link")
	spatial := measles.SpatialSkillDecomposition(surf.Smoothed, pop, obs.Counts,
		measles.DefaultSuppressionThreshold, post.Alpha, post.Beta)

	// --- Temporal: national weekly series (drop right-truncated tail), decompose.
	weekly := nationalWeeklySeries(regionWeek)
	temporal := measles.TemporalSkillDecomposition(weekly)

	out := map[string]any{
		"report_year": *year,
		"spatial": map[string]any{
			"question":            "Does the susceptibility surface predict WHERE cases concentrate, beyond population?",
			"mcfadden_r2":         round3(spatial.McFaddenR2),
			"irreducible_frac":    round3(spatial.IrreducibleFrac),
			"lr_statistic":        round1(spatial.LRStatistic),
			"beta":                round2(post.Beta),
			"beta_ci95":           []float64{round2(post.BetaCI95[0]), round2(post.BetaCI95[1])},
			"n_listed":            spatial.NListed,
			"n_censored":          spatial.NCensored,
			"interpretation": "susceptibility explains this McFadden-R2 share of case-concentration " +
				"deviance beyond population; the remainder is irreducible spark placement",
		},
		"temporal": map[string]any{
			"question":      "How predictable is national weekly activity short-horizon?",
			"n_weeks":       temporal.NWeeks,
			"lag1_autocorr": round3(temporal.Lag1Autocorr),
			"momentum_r2":   round3(temporal.MomentumR2),
			"caveat":        temporal.Caveat,
		},
		"headline": fmt.Sprintf(
			"Susceptibility explains ~%.0f%% of where cases concentrate (McFadden R2=%.2f); "+
				"~%.0f%% is irreducible spark placement. National activity is %.0f%% "+
				"short-horizon predictable, but that is transmission momentum, not importation skill.",
			100*spatial.McFaddenR2, spatial.McFaddenR2,
			100*spatial.IrreducibleFrac, 100*temporal.MomentumR2),
	}

	must(os.MkdirAll(*outdir, 0o755), "mkdir outdir")
	outpath := filepath.Join(*outdir, "skill_ceiling.json")
	f, err := os.Create(outpath)
	must(err, "create output")
	enc := json.NewEncoder(f)
	enc.SetIndent("", "  ")
	must(enc.Encode(out), "encode")
	must(f.Close(), "close")

	fmt.Printf("Skill-ceiling decomposition -> %s\n", outpath)
	fmt.Printf("  SPATIAL (%d): susceptibility McFadden R2=%.3f  irreducible=%.3f  "+
		"LR=%.0f  beta=%.2f CI[%.2f,%.2f]  (%d listed, %d censored)\n",
		*year, spatial.McFaddenR2, spatial.IrreducibleFrac, spatial.LRStatistic,
		post.Beta, post.BetaCI95[0], post.BetaCI95[1], spatial.NListed, spatial.NCensored)
	fmt.Printf("  TEMPORAL: %d weeks  lag1=%.2f  momentum R2=%.2f\n",
		temporal.NWeeks, temporal.Lag1Autocorr, temporal.MomentumR2)
	fmt.Printf("  CAVEAT: %s\n", temporal.Caveat)
	fmt.Printf("  HEADLINE: %s\n", out["headline"])
}

// nationalWeeklySeries aggregates Region x week to a national onset-ordered
// series, dropping trailing reporting-delay-flagged (right-truncated) weeks.
func nationalWeeklySeries(recs []measles.RegionWeekCase) []float64 {
	type wk struct {
		year, week int
		cases      float64
		inDelay    bool
	}
	byWeek := map[[2]int]*wk{}
	for _, r := range recs {
		k := [2]int{r.Year, r.Epiweek}
		w := byWeek[k]
		if w == nil {
			w = &wk{year: r.Year, week: r.Epiweek}
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
	for len(weeks) > 0 && weeks[len(weeks)-1].inDelay {
		weeks = weeks[:len(weeks)-1]
	}
	out := make([]float64, len(weeks))
	for i, w := range weeks {
		out[i] = w.cases
	}
	return out
}

func must(err error, ctx string) {
	if err != nil {
		log.Fatalf("%s: %v", ctx, err)
	}
}
func round1(x float64) float64 { return float64(int(x*10+0.5)) / 10 }
func round2(x float64) float64 { return float64(int(x*100+0.5)) / 100 }
func round3(x float64) float64 { return float64(int(x*1000+0.5)) / 1000 }
