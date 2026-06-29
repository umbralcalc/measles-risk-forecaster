package measles

import (
	"math"
	"math/rand/v2"
	"testing"

	"gonum.org/v1/gonum/stat/distuv"
)

const populationCSV = "../../dat/population_utla.csv"

// TestCaseLinkRecoversKnownParams validates seam 2 end to end: generate censored
// case counts from a known (alpha, beta) over the real susceptibility/population
// geometry, then confirm the stochadex SBI recovers them and that beta's credible
// interval excludes zero (the link is detected). This exercises the censored
// likelihood inside inference.ComputePosterior, not a hand-rolled fitter.
func TestCaseLinkRecoversKnownParams(t *testing.T) {
	g := loadRealGraph(t)
	cov, err := LoadCoverage(coverageCSV)
	if err != nil {
		t.Fatalf("LoadCoverage: %v", err)
	}
	surf, err := BuildSusceptibilitySurface(
		g, cov, "5y", DefaultMMR1Efficacy, DefaultMMR2Efficacy, 3.0, 50.0)
	if err != nil {
		t.Fatalf("surface: %v", err)
	}
	popMap, err := LoadPopulation(populationCSV)
	if err != nil {
		t.Fatalf("LoadPopulation: %v (run dat/fetch_population.py?)", err)
	}
	pop, err := PopulationVector(g, popMap)
	if err != nil {
		t.Fatalf("PopulationVector: %v", err)
	}

	const alphaTrue, betaTrue = 1.5, 1.2
	z := standardisedLogSusceptibility(surf.Smoothed)
	offset := logPopOffsets(pop)
	mu := caseLinkMeans(alphaTrue, betaTrue, z, offset)

	// Generate counts and apply the report's <10 row-omission suppression.
	rng := rand.New(rand.NewPCG(5, 5))
	counts := make([]float64, g.NumNodes())
	censoredFrac := 0
	for i := range counts {
		c := distuv.Poisson{Lambda: mu[i], Src: rng}.Rand()
		if c < DefaultSuppressionThreshold {
			counts[i] = Suppressed
			censoredFrac++
		} else {
			counts[i] = c
		}
	}

	post, err := FitCaseLink(surf.Smoothed, pop, counts,
		DefaultSuppressionThreshold, 20000, DefaultCaseLinkPriors(), 99)
	if err != nil {
		t.Fatalf("FitCaseLink: %v", err)
	}
	t.Logf("recovery: alpha=%.2f±%.2f (true %.1f)  beta=%.2f±%.2f (true %.1f)  "+
		"betaCI=[%.2f,%.2f] ESS=%.0f/%d  censored=%d/%d",
		post.Alpha, post.AlphaStd, alphaTrue, post.Beta, post.BetaStd, betaTrue,
		post.BetaCI95[0], post.BetaCI95[1], post.EffectiveSampleSize,
		post.NumParticles, censoredFrac, g.NumNodes())

	if math.Abs(post.Alpha-alphaTrue) > 3*post.AlphaStd+0.2 {
		t.Errorf("alpha %.2f not near true %.1f (sd %.2f)", post.Alpha, alphaTrue, post.AlphaStd)
	}
	if math.Abs(post.Beta-betaTrue) > 3*post.BetaStd+0.2 {
		t.Errorf("beta %.2f not near true %.1f (sd %.2f)", post.Beta, betaTrue, post.BetaStd)
	}
	if !post.BetaExcludesZero() {
		t.Errorf("beta CI [%.2f,%.2f] should exclude 0", post.BetaCI95[0], post.BetaCI95[1])
	}
	if post.EffectiveSampleSize < 100 {
		t.Errorf("ESS %.0f too low — posterior unreliable", post.EffectiveSampleSize)
	}
}

// TestCaseLinkRealData runs the fit on the real 2025 censored UTLA counts and
// asserts the scientific headline: susceptibility positively predicts where
// measles cases concentrated (beta > 0, CI excluding 0).
func TestCaseLinkRealData(t *testing.T) {
	g := loadRealGraph(t)
	cov, err := LoadCoverage(coverageCSV)
	if err != nil {
		t.Fatalf("LoadCoverage: %v", err)
	}
	surf, err := BuildSusceptibilitySurface(
		g, cov, "5y", DefaultMMR1Efficacy, DefaultMMR2Efficacy, 3.0, 50.0)
	if err != nil {
		t.Fatalf("surface: %v", err)
	}
	popMap, err := LoadPopulation(populationCSV)
	if err != nil {
		t.Fatalf("LoadPopulation: %v", err)
	}
	pop, err := PopulationVector(g, popMap)
	if err != nil {
		t.Fatalf("PopulationVector: %v", err)
	}
	cases, err := LoadUTLACases(utlaCasesCSV)
	if err != nil {
		t.Fatalf("LoadUTLACases: %v", err)
	}
	obs := BuildCensoredCaseObservations(g, cases, 2025)

	post, err := FitCaseLink(surf.Smoothed, pop, obs.Counts,
		DefaultSuppressionThreshold, 20000, DefaultCaseLinkPriors(), 7)
	if err != nil {
		t.Fatalf("FitCaseLink: %v", err)
	}
	t.Logf("2025 real fit: alpha=%.2f±%.2f  beta=%.2f±%.2f  betaCI=[%.2f,%.2f]  "+
		"logMargLik=%.1f  ESS=%.0f/%d",
		post.Alpha, post.AlphaStd, post.Beta, post.BetaStd,
		post.BetaCI95[0], post.BetaCI95[1], post.LogMarginalLik,
		post.EffectiveSampleSize, post.NumParticles)

	if post.EffectiveSampleSize < 50 {
		t.Errorf("ESS %.0f too low — tighten priors", post.EffectiveSampleSize)
	}
	if post.Beta <= 0 || !post.BetaExcludesZero() {
		t.Errorf("expected susceptibility to predict cases: beta=%.2f CI=[%.2f,%.2f]",
			post.Beta, post.BetaCI95[0], post.BetaCI95[1])
	}
}
