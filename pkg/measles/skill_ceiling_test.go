package measles

import (
	"math"
	"math/rand/v2"
	"testing"

	"gonum.org/v1/gonum/stat/distuv"
)

// buildRealSurfacePop is a shared fixture: real graph + susceptibility + population.
func buildRealSurfacePop(t *testing.T) (*AdjacencyGraph, *SusceptibilitySurface, []float64) {
	t.Helper()
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
	return g, surf, pop
}

// genCensoredCounts draws counts from the case-link model and applies suppression.
func genCensoredCounts(surf *SusceptibilitySurface, pop []float64, alpha, beta float64, seed uint64) []float64 {
	z := standardisedLogSusceptibility(surf.Smoothed)
	offset := logPopOffsets(pop)
	mu := caseLinkMeans(alpha, beta, z, offset)
	rng := rand.New(rand.NewPCG(seed, seed))
	counts := make([]float64, len(mu))
	for i := range counts {
		c := distuv.Poisson{Lambda: mu[i], Src: rng}.Rand()
		if c < DefaultSuppressionThreshold {
			counts[i] = Suppressed
		} else {
			counts[i] = c
		}
	}
	return counts
}

// TestSpatialSkillDiscriminates: when susceptibility genuinely drives case
// placement, the decomposition must report real skill (McFadden R2 well above 0,
// large likelihood-ratio); when cases follow population only, it must report ~no
// skill — i.e. the "irreducible" fraction is near 1.
func TestSpatialSkillDiscriminates(t *testing.T) {
	_, surf, pop := buildRealSurfacePop(t)

	t.Run("susceptibility-driven", func(t *testing.T) {
		counts := genCensoredCounts(surf, pop, 1.5, 1.2, 11)
		post, err := FitCaseLink(surf.Smoothed, pop, counts,
			DefaultSuppressionThreshold, 20000, DefaultCaseLinkPriors(), 3)
		if err != nil {
			t.Fatal(err)
		}
		sk := SpatialSkillDecomposition(surf.Smoothed, pop, counts,
			DefaultSuppressionThreshold, post.Alpha, post.Beta)
		t.Logf("driven: McFaddenR2=%.3f LR=%.1f irreducible=%.3f (listed=%d censored=%d)",
			sk.McFaddenR2, sk.LRStatistic, sk.IrreducibleFrac, sk.NListed, sk.NCensored)
		if sk.McFaddenR2 < 0.1 || sk.LRStatistic < 20 {
			t.Errorf("expected clear spatial skill, got R2=%.3f LR=%.1f",
				sk.McFaddenR2, sk.LRStatistic)
		}
	})

	t.Run("population-only", func(t *testing.T) {
		counts := genCensoredCounts(surf, pop, 2.0, 0.0, 22) // beta=0: susceptibility irrelevant
		post, err := FitCaseLink(surf.Smoothed, pop, counts,
			DefaultSuppressionThreshold, 20000, DefaultCaseLinkPriors(), 4)
		if err != nil {
			t.Fatal(err)
		}
		sk := SpatialSkillDecomposition(surf.Smoothed, pop, counts,
			DefaultSuppressionThreshold, post.Alpha, post.Beta)
		t.Logf("null: McFaddenR2=%.3f LR=%.1f irreducible=%.3f fittedBeta=%.2f",
			sk.McFaddenR2, sk.LRStatistic, sk.IrreducibleFrac, post.Beta)
		if sk.McFaddenR2 > 0.05 {
			t.Errorf("expected ~no spatial skill when susceptibility irrelevant, got R2=%.3f",
				sk.McFaddenR2)
		}
	})
}

// TestTemporalSkillDiscriminates: an AR(1) series is short-horizon predictable
// (high momentum), white noise is not.
func TestTemporalSkillDiscriminates(t *testing.T) {
	rng := rand.New(rand.NewPCG(9, 9))
	const n = 140

	ar := make([]float64, n)
	l := 3.0
	for i := range ar {
		l = 3.0 + 0.8*(l-3.0) + 0.3*rng.NormFloat64()
		ar[i] = math.Exp(l)
	}
	skAR := TemporalSkillDecomposition(ar)
	t.Logf("AR(1): lag1=%.2f momentumR2=%.2f", skAR.Lag1Autocorr, skAR.MomentumR2)
	if skAR.Lag1Autocorr < 0.5 || skAR.MomentumR2 < 0.3 {
		t.Errorf("AR(1) should be predictable: lag1=%.2f R2=%.2f",
			skAR.Lag1Autocorr, skAR.MomentumR2)
	}

	white := make([]float64, n)
	for i := range white {
		white[i] = math.Exp(3.0 + 0.3*rng.NormFloat64())
	}
	skW := TemporalSkillDecomposition(white)
	t.Logf("white: lag1=%.2f momentumR2=%.2f", skW.Lag1Autocorr, skW.MomentumR2)
	if skW.MomentumR2 > 0.15 {
		t.Errorf("white noise should not be predictable: R2=%.2f", skW.MomentumR2)
	}
}
