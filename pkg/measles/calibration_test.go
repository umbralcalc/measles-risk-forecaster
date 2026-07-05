package measles

import (
	"math"
	"math/rand/v2"
	"testing"

	"github.com/umbralcalc/stochadex/pkg/inference"
	"gonum.org/v1/gonum/stat/distuv"
)

// directMeanChainSize computes E[final chain size] by direct simulation (grid-free),
// used to generate synthetic data independent of the fitter's interpolation.
func directMeanChainSize(rLocal, disp, cluster float64, nSims int, rng *rand.Rand) float64 {
	total := 0.0
	for i := 0; i < nSims; i++ {
		total += float64(simulateCluster(rLocal, disp, cluster, 1<<30, rng))
	}
	return total / float64(nSims)
}

func calibrationInputs(t *testing.T) (susceptibility, receptivity []float64) {
	g, surf, pop := buildRealSurfacePop(t)
	_ = g
	return surf.Smoothed, ReceptivityWeights(pop)
}

// TestReachableClusterRecovery is the identifiability check: generate censored
// UTLA totals from a known reachable-cluster size (via the compound
// importations x depleting-chains model) and confirm the calibration recovers it.
func TestReachableClusterRecovery(t *testing.T) {
	susceptibility, receptivity := calibrationInputs(t)
	const r0, disp = 15.0, 0.5
	const trueScale, trueCluster = 40.0, 400.0

	// Generate expected totals mu_i, then Poisson counts, then suppress <10.
	rng := rand.New(rand.NewPCG(3, 4))
	n := len(susceptibility)
	counts := make([]float64, n)
	for i := 0; i < n; i++ {
		mc := directMeanChainSize(r0*susceptibility[i], disp, trueCluster, 3000, rng)
		mu := trueScale * receptivity[i] * mc
		c := distuv.Poisson{Lambda: mu, Src: rng}.Rand()
		if c < DefaultSuppressionThreshold {
			counts[i] = Suppressed
		} else {
			counts[i] = c
		}
	}

	post, err := CalibrateReachableCluster(susceptibility, receptivity, counts,
		DefaultSuppressionThreshold, r0, disp,
		inference.UniformPrior{Lo: 1, Hi: 300},
		inference.UniformPrior{Lo: 20, Hi: 3000},
		20000, 7)
	if err != nil {
		t.Fatalf("CalibrateReachableCluster: %v", err)
	}
	t.Logf("recovery: cluster=%.0f±%.0f (true %.0f) CI[%.0f,%.0f]  scale=%.1f ESS=%.0f/%d",
		post.ReachableCluster, post.ReachableClusterStd, trueCluster,
		post.ReachableClusterCI95[0], post.ReachableClusterCI95[1],
		post.ImportScale, post.EffectiveSampleSize, post.NumParticles)

	if post.ReachableClusterCI95[0] > trueCluster || post.ReachableClusterCI95[1] < trueCluster {
		t.Errorf("95%% CI [%.0f,%.0f] does not contain true cluster %.0f",
			post.ReachableClusterCI95[0], post.ReachableClusterCI95[1], trueCluster)
	}
	if r := post.ReachableCluster / trueCluster; r < 0.5 || r > 2.0 {
		t.Errorf("cluster estimate %.0f not within 2x of true %.0f", post.ReachableCluster, trueCluster)
	}
	if post.EffectiveSampleSize < 100 {
		t.Errorf("ESS %.0f too low", post.EffectiveSampleSize)
	}
}

// TestReachableClusterRealData calibrates to the real 2025 UTLA totals and reports
// the posterior reachable-cluster size — the honest, data-driven replacement for
// the default guess.
func TestReachableClusterRealData(t *testing.T) {
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

	post, err := CalibrateReachableCluster(surf.Smoothed, ReceptivityWeights(pop),
		obs.Counts, DefaultSuppressionThreshold, 15.0, 0.5,
		inference.UniformPrior{Lo: 1, Hi: 300},
		inference.UniformPrior{Lo: 20, Hi: 3000},
		20000, 11)
	if err != nil {
		t.Fatalf("CalibrateReachableCluster: %v", err)
	}
	t.Logf("2025 real: reachable cluster=%.0f±%.0f CI[%.0f,%.0f]  importScale=%.1f  ESS=%.0f/%d",
		post.ReachableCluster, post.ReachableClusterStd,
		post.ReachableClusterCI95[0], post.ReachableClusterCI95[1],
		post.ImportScale, post.EffectiveSampleSize, post.NumParticles)

	if post.EffectiveSampleSize < 50 {
		t.Errorf("ESS %.0f too low", post.EffectiveSampleSize)
	}
	if post.ReachableCluster <= 20 || post.ReachableCluster >= 3000 {
		t.Errorf("cluster %.0f railed to a prior edge — check identifiability", post.ReachableCluster)
	}
	if math.IsNaN(post.ReachableCluster) {
		t.Errorf("cluster is NaN")
	}
}
