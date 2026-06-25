package measles

import (
	"math"
	"math/rand/v2"
	"testing"

	"github.com/umbralcalc/stochadex/pkg/simulator"
	"gonum.org/v1/gonum/stat/distuv"
)

// TestSuppressedContributionIsIntervalMass checks the core censoring identity:
// a suppressed cell's log-likelihood must equal the log of the total Poisson
// probability mass over the censored interval {0, ..., threshold-1}, i.e. the
// sum of the individual point masses — and equivalently the log CDF at
// threshold-1. If this is wrong, every downstream estimate is wrong.
func TestSuppressedContributionIsIntervalMass(t *testing.T) {
	const threshold = DefaultSuppressionThreshold
	for _, mean := range []float64{0.3, 2.0, 6.0, 9.5, 20.0} {
		dist := distuv.Poisson{Lambda: mean}
		intervalMass := 0.0
		for k := 0; k < threshold; k++ {
			intervalMass += dist.Prob(float64(k))
		}
		want := math.Log(intervalMass)
		got := CensoredPoissonLogProb(Suppressed, mean, threshold)
		if math.Abs(got-want) > 1e-9 {
			t.Errorf(
				"mean=%v: suppressed logprob = %v, want log(sum P(0..9)) = %v",
				mean, got, want,
			)
		}
	}
}

// TestObservedContributionIsPointMass checks that a non-suppressed count is
// scored as the ordinary Poisson point mass (the censoring must not perturb
// observed cells).
func TestObservedContributionIsPointMass(t *testing.T) {
	dist := distuv.Poisson{Lambda: 7.0}
	for _, k := range []float64{0, 5, 9, 10, 25} {
		want := dist.LogProb(k)
		got := CensoredPoissonLogProb(k, 7.0, DefaultSuppressionThreshold)
		if math.Abs(got-want) > 1e-12 {
			t.Errorf("count=%v: got %v, want %v", k, got, want)
		}
	}
}

// TestCensoredRateRecovery is gating check #3: recover a known rate from
// synthetic, suppression-censored data, and demonstrate that the naive
// zero-fill fix is biased.
//
// We place the true rate at 6.0 — below the suppression threshold — so ~92% of
// cells are suppressed. This is the measles UTLA reality: most UTLAs sit under
// the <10 cut. An honest method must still recover the rate from the small
// observed tail plus the censored-interval mass; the naive method (suppressed →
// 0) cannot.
func TestCensoredRateRecovery(t *testing.T) {
	const (
		trueRate  = 6.0
		nCells    = 5000
		threshold = DefaultSuppressionThreshold
	)
	src := rand.NewPCG(42, 42)
	dist := distuv.Poisson{Lambda: trueRate, Src: src}

	censoredData := make([]float64, nCells) // suppressed cells carry NaN
	zeroFilled := make([]float64, nCells)   // the naive fix: suppressed -> 0
	nSuppressed := 0
	for i := range censoredData {
		draw := dist.Rand()
		if int(draw) < threshold {
			censoredData[i] = Suppressed
			zeroFilled[i] = 0
			nSuppressed++
		} else {
			censoredData[i] = draw
			zeroFilled[i] = draw
		}
	}

	suppFrac := float64(nSuppressed) / float64(nCells)
	if suppFrac < 0.85 {
		t.Fatalf("expected a heavily-suppressed regime, got %.2f suppressed", suppFrac)
	}

	censoredMLE := FitSharedRateMLE(censoredData, threshold, 0.1, 30.0)
	naiveMLE := FitSharedRateMLE(zeroFilled, threshold, 0.1, 30.0)

	t.Logf(
		"suppressed=%.1f%%  trueRate=%.2f  censoredMLE=%.3f  naiveZeroFillMLE=%.3f",
		100*suppFrac, trueRate, censoredMLE, naiveMLE,
	)

	// The honest estimator recovers the true rate.
	if relErr := math.Abs(censoredMLE-trueRate) / trueRate; relErr > 0.05 {
		t.Errorf(
			"censored MLE %.3f not within 5%% of true rate %.2f (rel err %.3f)",
			censoredMLE, trueRate, relErr,
		)
	}

	// The naive zero-fill estimator is biased substantially downward — this is
	// the bias cmd/suppression-ablation will quantify on real data.
	if naiveMLE > 0.6*trueRate {
		t.Errorf(
			"expected naive zero-fill to be badly biased low, got %.3f (true %.2f)",
			naiveMLE, trueRate,
		)
	}
	if naiveMLE >= censoredMLE {
		t.Errorf(
			"naive zero-fill (%.3f) should underestimate vs censored (%.3f)",
			naiveMLE, censoredMLE,
		)
	}
}

// TestCensoredPoissonDistributionInterface exercises the stochadex
// LikelihoodDistribution adapter the way stochadex's own poisson_test does:
// configure a mean via params, generate samples, and score them. It confirms
// the adapter is SBI-ready and that generated suppressed draws round-trip as
// the censoring sentinel and are scored to a finite log-likelihood.
func TestCensoredPoissonDistributionInterface(t *testing.T) {
	d := &CensoredPoissonLikelihoodDistribution{Src: rand.NewPCG(7, 7)}
	params := simulator.NewParams(make(map[string][]float64))
	// One small mean (mostly suppressed) and one large (mostly observed).
	params.Set("mean", []float64{3.0, 40.0})
	d.SetParams(&params, 0, nil, nil)

	sawSuppressed, sawObserved := false, false
	for range 200 {
		s := d.GenerateNewSamples()
		if len(s) != 2 {
			t.Fatalf("expected 2 samples, got %d", len(s))
		}
		if IsSuppressed(s[0]) {
			sawSuppressed = true
		}
		if !IsSuppressed(s[1]) {
			sawObserved = true
		}
		if ll := d.EvaluateLogLike(s); math.IsInf(ll, 0) || math.IsNaN(ll) {
			t.Fatalf("non-finite log-like %v for samples %v", ll, s)
		}
	}
	if !sawSuppressed {
		t.Error("expected some suppressed draws from the low-mean dimension")
	}
	if !sawObserved {
		t.Error("expected some observed draws from the high-mean dimension")
	}
}
