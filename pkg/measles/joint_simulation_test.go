package measles

import (
	"math"
	"testing"
)

func jointTestParams(low, high float64) JointParams {
	return JointParams{
		R0: 15, Dispersion: 0.5,
		SeedLow: low, SeedHigh: high, MaxGenerations: 30,
	}
}

// uniformPool builds a per-UTLA susceptible pool s_i * N for a synthetic N.
func uniformPool(susceptibility []float64, population float64) []float64 {
	pool := make([]float64, len(susceptibility))
	for i, s := range susceptibility {
		pool[i] = s * population
	}
	return pool
}

// nationalTotals sums each scenario's per-UTLA cumulative cases.
func nationalTotals(sc []JointScenario) []float64 {
	out := make([]float64, len(sc))
	for k, s := range sc {
		for _, c := range s.Cumulative {
			out[k] += c
		}
	}
	return out
}

func variance(x []float64) float64 {
	m, n := 0.0, float64(len(x))
	for _, v := range x {
		m += v
	}
	m /= n
	s := 0.0
	for _, v := range x {
		s += (v - m) * (v - m)
	}
	return s / n
}

// meanOutbreakCorrelation is the mean pairwise Pearson correlation of per-UTLA
// outbreak indicators (>=5 cases) across scenarios — the co-occurrence metric.
func meanOutbreakCorrelation(sc []JointScenario, nUTLA int) float64 {
	ind := make([][]float64, nUTLA) // [utla][scenario]
	for i := range ind {
		ind[i] = make([]float64, len(sc))
	}
	for k, s := range sc {
		for i := 0; i < nUTLA; i++ {
			if s.Cumulative[i] >= 5 {
				ind[i][k] = 1
			}
		}
	}
	sum, pairs := 0.0, 0
	for i := 0; i < nUTLA; i++ {
		for j := i + 1; j < nUTLA; j++ {
			if c, ok := pearson(ind[i], ind[j]); ok {
				sum += c
				pairs++
			}
		}
	}
	if pairs == 0 {
		return 0
	}
	return sum / float64(pairs)
}

func pearson(a, b []float64) (float64, bool) {
	n := float64(len(a))
	ma, mb := 0.0, 0.0
	for i := range a {
		ma += a[i]
		mb += b[i]
	}
	ma /= n
	mb /= n
	var sab, saa, sbb float64
	for i := range a {
		da, db := a[i]-ma, b[i]-mb
		sab += da * db
		saa += da * da
		sbb += db * db
	}
	if saa == 0 || sbb == 0 {
		return 0, false // a constant series has no correlation
	}
	return sab / math.Sqrt(saa*sbb), true
}

// TestJointSimulationRuns is the integration check: the 2-partition joint model
// runs under the stochadex coordinator and yields sane per-UTLA totals.
func TestJointSimulationRuns(t *testing.T) {
	n := 12
	suscept := make([]float64, n)
	recept := make([]float64, n)
	for i := range suscept {
		suscept[i] = 0.12
		recept[i] = 1.0 / float64(n)
	}
	pool := uniformPool(suscept, 100000)
	sc := RunJointEnsemble(suscept, recept, pool, jointTestParams(40, 40), 20, 1)
	if len(sc) != 20 {
		t.Fatalf("got %d scenarios, want 20", len(sc))
	}
	for _, s := range sc {
		if len(s.Cumulative) != n {
			t.Fatalf("cumulative width %d, want %d", len(s.Cumulative), n)
		}
		for _, c := range s.Cumulative {
			if c < 0 {
				t.Errorf("negative cumulative %v", c)
			}
		}
	}
}

// TestJointCoOccurrence is the whole point: a shared national importation latent
// makes UTLA outbreaks CO-OCCUR (positively correlated) and fattens the
// national-total tail, versus an independent baseline with the national total
// fixed at the same central level.
func TestJointCoOccurrence(t *testing.T) {
	n := 15
	suscept := make([]float64, n)
	recept := make([]float64, n)
	for i := range suscept {
		suscept[i] = 0.13 // moderately supercritical: outbreaks are seeding-limited
		recept[i] = 1.0 / float64(n)
	}
	const reps = 600

	pool := uniformPool(suscept, 100000)

	// Shared latent: national total drawn log-uniform over a wide band.
	low, high := 10.0, 300.0
	drawn := RunJointEnsemble(suscept, recept, pool, jointTestParams(low, high), reps, 100)

	// Independent baseline: national total fixed at the band's geometric mean.
	mid := math.Exp(0.5 * (math.Log(low) + math.Log(high)))
	fixed := RunJointEnsemble(suscept, recept, pool, jointTestParams(mid, mid), reps, 500)

	corrDrawn := meanOutbreakCorrelation(drawn, n)
	corrFixed := meanOutbreakCorrelation(fixed, n)
	varDrawn := variance(nationalTotals(drawn))
	varFixed := variance(nationalTotals(fixed))

	t.Logf("shared-latent : mean outbreak corr=%.3f  national-total var=%.3g", corrDrawn, varDrawn)
	t.Logf("independent   : mean outbreak corr=%.3f  national-total var=%.3g", corrFixed, varFixed)

	if corrDrawn < 0.1 {
		t.Errorf("shared latent should induce co-occurrence, corr=%.3f", corrDrawn)
	}
	if corrDrawn <= corrFixed+0.05 {
		t.Errorf("shared latent corr (%.3f) should exceed independent (%.3f)", corrDrawn, corrFixed)
	}
	if varDrawn <= varFixed {
		t.Errorf("shared latent should fatten national-total tail: var %.3g vs %.3g",
			varDrawn, varFixed)
	}
}
