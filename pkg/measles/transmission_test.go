package measles

import (
	"math"
	"math/rand/v2"
	"testing"
)

func TestProbSelfSustainingMonotone(t *testing.T) {
	// Below the herd-immunity threshold (s small) -> 0; well above -> 1; and
	// monotonically increasing in susceptibility between.
	if got := ProbSelfSustaining(0.04, 12, 18); got != 0 {
		t.Errorf("s=0.04 (1/s=25 > R0Max): got %v, want 0", got)
	}
	if got := ProbSelfSustaining(0.20, 12, 18); got != 1 {
		t.Errorf("s=0.20 (1/s=5 < R0Min): got %v, want 1", got)
	}
	prev := -1.0
	for s := 0.0; s <= 0.25; s += 0.01 {
		p := ProbSelfSustaining(s, 12, 18)
		if p < prev-1e-12 {
			t.Errorf("not monotone at s=%.2f: %v < %v", s, p, prev)
		}
		prev = p
	}
	// The discriminating band is s in (1/18, 1/12) ~ (0.056, 0.083).
	if p := ProbSelfSustaining(0.07, 12, 18); p <= 0 || p >= 1 {
		t.Errorf("s=0.07 should be partial, got %v", p)
	}
}

// TestSubcriticalMeanClusterSize validates the branching process against theory:
// for a subcritical process (R_local < 1) with near-Poisson offspring, the mean
// total cluster size seeded by one case is 1/(1-R_local).
func TestSubcriticalMeanClusterSize(t *testing.T) {
	rng := rand.New(rand.NewPCG(20, 26))
	const rLocal = 0.6
	const k = 1e6 // huge dispersion -> offspring ~ Poisson (near-deterministic mixing)
	const n = 40000

	total := 0
	for i := 0; i < n; i++ {
		total += simulateCluster(rLocal, k, 100000, rng)
	}
	got := float64(total) / float64(n)
	want := expectedTotalProgeny(rLocal) // = 1/(1-0.6) = 2.5
	t.Logf("mean cluster size = %.4f, theory 1/(1-m) = %.4f", got, want)
	if math.Abs(got-want)/want > 0.05 {
		t.Errorf("mean cluster size %.4f not within 5%% of %.4f", got, want)
	}
}

// TestSupercriticalProducesLargeOutbreaks checks the tail: a supercritical
// process (R_local > 1) reaches the outbreak cap a substantial fraction of the
// time, while a subcritical one essentially never does.
func TestSupercriticalProducesLargeOutbreaks(t *testing.T) {
	p := TransmissionParams{R0Min: 12, R0Max: 18, Dispersion: 0.5,
		NumSims: 4000, OutbreakCap: 2000}
	rng := rand.New(rand.NewPCG(1, 1))

	// High susceptibility -> R_local well above 1 -> heavy outbreak tail.
	hi := ForecastUTLARisk(0.30, p, rng)
	if hi.ProbLargeOutbreak < 0.3 {
		t.Errorf("high-susceptibility P(large outbreak) = %.3f, expected substantial",
			hi.ProbLargeOutbreak)
	}
	if hi.ClusterP99 < p.OutbreakCap {
		t.Errorf("high-susceptibility p99 cluster = %d, expected to hit cap %d",
			hi.ClusterP99, p.OutbreakCap)
	}

	// Below herd-immunity threshold -> R_local < 1 -> tiny clusters only.
	lo := ForecastUTLARisk(0.04, p, rng)
	if lo.ProbLargeOutbreak > 0.01 {
		t.Errorf("low-susceptibility P(large outbreak) = %.3f, expected ~0",
			lo.ProbLargeOutbreak)
	}
	if lo.ProbRLocalGt1 != 0 {
		t.Errorf("s=0.04 P(R_local>1) = %.3f, want 0", lo.ProbRLocalGt1)
	}
	t.Logf("hi: P(R>1)=%.2f P(large)=%.2f p50=%d p99=%d | lo: P(large)=%.3f p99=%d",
		hi.ProbRLocalGt1, hi.ProbLargeOutbreak, hi.ClusterP50, hi.ClusterP99,
		lo.ProbLargeOutbreak, lo.ClusterP99)
}
