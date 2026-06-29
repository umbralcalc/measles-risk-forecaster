package measles

import (
	"math/rand/v2"
	"testing"
)

// TestEndToEndRiskDiscrimination exercises the full headline pipeline on real
// data — coverage -> susceptibility surface (A) -> transmission risk (B) — and
// asserts the discrimination the project exists to produce: a low-coverage London
// borough is primed for self-sustaining transmission, while a high-coverage area
// is below the herd-immunity threshold (importations fizzle). This is the
// committed risk map's core claim, pinned against regression.
func TestEndToEndRiskDiscrimination(t *testing.T) {
	g := loadRealGraph(t)
	cov, err := LoadCoverage(coverageCSV)
	if err != nil {
		t.Fatalf("LoadCoverage: %v", err)
	}
	surf, err := BuildSusceptibilitySurface(
		g, cov, "5y", DefaultMMR1Efficacy, DefaultMMR2Efficacy, 3.0, 50.0)
	if err != nil {
		t.Fatalf("BuildSusceptibilitySurface: %v", err)
	}

	tp := DefaultTransmissionParams()
	tp.NumSims = 2000 // lighter for the test
	rng := rand.New(rand.NewPCG(1, 1))

	riskAt := func(name string) UTLARisk {
		i := nodeByName(g, name)
		if i < 0 {
			t.Fatalf("%s not in graph", name)
		}
		return ForecastUTLARisk(surf.Smoothed[i], tp, rng)
	}

	hackney := riskAt("Hackney")       // 2025 hotspot, low coverage
	tyneside := riskAt("South Tyneside") // high coverage

	t.Logf("Hackney:      s=%.3f P(R>1)=%.2f P(large)=%.2f",
		hackney.Susceptibility, hackney.ProbRLocalGt1, hackney.ProbLargeOutbreak)
	t.Logf("S. Tyneside:  s=%.3f P(R>1)=%.2f P(large)=%.2f",
		tyneside.Susceptibility, tyneside.ProbRLocalGt1, tyneside.ProbLargeOutbreak)

	if hackney.ProbRLocalGt1 < 0.99 {
		t.Errorf("Hackney P(R_local>1) = %.2f, expected ~1", hackney.ProbRLocalGt1)
	}
	if hackney.ProbLargeOutbreak < 0.3 {
		t.Errorf("Hackney P(large outbreak) = %.2f, expected substantial",
			hackney.ProbLargeOutbreak)
	}
	if hackney.Susceptibility <= tyneside.Susceptibility {
		t.Errorf("Hackney susceptibility (%.3f) should exceed South Tyneside (%.3f)",
			hackney.Susceptibility, tyneside.Susceptibility)
	}
	// The honest decomposition: even the safer area carries irreducible
	// importation risk (a non-trivial cluster tail), it just rarely self-sustains.
	if tyneside.ProbLargeOutbreak > 0.05 {
		t.Errorf("South Tyneside P(large) = %.2f, expected small", tyneside.ProbLargeOutbreak)
	}
}
