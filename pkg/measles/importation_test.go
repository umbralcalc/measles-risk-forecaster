package measles

import (
	"math"
	"testing"
)

const europeanCSV = "../../dat/european_measles_annual.csv"

func TestImportationPressureRealData(t *testing.T) {
	series, err := LoadEuropeanIncidence(europeanCSV)
	if err != nil {
		t.Fatalf("LoadEuropeanIncidence: %v (run dat/fetch_european_measles.py?)", err)
	}

	// 2024 was a major European surge -> index well above 1, with a wide band.
	base := RecentBaselineYears(series, 2024, 8)
	p, err := ComputeImportationPressure(series, 2024, base)
	if err != nil {
		t.Fatalf("ComputeImportationPressure: %v", err)
	}
	t.Logf("2024: E=%.0f baseline=%.0f  relative central=%.2f band[%.2f,%.2f]  base years=%v",
		p.EuropeanCases, p.BaselineGeomean, p.RelativeCentral, p.RelativeLow, p.RelativeHigh, base)

	if p.RelativeCentral < 2 {
		t.Errorf("2024 index central=%.2f, expected well above baseline", p.RelativeCentral)
	}
	if !(p.RelativeLow < p.RelativeCentral && p.RelativeCentral < p.RelativeHigh) {
		t.Errorf("band not ordered: [%.2f, %.2f, %.2f]", p.RelativeLow, p.RelativeCentral, p.RelativeHigh)
	}
	// The band must be genuinely wide (elasticity uncertainty), not cosmetic.
	if p.RelativeHigh/p.RelativeLow < 1.8 {
		t.Errorf("band too narrow (%.2f/%.2f) — should reflect weak identification",
			p.RelativeHigh, p.RelativeLow)
	}

	// 2021 (pandemic collapse) -> index far below 1.
	p21, err := ComputeImportationPressure(series, 2021, base)
	if err != nil {
		t.Fatalf("2021: %v", err)
	}
	if p21.RelativeCentral > 0.2 {
		t.Errorf("2021 index=%.3f, expected near zero (collapse)", p21.RelativeCentral)
	}
	t.Logf("2021: relative central=%.4f (pandemic collapse)", p21.RelativeCentral)

	// Annual-lag fallback: a future year with no data uses the latest available.
	pFuture, err := ComputeImportationPressure(series, 2030, base)
	if err != nil {
		t.Fatalf("future fallback: %v", err)
	}
	if pFuture.Year > 2024 {
		t.Errorf("expected fallback to latest available year, got %d", pFuture.Year)
	}
}

// TestImportationCoOccurrence checks the mechanism: raising the national index
// scales every UTLA's expected pressure together (co-occurrence), and the score
// concentrates where receptivity (population) and amplification (susceptibility)
// are both high.
func TestImportationCoOccurrence(t *testing.T) {
	receptivity := []float64{0.1, 0.5, 0.4}
	pSelf := []float64{0.2, 0.9, 0.1}

	low := ImportationAdjustedPressure(1.0, receptivity, pSelf)
	high := ImportationAdjustedPressure(3.0, receptivity, pSelf)

	for i := range low {
		if math.Abs(high[i]-3*low[i]) > 1e-9 {
			t.Errorf("co-occurrence broken at %d: high=%.4f low=%.4f", i, high[i], low[i])
		}
	}
	// UTLA 1 (high receptivity AND high amplification) must dominate.
	if !(low[1] > low[0] && low[1] > low[2]) {
		t.Errorf("expected UTLA 1 to dominate: %v", low)
	}
}

func TestReceptivityWeightsSumToOne(t *testing.T) {
	w := ReceptivityWeights([]float64{100, 300, 600})
	sum := 0.0
	for _, x := range w {
		sum += x
	}
	if math.Abs(sum-1) > 1e-12 {
		t.Errorf("weights sum = %.6f, want 1", sum)
	}
	if math.Abs(w[2]-0.6) > 1e-12 {
		t.Errorf("largest weight = %.3f, want 0.6", w[2])
	}
}
