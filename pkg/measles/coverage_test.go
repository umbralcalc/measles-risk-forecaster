package measles

import (
	"math"
	"testing"
)

const coverageCSV = "../../dat/cover_mmr_coverage.csv"

func TestSusceptibilityFromCoverage(t *testing.T) {
	// Full two-dose coverage -> susceptibility = 1 - e2.
	if got := SusceptibilityFromCoverage(1.0, 1.0, 0.93, 0.97); math.Abs(got-0.03) > 1e-12 {
		t.Errorf("full coverage: got %v, want 0.03", got)
	}
	// Zero coverage -> everyone susceptible.
	if got := SusceptibilityFromCoverage(0, 0, 0.93, 0.97); got != 1.0 {
		t.Errorf("zero coverage: got %v, want 1.0", got)
	}
	// One-dose-only cohort uses e1 for the gap above two-dose coverage.
	// c1=0.95, c2=0.85: immune = 0.85*0.97 + 0.10*0.93 = 0.9175 -> s=0.0825.
	if got := SusceptibilityFromCoverage(0.95, 0.85, 0.93, 0.97); math.Abs(got-0.0825) > 1e-9 {
		t.Errorf("partial coverage: got %v, want 0.0825", got)
	}
}

func TestLoadCoverageReal(t *testing.T) {
	recs, err := LoadCoverage(coverageCSV)
	if err != nil {
		t.Fatalf("LoadCoverage: %v (run dat/fetch_coverage.py?)", err)
	}
	// 152 UTLAs x 2 doses x 2 strata x 12 years = 7296 max; expect a big table.
	if len(recs) < 5000 {
		t.Errorf("only %d coverage rows, expected thousands", len(recs))
	}
	codes, strata, doses := map[string]bool{}, map[string]bool{}, map[string]bool{}
	minYear, maxYear := 9999, 0
	for _, r := range recs {
		codes[r.Code] = true
		strata[r.Stratum] = true
		doses[r.Dose] = true
		if r.Year < minYear {
			minYear = r.Year
		}
		if r.Year > maxYear {
			maxYear = r.Year
		}
		if r.CoveragePct < 0 || r.CoveragePct > 100 {
			t.Errorf("coverage %v out of [0,100] for %s", r.CoveragePct, r.Name)
		}
	}
	if len(codes) != 152 {
		t.Errorf("UTLA codes = %d, want 152", len(codes))
	}
	if len(doses) != 2 || len(strata) != 2 {
		t.Errorf("doses=%v strata=%v, want 2 each", doses, strata)
	}
	t.Logf("loaded %d rows, %d UTLAs, years %d-%d", len(recs), len(codes), minYear, maxYear)
}

// TestSusceptibilitySurfaceReal is the real-data end-to-end check #4 deferred:
// build the susceptibility surface from real COVER coverage on the real UTLA
// graph, and verify (a) it is a valid field, (b) the documented COVER-absent
// UTLAs are imputed rather than dropped, and (c) it reproduces the known
// epidemiological signal — London's susceptibility sits above the national
// median (the low-MMR-coverage concentration the 2024 outbreaks tracked).
func TestSusceptibilitySurfaceReal(t *testing.T) {
	g := loadRealGraph(t)
	recs, err := LoadCoverage(coverageCSV)
	if err != nil {
		t.Fatalf("LoadCoverage: %v", err)
	}

	surf, err := BuildSusceptibilitySurface(
		g, recs, "5y", DefaultMMR1Efficacy, DefaultMMR2Efficacy,
		/*tau=*/ 3.0, /*obsPrecision=*/ 50.0,
	)
	if err != nil {
		t.Fatalf("BuildSusceptibilitySurface: %v", err)
	}

	// (a) Valid field: every smoothed value finite and a plausible fraction.
	for i, v := range surf.Smoothed {
		if math.IsNaN(v) || math.IsInf(v, 0) || v < 0 || v > 1 {
			t.Fatalf("node %s: implausible susceptibility %v", g.Names[i], v)
		}
	}

	// (b) The documented COVER-absent UTLAs are imputed, not missing.
	for _, name := range []string{"Rutland", "City of London", "Isles of Scilly"} {
		i := nodeByName(g, name)
		if i < 0 {
			t.Fatalf("%s not in graph", name)
		}
		if !surf.Imputed[i] {
			t.Errorf("%s should be imputed (no COVER row)", name)
		}
		if math.IsNaN(surf.Smoothed[i]) {
			t.Errorf("%s imputed value is NaN", name)
		}
	}
	// Isles of Scilly (island, combined) must equal Cornwall's value.
	ios, corn := nodeByName(g, "Isles of Scilly"), nodeByName(g, "Cornwall")
	if math.Abs(surf.Raw[ios]-surf.Raw[corn]) > 1e-9 {
		t.Errorf("Isles of Scilly raw %.4f != Cornwall %.4f (combination)",
			surf.Raw[ios], surf.Raw[corn])
	}

	// (c) Real-data sanity: London susceptibility > national median.
	londonBoroughs := []string{
		"Hackney", "Newham", "Haringey", "Camden", "Islington", "Brent",
		"Tower Hamlets", "Lambeth", "Lewisham", "Enfield",
	}
	var londonVals []float64
	for _, b := range londonBoroughs {
		if i := nodeByName(g, b); i >= 0 {
			londonVals = append(londonVals, surf.Smoothed[i])
		}
	}
	londonMean := mean(londonVals)
	natMedian := median(surf.Smoothed)
	t.Logf("London mean susceptibility=%.3f  national median=%.3f", londonMean, natMedian)
	if londonMean <= natMedian {
		t.Errorf("expected London (%.3f) above national median (%.3f)", londonMean, natMedian)
	}
}

func mean(xs []float64) float64 {
	s := 0.0
	for _, x := range xs {
		s += x
	}
	return s / float64(len(xs))
}

func median(xs []float64) float64 {
	c := append([]float64(nil), xs...)
	for i := 1; i < len(c); i++ { // insertion sort (n=153, fine)
		for j := i; j > 0 && c[j-1] > c[j]; j-- {
			c[j-1], c[j] = c[j], c[j-1]
		}
	}
	return c[len(c)/2]
}
