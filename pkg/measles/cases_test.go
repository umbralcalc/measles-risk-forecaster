package measles

import (
	"math"
	"testing"
)

const (
	regionWeekCSV = "../../dat/measles_cases_region_week.csv"
	utlaCasesCSV  = "../../dat/measles_cases_utla.csv"
)

func TestLoadRegionWeekCases(t *testing.T) {
	recs, err := LoadRegionWeekCases(regionWeekCSV)
	if err != nil {
		t.Fatalf("LoadRegionWeekCases: %v (run dat/fetch_cases.py?)", err)
	}
	regions := map[string]bool{}
	delayFlagged := 0
	for _, r := range recs {
		regions[r.Region] = true
		if r.InReportingDelay {
			delayFlagged++
		}
		if r.Cases < 0 {
			t.Errorf("negative cases for %s %d-w%d", r.Region, r.Year, r.Epiweek)
		}
	}
	if len(regions) != 9 {
		t.Errorf("regions = %d, want 9", len(regions))
	}
	t.Logf("region-week rows=%d regions=%d delay-flagged weeks=%d",
		len(recs), len(regions), delayFlagged)
}

// TestCensoredCaseObservations checks the row-omission reconstruction: for the
// 2025 report, listed UTLAs carry their count and every other boundary UTLA is
// marked censored (Suppressed). It also asserts every report name matched a
// boundary node (no silent drops).
func TestCensoredCaseObservations(t *testing.T) {
	g := loadRealGraph(t)
	cases, err := LoadUTLACases(utlaCasesCSV)
	if err != nil {
		t.Fatalf("LoadUTLACases: %v", err)
	}

	obs := BuildCensoredCaseObservations(g, cases, 2025)
	if len(obs.Unmatched) > 0 {
		t.Errorf("unmatched report UTLA names (add aliases): %v", obs.Unmatched)
	}

	listed, censored := 0, 0
	for i := range obs.Counts {
		if obs.Censored[i] {
			censored++
			if !math.IsNaN(obs.Counts[i]) {
				t.Errorf("node %s censored but count not Suppressed", g.Names[i])
			}
		} else {
			listed++
			if math.IsNaN(obs.Counts[i]) || obs.Counts[i] < 10 {
				t.Errorf("node %s listed but count = %v (expected >=10)",
					g.Names[i], obs.Counts[i])
			}
		}
	}
	if listed+censored != g.NumNodes() {
		t.Errorf("listed+censored = %d, want %d", listed+censored, g.NumNodes())
	}
	t.Logf("2025: %d UTLAs listed (>=10), %d censored [0,9]", listed, censored)

	// Hackney was the 2025 hotspot — must be listed with a large count.
	hk := nodeByName(g, "Hackney")
	if obs.Censored[hk] || obs.Counts[hk] < 100 {
		t.Errorf("Hackney 2025 count = %v (censored=%v), expected large",
			obs.Counts[hk], obs.Censored[hk])
	}
}
