package measles

import (
	"math"
	"testing"
)

func TestAUC(t *testing.T) {
	// Perfect ranking: all positives scored above all negatives -> AUC 1.
	if got := AUC([]float64{4, 3, 2, 1}, []bool{true, true, false, false}); math.Abs(got-1) > 1e-9 {
		t.Errorf("perfect ranking AUC = %v, want 1", got)
	}
	// Reversed -> 0.
	if got := AUC([]float64{1, 2, 3, 4}, []bool{true, true, false, false}); math.Abs(got-0) > 1e-9 {
		t.Errorf("reversed AUC = %v, want 0", got)
	}
	// All tied -> 0.5.
	if got := AUC([]float64{2, 2, 2, 2}, []bool{true, false, true, false}); math.Abs(got-0.5) > 1e-9 {
		t.Errorf("tied AUC = %v, want 0.5", got)
	}
	// Single class -> 0.5 (undefined, safe default).
	if got := AUC([]float64{1, 2, 3}, []bool{false, false, false}); got != 0.5 {
		t.Errorf("single-class AUC = %v, want 0.5", got)
	}
	// Known partial-overlap case: scores 3,2,2,1 labels T,T,F,F.
	// positive ranks: 3->4, 2->avg(2,3)=2.5 ; nPos=2,nNeg=2
	// AUC = (4+2.5 - 2*3/2)/(2*2) = (6.5-3)/4 = 0.875
	if got := AUC([]float64{3, 2, 2, 1}, []bool{true, true, false, false}); math.Abs(got-0.875) > 1e-9 {
		t.Errorf("tie-handling AUC = %v, want 0.875", got)
	}
}

func TestClusterLabels(t *testing.T) {
	g := loadRealGraph(t)
	cases, err := LoadUTLACases(utlaCasesCSV)
	if err != nil {
		t.Fatalf("LoadUTLACases: %v", err)
	}
	lab := ClusterLabels(g, cases, 2025)
	n := 0
	for _, b := range lab {
		if b {
			n++
		}
	}
	if n != 26 { // 2025 report lists 26 UTLAs with >=10 cases
		t.Errorf("2025 clusters = %d, want 26", n)
	}
}
