package measles

import (
	"math"
	"math/rand/v2"
	"testing"

	"gonum.org/v1/gonum/stat/distuv"
)

func TestExponentialReportingDelayShape(t *testing.T) {
	d := ExponentialReportingDelay(2.0, 8)
	prev := 0.0
	for a, f := range d.Completeness {
		if f < prev-1e-12 {
			t.Errorf("completeness not monotone at age %d: %v < %v", a, f, prev)
		}
		if f < 0 || f > 1.0000001 {
			t.Errorf("completeness out of [0,1] at age %d: %v", a, f)
		}
		prev = f
	}
	if math.Abs(d.Completeness[len(d.Completeness)-1]-1) > 1e-9 {
		t.Errorf("oldest week not fully reported: %v", d.Completeness[len(d.Completeness)-1])
	}
}

// TestNowcastBeatsRawTruncation: on a right-truncated series the nowcast point
// estimates must be closer to the true totals than the raw reported counts, and
// the most recent (most truncated) week must be inflated the most.
func TestNowcastBeatsRawTruncation(t *testing.T) {
	rng := rand.New(rand.NewPCG(3, 9))
	delay := ExponentialReportingDelay(2.0, 8)

	const weeks = 60
	const tail = 6
	trueTotals := make([]float64, weeks)
	observed := make([]float64, weeks)
	for i := 0; i < weeks; i++ {
		trueTotals[i] = 80 // steady true weekly incidence
		age := weeks - 1 - i
		f := delay.completenessAt(age)
		observed[i] = distuv.Poisson{Lambda: trueTotals[i] * f, Src: rng}.Rand()
	}

	res := NowcastTail(observed, delay, tail)
	rawErr, nowErr := 0.0, 0.0
	for k, r := range res {
		truth := trueTotals[weeks-tail+k]
		rawErr += math.Abs(r.Observed - truth)
		nowErr += math.Abs(r.Point - truth)
	}
	t.Logf("tail=%d  mean|raw-truth|=%.1f  mean|nowcast-truth|=%.1f",
		tail, rawErr/tail, nowErr/tail)
	if nowErr >= rawErr {
		t.Errorf("nowcast (%.1f) did not beat raw truncated counts (%.1f)", nowErr, rawErr)
	}
	last := res[len(res)-1]
	if last.AgeWeeks != 0 || last.Point <= last.Observed {
		t.Errorf("current week not inflated: age=%d point=%.1f observed=%.1f",
			last.AgeWeeks, last.Point, last.Observed)
	}
}

// TestNowcastIntervalCalibration checks the 95% predictive interval is well
// calibrated at the hardest (most truncated) week, age 0: over many independent
// replicates, the truth should fall inside ~95% of the time. A single 6-week
// realisation is far too small a sample to assert this on; this is the proper
// frequentist calibration check.
func TestNowcastIntervalCalibration(t *testing.T) {
	rng := rand.New(rand.NewPCG(101, 202))
	delay := ExponentialReportingDelay(2.0, 8)
	f0 := delay.completenessAt(0)

	const truth = 80.0
	const reps = 5000
	covered := 0
	for r := 0; r < reps; r++ {
		y := distuv.Poisson{Lambda: truth * f0, Src: rng}.Rand()
		// One-week nowcast at age 0: replicate NowcastTail's single-week logic.
		res := NowcastTail([]float64{y}, ReportingDelay{Completeness: []float64{f0}}, 1)[0]
		if truth >= res.Lower && truth <= res.Upper {
			covered++
		}
	}
	cov := float64(covered) / reps
	t.Logf("age-0 completeness f=%.3f  empirical 95%% coverage = %.3f over %d reps",
		f0, cov, reps)
	if cov < 0.92 || cov > 0.98 {
		t.Errorf("interval miscalibrated: coverage %.3f, want ~0.95", cov)
	}
}
