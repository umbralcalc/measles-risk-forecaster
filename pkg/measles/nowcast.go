// This file is sub-model C (gating check #5): the reporting-lag nowcast. Recent
// onset weeks are right-truncated — not all their cases have been reported yet —
// so the raw tail of the Region x week series understates current activity. This
// de-biases it before it informs the risk surface (a supporting sub-model, not
// the headline; PLAN.md sub-model C).
//
// Method (a generative right-truncation correction, the NobBS / renewal
// nowcasting lineage the plan cites rather than invents): a week whose onset was
// `age` weeks ago has had a fraction Completeness[age] of its eventual cases
// reported. Modelling the reported count as Poisson-thinned by that fraction, the
// as-yet-unreported remainder is Poisson with mean y*(1-F)/F, giving both a point
// estimate and an honest predictive interval that widens for the most recent
// (most truncated) weeks.
package measles

import (
	"math"

	"gonum.org/v1/gonum/stat/distuv"
)

// ReportingDelay describes reporting completeness as a function of weeks since
// onset: Completeness[a] is the expected fraction of an onset-week's eventual
// cases reported within a weeks. It must be non-decreasing and approach 1.
type ReportingDelay struct {
	Completeness []float64
}

// ExponentialReportingDelay builds a saturating completeness curve
// F[a] = 1 - exp(-(a+1)/meanWeeks), truncated at `horizon` weeks (beyond which
// reporting is treated as complete). A reasonable default shape when an empirical
// delay distribution is not yet available from harvested report snapshots.
func ExponentialReportingDelay(meanWeeks float64, horizon int) ReportingDelay {
	c := make([]float64, horizon+1)
	for a := 0; a <= horizon; a++ {
		c[a] = 1 - math.Exp(-float64(a+1)/meanWeeks)
	}
	// Normalise so the last (most complete) entry is treated as fully reported.
	last := c[horizon]
	for a := range c {
		c[a] /= last
	}
	return ReportingDelay{Completeness: c}
}

// completenessAt returns F for a week whose onset was `age` weeks ago (clamped
// to the curve; ages beyond the horizon are treated as fully reported).
func (d ReportingDelay) completenessAt(age int) float64 {
	if age < 0 {
		age = 0
	}
	if age >= len(d.Completeness) {
		return 1
	}
	f := d.Completeness[age]
	if f <= 0 {
		return 1e-6 // guard against division blow-up for an all-truncated week
	}
	return f
}

// NowcastResult is the de-biased estimate for one onset week.
type NowcastResult struct {
	AgeWeeks     int     // weeks since onset (0 = current week)
	Observed     float64 // raw reported count so far
	Completeness float64 // F(age)
	Point        float64 // de-biased point estimate of the eventual total
	Lower, Upper float64 // 95% predictive interval
}

// NowcastTail de-biases the most recent `tail` weeks of an onset-ordered count
// series (index 0 = oldest, last = current week). Older, fully-reported weeks are
// returned unchanged. counts must be in onset order.
func NowcastTail(counts []float64, delay ReportingDelay, tail int) []NowcastResult {
	n := len(counts)
	out := make([]NowcastResult, 0, tail)
	start := n - tail
	if start < 0 {
		start = 0
	}
	for i := start; i < n; i++ {
		age := n - 1 - i // current week has age 0
		f := delay.completenessAt(age)
		y := counts[i]
		// Estimate the eventual total by de-biasing for completeness. Treating the
		// reported count as a Poisson thinning of the true weekly rate lambda,
		// y ~ Poisson(lambda*f), the posterior on lambda under an improper prior is
		// Gamma(shape=y, rate=f) — mean y/f, and a credible interval that widens
		// symmetrically as f shrinks (the most recent, most-truncated weeks). This
		// captures the rate-estimation uncertainty a plug-in y/f point would hide.
		lower, upper := 0.0, 0.0
		if y > 0 {
			post := distuv.Gamma{Alpha: y, Beta: f}
			lower, upper = post.Quantile(0.025), post.Quantile(0.975)
		}
		out = append(out, NowcastResult{
			AgeWeeks:     age,
			Observed:     y,
			Completeness: f,
			Point:        y / f,
			Lower:        lower,
			Upper:        upper,
		})
	}
	return out
}

