// Package measles implements the statistical kernels for the
// measles-risk-forecaster. This file holds the censored count likelihood that
// makes the UTLA case data usable.
//
// Why this exists (gating check #3, see PLAN.md "Validation & honesty rules"):
// the GOV.UK measles epidemiology report suppresses any UTLA with fewer than 10
// confirmed cases — and, crucially, suppression is *row omission*: a suppressed
// UTLA simply does not appear in the table (it is NOT shown as 0, blank, or
// "<10"). See SOURCES.md §2. So a suppressed cell carries the information
// "true count is somewhere in {0,1,...,9}" — interval-censored data, not a zero.
//
// Treating those cells as zero (the naive fix) biases every rate estimate
// downward, exactly where most UTLAs sit. This file provides the *correct*
// likelihood: an observed count contributes its point mass; a suppressed count
// contributes the probability that the count fell anywhere in the censored
// interval, log P(0 <= X <= threshold-1).
package measles

import (
	"math"
	"math/rand/v2"

	"github.com/umbralcalc/stochadex/pkg/inference"
	"github.com/umbralcalc/stochadex/pkg/simulator"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/stat/distuv"
)

// DefaultSuppressionThreshold is the disclosure-control cut applied by the
// UKHSA measles epidemiology tables: counts strictly below this are suppressed.
// A suppressed cell is therefore interval-censored on {0, ..., threshold-1}.
const DefaultSuppressionThreshold = 10

// Suppressed is the sentinel used inside a data vector to mark a cell whose
// true count is unknown because it was below the suppression threshold. We use
// NaN so a suppressed cell can never be confused with a genuine observed count
// (counts are non-negative integers; NaN is unmistakable and propagates loudly
// if accidentally treated as a number).
var Suppressed = math.NaN()

// IsSuppressed reports whether a data value is the censoring sentinel.
func IsSuppressed(x float64) bool { return math.IsNaN(x) }

// CensoredPoissonLogProb returns the log-likelihood contribution of a single
// cell under a Poisson(mean) model with low-count suppression.
//
//   - An observed (non-suppressed) count k contributes log P(X = k).
//   - A suppressed cell contributes log P(0 <= X <= threshold-1), i.e. the log
//     of the Poisson CDF evaluated at threshold-1 (inclusive) — the probability
//     mass of the whole censored interval.
//
// threshold is the suppression cut (counts < threshold are suppressed); pass
// DefaultSuppressionThreshold for the UKHSA tables.
func CensoredPoissonLogProb(value, mean float64, threshold int) float64 {
	dist := distuv.Poisson{Lambda: mean}
	if IsSuppressed(value) {
		// log P(X <= threshold-1). distuv CDF is inclusive of the argument.
		return math.Log(dist.CDF(float64(threshold - 1)))
	}
	return dist.LogProb(value)
}

// CensoredPoissonLogLike sums CensoredPoissonLogProb over a vector of cells,
// each with its own mean. data and mean must be the same length. Suppressed
// entries in data (NaN) are handled as interval-censored on {0,...,threshold-1}.
func CensoredPoissonLogLike(data, mean []float64, threshold int) float64 {
	logLike := 0.0
	for i := range data {
		logLike += CensoredPoissonLogProb(data[i], mean[i], threshold)
	}
	return logLike
}

// FitSharedRateMLE returns the maximum-likelihood Poisson rate for a set of
// cells assumed to share a single rate, under low-count suppression. It is the
// honest estimator the suppression-ablation command contrasts against naive
// zero-filling. Maximisation is a golden-section search on [lo, hi]; the
// censored log-likelihood is unimodal in the rate, so this converges reliably.
func FitSharedRateMLE(data []float64, threshold int, lo, hi float64) float64 {
	ll := func(rate float64) float64 {
		total := 0.0
		for _, v := range data {
			total += CensoredPoissonLogProb(v, rate, threshold)
		}
		return total
	}
	const invPhi = 0.6180339887498949 // 1/phi
	a, b := lo, hi
	c := b - invPhi*(b-a)
	d := a + invPhi*(b-a)
	fc, fd := ll(c), ll(d)
	for i := 0; i < 200 && (b-a) > 1e-9; i++ {
		if fc > fd {
			b, d, fd = d, c, fc
			c = b - invPhi*(b-a)
			fc = ll(c)
		} else {
			a, c, fc = c, d, fd
			d = a + invPhi*(b-a)
			fd = ll(d)
		}
	}
	return 0.5 * (a + b)
}

// CensoredPoissonLikelihoodDistribution adapts the censored Poisson likelihood
// to the stochadex inference.LikelihoodDistribution interface, so it drops into
// the simulation-based inference machinery exactly like the built-in
// PoissonLikelihoodDistribution. The mean is supplied via params or an upstream
// partition, identically to the stochadex Poisson model; the only difference is
// that NaN-marked observations are scored as interval-censored counts.
type CensoredPoissonLikelihoodDistribution struct {
	// Threshold is the suppression cut; zero-value falls back to the default.
	Threshold int
	// Src is the RNG source for GenerateNewSamples; set directly or via SetSeed.
	Src rand.Source

	mean *mat.VecDense
}

func (c *CensoredPoissonLikelihoodDistribution) threshold() int {
	if c.Threshold == 0 {
		return DefaultSuppressionThreshold
	}
	return c.Threshold
}

func (c *CensoredPoissonLikelihoodDistribution) SetSeed(
	partitionIndex int,
	settings *simulator.Settings,
) {
	c.Src = rand.NewPCG(
		settings.Iterations[partitionIndex].Seed,
		settings.Iterations[partitionIndex].Seed,
	)
}

func (c *CensoredPoissonLikelihoodDistribution) SetParams(
	params *simulator.Params,
	partitionIndex int,
	stateHistories []*simulator.StateHistory,
	timestepsHistory *simulator.CumulativeTimestepsHistory,
) {
	c.mean = inference.MeanFromParamsOrPartition(
		params, partitionIndex, stateHistories,
	)
}

func (c *CensoredPoissonLikelihoodDistribution) EvaluateLogLike(
	data []float64,
) float64 {
	threshold := c.threshold()
	logLike := 0.0
	for i := 0; i < c.mean.Len(); i++ {
		logLike += CensoredPoissonLogProb(data[i], c.mean.AtVec(i), threshold)
	}
	return logLike
}

// GenerateNewSamples draws Poisson variates per dimension and applies the same
// suppression rule the real tables apply: any draw below the threshold is
// returned as the Suppressed sentinel. This keeps the generator and the
// likelihood mutually consistent — synthetic data generated here is scored
// correctly by EvaluateLogLike — which is what the recovery test relies on.
func (c *CensoredPoissonLikelihoodDistribution) GenerateNewSamples() []float64 {
	threshold := c.threshold()
	dist := distuv.Poisson{Lambda: 1.0, Src: c.Src}
	samples := make([]float64, 0, c.mean.Len())
	for i := 0; i < c.mean.Len(); i++ {
		dist.Lambda = c.mean.AtVec(i)
		draw := dist.Rand()
		if int(draw) < threshold {
			samples = append(samples, Suppressed)
		} else {
			samples = append(samples, draw)
		}
	}
	return samples
}
