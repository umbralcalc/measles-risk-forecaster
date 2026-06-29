// This file is seam 2 of the stochadex integration: fitting the
// susceptibility->cases observation model by simulation-based inference, using
// the stochadex inference machinery (inference.ComputePosterior + Prior types)
// and the censored Poisson likelihood (#3) on the real, row-omission-censored
// UTLA case counts.
//
// Model. For UTLA i with population pop_i and susceptibility s_i, the cumulative
// case count is
//
//	cases_i ~ Poisson(mu_i),   log mu_i = alpha + log(pop_i / popBar) + beta * z_i,
//
// where z_i is standardised log-susceptibility. beta is the headline inferential
// output: how strongly the susceptibility surface predicts where cases actually
// land, with full posterior uncertainty. UTLAs the report suppressed (<10 cases)
// enter the likelihood as interval-censored [0,9], not as zero.
//
// Inference. We importance-sample (alpha, beta) from the prior and weight each
// particle by the censored-data likelihood; since the proposal IS the prior, the
// importance weights are exactly the likelihood, which is what
// inference.ComputePosterior consumes. Effective sample size is reported so the
// posterior's reliability is explicit, not assumed.
package measles

import (
	"fmt"
	"math"
	"math/rand/v2"

	"github.com/umbralcalc/stochadex/pkg/inference"
)

// CaseLinkPosterior summarises the fitted susceptibility->cases link.
type CaseLinkPosterior struct {
	Alpha, Beta         float64    // posterior means
	AlphaStd, BetaStd   float64    // posterior standard deviations
	BetaCI95            [2]float64 // 95% credible interval for beta (weighted quantiles)
	LogMarginalLik      float64
	EffectiveSampleSize float64 // 1 / sum(w_i^2); compare to NumParticles
	NumParticles        int
}

// BetaExcludesZero reports whether the 95% credible interval for beta lies wholly
// above zero — i.e. susceptibility significantly predicts case concentration.
func (c *CaseLinkPosterior) BetaExcludesZero() bool { return c.BetaCI95[0] > 0 }

// standardisedLogSusceptibility returns z_i = (log s_i - mean) / sd across UTLAs.
func standardisedLogSusceptibility(s []float64) []float64 {
	n := len(s)
	logs := make([]float64, n)
	mean := 0.0
	for i, v := range s {
		logs[i] = math.Log(math.Max(v, 1e-9))
		mean += logs[i]
	}
	mean /= float64(n)
	varSum := 0.0
	for _, v := range logs {
		varSum += (v - mean) * (v - mean)
	}
	sd := math.Sqrt(varSum / float64(n))
	if sd == 0 {
		sd = 1
	}
	z := make([]float64, n)
	for i, v := range logs {
		z[i] = (v - mean) / sd
	}
	return z
}

// logPopOffsets returns log(pop_i / popBar), the centred population offset.
func logPopOffsets(pop []float64) []float64 {
	mean := 0.0
	for _, p := range pop {
		mean += p
	}
	mean /= float64(len(pop))
	out := make([]float64, len(pop))
	for i, p := range pop {
		out[i] = math.Log(p / mean)
	}
	return out
}

// caseLinkMeans computes mu_i = exp(alpha + offset_i + beta*z_i).
func caseLinkMeans(alpha, beta float64, z, offset []float64) []float64 {
	mu := make([]float64, len(z))
	for i := range z {
		mu[i] = math.Exp(alpha + offset[i] + beta*z[i])
	}
	return mu
}

// CaseLinkPriors are the priors over (alpha, beta). alpha is the log expected
// count at mean population and mean susceptibility; beta the susceptibility
// elasticity. Defaults are weakly informative and centred on plausible values.
type CaseLinkPriors struct {
	Alpha, Beta inference.UniformPrior
}

// DefaultCaseLinkPriors: alpha ~ U(-2, 6), beta ~ U(-3, 5).
func DefaultCaseLinkPriors() CaseLinkPriors {
	return CaseLinkPriors{
		Alpha: inference.UniformPrior{Lo: -2, Hi: 6},
		Beta:  inference.UniformPrior{Lo: -3, Hi: 5},
	}
}

// caseLinkRounds is the number of iterated importance-sampling rounds: round 1
// proposes from the prior, later rounds from a normal centred on the running
// posterior (population Monte Carlo). The data is highly informative, so a single
// prior round has a tiny effective sample size; refining the proposal fixes that.
const caseLinkRounds = 4

// FitCaseLink runs the SBI. susceptibility, population and censoredCounts must be
// aligned to the same UTLA order; censoredCounts carries observed counts or
// Suppressed (NaN) for row-omitted UTLAs. threshold is the suppression cut (10).
// numParticles is the per-round particle count.
func FitCaseLink(
	susceptibility, population, censoredCounts []float64,
	threshold, numParticles int,
	priors CaseLinkPriors,
	seed uint64,
) (*CaseLinkPosterior, error) {
	n := len(susceptibility)
	if len(population) != n || len(censoredCounts) != n {
		return nil, fmt.Errorf("susceptibility, population, counts must align (%d)", n)
	}
	z := standardisedLogSusceptibility(susceptibility)
	offset := logPopOffsets(population)
	rng := rand.New(rand.NewPCG(seed, seed))

	logPriorConst := priors.Alpha.LogPDF(0.5*(priors.Alpha.Lo+priors.Alpha.Hi)) +
		priors.Beta.LogPDF(0.5*(priors.Beta.Lo+priors.Beta.Hi)) // uniform: constant in support

	var res *inference.SMCResult
	var particles [][]float64
	var proposal *mvn2D // nil on round 1 (propose from prior)

	for round := 0; round < caseLinkRounds; round++ {
		particles = make([][]float64, numParticles)
		logW := make([]float64, numParticles)
		for p := 0; p < numParticles; p++ {
			var a, b, logProp float64
			if proposal == nil {
				a, b = priors.Alpha.Sample(rng), priors.Beta.Sample(rng)
			} else {
				// Propose from the running posterior, rejecting outside support.
				for {
					a, b = proposal.sample(rng)
					if priors.Alpha.InSupport(a) && priors.Beta.InSupport(b) {
						break
					}
				}
				logProp = proposal.logPDF(a, b)
			}
			particles[p] = []float64{a, b}
			mu := caseLinkMeans(a, b, z, offset)
			loglik := CensoredPoissonLogLike(censoredCounts, mu, threshold)
			if proposal == nil {
				// Proposal IS the prior: importance weight is the likelihood.
				logW[p] = loglik
			} else {
				// Importance weight = likelihood * prior / proposal.
				logW[p] = loglik + logPriorConst - logProp
			}
		}
		res = inference.ComputePosterior([]string{"alpha", "beta"}, particles, logW, nil)
		// Build next round's proposal: posterior, covariance inflated for coverage.
		next, ok := newMVN2D(res.PosteriorMean, res.PosteriorCov, 2.5)
		if !ok {
			break // degenerate covariance — keep the current posterior
		}
		proposal = next
	}

	sumSq := 0.0
	betaValues := make([]float64, len(particles))
	for i, w := range res.Weights {
		sumSq += w * w
		betaValues[i] = particles[i][1]
	}
	ess := 0.0
	if sumSq > 0 {
		ess = 1.0 / sumSq
	}
	ci := inference.WeightedQuantiles(betaValues, res.Weights, []float64{0.025, 0.975})

	return &CaseLinkPosterior{
		Alpha: res.PosteriorMean[0], Beta: res.PosteriorMean[1],
		AlphaStd: res.PosteriorStd[0], BetaStd: res.PosteriorStd[1],
		BetaCI95:            [2]float64{ci[0], ci[1]},
		LogMarginalLik:      res.LogMarginalLik,
		EffectiveSampleSize: ess, NumParticles: numParticles,
	}, nil
}

// mvn2D is a 2-dimensional multivariate normal proposal (the model has exactly
// two parameters), used for the population-Monte-Carlo refinement.
type mvn2D struct {
	mean         [2]float64
	l00, l10, l11 float64 // lower Cholesky factor (cov = L Lᵀ)
	invDet       float64
	c00, c01, c11 float64 // (inflated) covariance entries
	logNorm      float64  // -log(2π) - 0.5 log det
}

func newMVN2D(mean, cov []float64, inflate float64) (*mvn2D, bool) {
	c00, c01, c11 := cov[0]*inflate, cov[1]*inflate, cov[3]*inflate
	det := c00*c11 - c01*c01
	if det <= 0 || c00 <= 0 {
		return nil, false
	}
	l00 := math.Sqrt(c00)
	l10 := c01 / l00
	rem := c11 - l10*l10
	if rem <= 0 {
		return nil, false
	}
	return &mvn2D{
		mean: [2]float64{mean[0], mean[1]},
		l00:  l00, l10: l10, l11: math.Sqrt(rem),
		invDet: 1 / det, c00: c00, c01: c01, c11: c11,
		logNorm: -math.Log(2*math.Pi) - 0.5*math.Log(det),
	}, true
}

func (m *mvn2D) sample(rng *rand.Rand) (float64, float64) {
	z0, z1 := rng.NormFloat64(), rng.NormFloat64()
	return m.mean[0] + m.l00*z0, m.mean[1] + m.l10*z0 + m.l11*z1
}

func (m *mvn2D) logPDF(a, b float64) float64 {
	da, db := a-m.mean[0], b-m.mean[1]
	// quadratic form dᵀ Σ⁻¹ d with Σ⁻¹ = 1/det [[c11,-c01],[-c01,c00]]
	q := m.invDet * (m.c11*da*da - 2*m.c01*da*db + m.c00*db*db)
	return m.logNorm - 0.5*q
}
