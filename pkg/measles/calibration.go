// This file calibrates the reachable-cluster size (the well-mixed susceptible
// community reachable from one introduction) to the observed UTLA outbreak-size
// data, so the transmission model's absolute magnitudes stop being illustrative.
//
// The key idea: an observed UTLA annual total is a COMPOUND of importations x
// depleting chains, so its expected value factorises as
//
//	mu_i = importScale * receptivity_i * meanChainSize(R_local,i, reachableCluster),
//
// where meanChainSize is the mean final size of one susceptible-depleting chain.
// For supercritical (high-susceptibility) UTLAs meanChainSize saturates near the
// reachable-cluster size; for subcritical UTLAs it stays small and is independent
// of it. So the CONTRAST between high- and low-susceptibility totals identifies
// the reachable-cluster size, while the overall level identifies the importation
// scale. Both are fitted jointly by the same SBI as seam 2 (censored likelihood +
// inference.ComputePosterior), reusing the iterated-importance-sampling machinery.
package measles

import (
	"fmt"
	"math"
	"math/rand/v2"

	"github.com/umbralcalc/stochadex/pkg/inference"
)

// meanChainSizeGrid precomputes E[final chain size] over a grid of (R_local,
// reachableCluster), for fast lookup during the SBI. The susceptible pool for a
// chain is taken as the reachable cluster itself (s*N >> cluster for essentially
// all UTLAs), so the grid needs only these two axes.
type meanChainSizeGrid struct {
	rLocals, clusters []float64
	values            [][]float64 // values[ri][ci]
}

func linspace(lo, hi float64, n int) []float64 {
	out := make([]float64, n)
	for i := range out {
		out[i] = lo + (hi-lo)*float64(i)/float64(n-1)
	}
	return out
}

func buildMeanChainSizeGrid(
	rLocalLo, rLocalHi, clusterLo, clusterHi float64,
	nR, nC, nSims int, dispersion float64, seed uint64,
) *meanChainSizeGrid {
	g := &meanChainSizeGrid{
		rLocals:  linspace(rLocalLo, rLocalHi, nR),
		clusters: linspace(clusterLo, clusterHi, nC),
		values:   make([][]float64, nR),
	}
	rng := rand.New(rand.NewPCG(seed, seed+1))
	for ri, rl := range g.rLocals {
		g.values[ri] = make([]float64, nC)
		for ci, cl := range g.clusters {
			total := 0.0
			for s := 0; s < nSims; s++ {
				total += float64(simulateCluster(rl, dispersion, cl, 1<<30, rng))
			}
			g.values[ri][ci] = total / float64(nSims)
		}
	}
	return g
}

// lookup returns a bilinearly-interpolated mean chain size, clamped to the grid.
func (g *meanChainSizeGrid) lookup(rLocal, cluster float64) float64 {
	ri, rf := gridPos(g.rLocals, rLocal)
	ci, cf := gridPos(g.clusters, cluster)
	v00 := g.values[ri][ci]
	v10 := g.values[min(ri+1, len(g.rLocals)-1)][ci]
	v01 := g.values[ri][min(ci+1, len(g.clusters)-1)]
	v11 := g.values[min(ri+1, len(g.rLocals)-1)][min(ci+1, len(g.clusters)-1)]
	return (1-rf)*(1-cf)*v00 + rf*(1-cf)*v10 + (1-rf)*cf*v01 + rf*cf*v11
}

// gridPos returns the lower index and the fractional offset for x on axis.
func gridPos(axis []float64, x float64) (int, float64) {
	if x <= axis[0] {
		return 0, 0
	}
	if x >= axis[len(axis)-1] {
		return len(axis) - 1, 0
	}
	for i := 0; i < len(axis)-1; i++ {
		if x < axis[i+1] {
			return i, (x - axis[i]) / (axis[i+1] - axis[i])
		}
	}
	return len(axis) - 1, 0
}

// ReachableClusterPosterior summarises the calibration.
type ReachableClusterPosterior struct {
	ImportScale         float64
	ReachableCluster    float64
	ReachableClusterStd float64
	ReachableClusterCI95 [2]float64
	R0Used              float64
	EffectiveSampleSize float64
	NumParticles        int
}

// CalibrateReachableCluster fits (importScale, reachableCluster) to the observed
// censored UTLA totals. r0 is the (fixed) basic reproduction number used to map
// susceptibility to R_local; the reachable-cluster prior bounds the grid.
func CalibrateReachableCluster(
	susceptibility, receptivity, observedCounts []float64,
	threshold int, r0, dispersion float64,
	scalePrior, clusterPrior inference.UniformPrior,
	numParticles int, seed uint64,
) (*ReachableClusterPosterior, error) {
	n := len(susceptibility)
	if len(receptivity) != n || len(observedCounts) != n {
		return nil, fmt.Errorf("inputs must align (%d)", n)
	}

	rLocal := make([]float64, n)
	rMin, rMax := math.Inf(1), math.Inf(-1)
	for i, s := range susceptibility {
		rLocal[i] = r0 * s
		rMin, rMax = math.Min(rMin, rLocal[i]), math.Max(rMax, rLocal[i])
	}
	grid := buildMeanChainSizeGrid(
		math.Max(0, rMin*0.9), rMax*1.1,
		math.Max(1, clusterPrior.Lo), clusterPrior.Hi,
		24, 24, 3000, dispersion, seed)

	loglik := func(importScale, cluster float64) float64 {
		mu := make([]float64, n)
		for i := range mu {
			mu[i] = importScale * receptivity[i] * grid.lookup(rLocal[i], cluster)
		}
		return CensoredPoissonLogLike(observedCounts, mu, threshold)
	}

	res, clusterVals, ess := fitTwoParamSBI(
		scalePrior, clusterPrior, loglik, numParticles, seed+13)
	ci := inference.WeightedQuantiles(clusterVals, res.Weights, []float64{0.025, 0.975})

	return &ReachableClusterPosterior{
		ImportScale:          res.PosteriorMean[0],
		ReachableCluster:     res.PosteriorMean[1],
		ReachableClusterStd:  res.PosteriorStd[1],
		ReachableClusterCI95: [2]float64{ci[0], ci[1]},
		R0Used:               r0,
		EffectiveSampleSize:  ess,
		NumParticles:         numParticles,
	}, nil
}

// fitTwoParamSBI runs iterated importance sampling (population Monte Carlo) over
// two uniform-prior parameters, returning the posterior, the second parameter's
// per-particle values (for weighted quantiles), and the effective sample size.
// Shared machinery with the seam-2 case-link fit.
func fitTwoParamSBI(
	priorA, priorB inference.UniformPrior,
	loglik func(a, b float64) float64,
	numParticles int, seed uint64,
) (*inference.SMCResult, []float64, float64) {
	rng := rand.New(rand.NewPCG(seed, seed))
	logPriorConst := priorA.LogPDF(0.5*(priorA.Lo+priorA.Hi)) +
		priorB.LogPDF(0.5*(priorB.Lo+priorB.Hi))

	var res *inference.SMCResult
	var particles [][]float64
	var proposal *mvn2D
	const rounds = 4

	for round := 0; round < rounds; round++ {
		particles = make([][]float64, numParticles)
		logW := make([]float64, numParticles)
		for p := 0; p < numParticles; p++ {
			var a, b, logProp float64
			if proposal == nil {
				a, b = priorA.Sample(rng), priorB.Sample(rng)
			} else {
				for {
					a, b = proposal.sample(rng)
					if priorA.InSupport(a) && priorB.InSupport(b) {
						break
					}
				}
				logProp = proposal.logPDF(a, b)
			}
			particles[p] = []float64{a, b}
			if proposal == nil {
				logW[p] = loglik(a, b)
			} else {
				logW[p] = loglik(a, b) + logPriorConst - logProp
			}
		}
		res = inference.ComputePosterior([]string{"scale", "cluster"}, particles, logW, nil)
		next, ok := newMVN2D(res.PosteriorMean, res.PosteriorCov, 2.5)
		if !ok {
			break
		}
		proposal = next
	}

	sumSq := 0.0
	bVals := make([]float64, len(particles))
	for i, w := range res.Weights {
		sumSq += w * w
		bVals[i] = particles[i][1]
	}
	ess := 0.0
	if sumSq > 0 {
		ess = 1 / sumSq
	}
	return res, bVals, ess
}
