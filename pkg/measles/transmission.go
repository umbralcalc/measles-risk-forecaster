// This file is sub-model B (gating check #5): local transmission. It turns the
// susceptibility surface (sub-model A) into the headline risk quantities by
// simulating the local epidemic process many times and reading risk off the
// ensemble (PLAN.md "What we forecast"):
//
//   - R_local,i = R0 * s_i, with R0 an uncertain measles basic reproduction
//     number (~12-18). P(R_local > 1) = P(self-sustaining transmission given an
//     importation) is computed by integrating over R0.
//   - A stochastic branching process per UTLA, seeded by one importation, gives
//     the cluster-size predictive distribution and the large-outbreak tail
//     (offspring are negative-binomial to capture measles superspreading).
package measles

import (
	"math/rand/v2"
	"sort"

	"gonum.org/v1/gonum/stat/distuv"
)

// TransmissionParams configures sub-model B.
type TransmissionParams struct {
	R0Min, R0Max float64 // uncertain basic reproduction number range (e.g. 12,18)
	Dispersion   float64 // negative-binomial dispersion k (smaller = more superspreading)
	NumSims      int     // ensemble size per UTLA
	OutbreakCap  int     // cluster size at which we call it a "large outbreak" and stop
}

// DefaultTransmissionParams are the measles defaults used by the forecast.
func DefaultTransmissionParams() TransmissionParams {
	return TransmissionParams{
		R0Min: 12, R0Max: 18, Dispersion: 0.5, NumSims: 4000, OutbreakCap: 5000,
	}
}

// ProbSelfSustaining returns P(R_local > 1) = P(R0 * s > 1) with R0 ~ Uniform
// (R0Min, R0Max). This is the calibrated transmission-risk headline: the
// probability an importation triggers self-sustaining local transmission.
func ProbSelfSustaining(s, r0Min, r0Max float64) float64 {
	if s <= 0 {
		return 0
	}
	threshold := 1.0 / s // R0 must exceed this for R_local > 1
	if threshold <= r0Min {
		return 1
	}
	if threshold >= r0Max {
		return 0
	}
	return (r0Max - threshold) / (r0Max - r0Min)
}

// nextGeneration draws the total number of secondary cases produced by the
// current generation of `infectious` cases, each with negative-binomial offspring
// (mean rLocal, dispersion k). Because the sum of n iid NegBin(mean rLocal,
// dispersion k) offspring counts is itself Poisson(Gamma(shape = n*k, rate =
// k/rLocal)), the whole generation is one Gamma draw plus one Poisson draw —
// exact, and O(1) in the generation size rather than O(n). This single step is
// shared by the standalone kernel and the stochadex BranchingProcessIteration,
// so the two cannot drift apart.
func nextGeneration(infectious int, rLocal, dispersion float64, rng *rand.Rand) int {
	if infectious <= 0 || rLocal <= 0 {
		return 0
	}
	shape := float64(infectious) * dispersion
	rate := dispersion / rLocal
	lambda := distuv.Gamma{Alpha: shape, Beta: rate, Src: rng}.Rand()
	return int(distuv.Poisson{Lambda: lambda, Src: rng}.Rand())
}

// simulateCluster runs one Galton-Watson branching process seeded by a single
// importation. It returns the final cluster size, capped at cap (a supercritical
// process is stopped once it reaches cap and reported as cap).
func simulateCluster(rLocal, dispersion float64, cap int, rng *rand.Rand) int {
	size := 1
	infectious := 1
	for infectious > 0 && size < cap {
		infectious = nextGeneration(infectious, rLocal, dispersion, rng)
		size += infectious
		if size >= cap {
			return cap
		}
	}
	return size
}

// UTLARisk is the per-UTLA risk readout for the committed map.
type UTLARisk struct {
	Susceptibility    float64 // s_i from sub-model A
	ProbRLocalGt1     float64 // P(R_local > 1)
	MedianRLocal      float64 // R0Mid * s (median R_local)
	ClusterP50        int     // cluster-size quantiles (given one importation)
	ClusterP90        int
	ClusterP99        int
	ProbLargeOutbreak float64 // P(cluster reaches OutbreakCap)
	MeanClusterSize   float64
}

// ForecastUTLARisk runs the transmission ensemble for one UTLA at susceptibility
// s and reads the risk quantities off the simulated cluster sizes.
func ForecastUTLARisk(s float64, p TransmissionParams, rng *rand.Rand) UTLARisk {
	sizes := make([]int, p.NumSims)
	large, total := 0, 0.0
	for i := 0; i < p.NumSims; i++ {
		r0 := p.R0Min + rng.Float64()*(p.R0Max-p.R0Min)
		rLocal := r0 * s
		sz := simulateCluster(rLocal, p.Dispersion, p.OutbreakCap, rng)
		sizes[i] = sz
		total += float64(sz)
		if sz >= p.OutbreakCap {
			large++
		}
	}
	sort.Ints(sizes)
	q := func(f float64) int { return sizes[int(f*float64(len(sizes)-1))] }
	return UTLARisk{
		Susceptibility:    s,
		ProbRLocalGt1:     ProbSelfSustaining(s, p.R0Min, p.R0Max),
		MedianRLocal:      0.5 * (p.R0Min + p.R0Max) * s,
		ClusterP50:        q(0.50),
		ClusterP90:        q(0.90),
		ClusterP99:        q(0.99),
		ProbLargeOutbreak: float64(large) / float64(p.NumSims),
		MeanClusterSize:   total / float64(p.NumSims),
	}
}

// expectedTotalProgeny returns the analytic mean final size of a subcritical
// branching process seeded by one case: 1/(1-m) for m < 1. Used in tests.
func expectedTotalProgeny(m float64) float64 { return 1.0 / (1.0 - m) }
