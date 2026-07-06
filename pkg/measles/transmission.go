// This file is sub-model B (gating check #5): local transmission. It turns the
// susceptibility surface (sub-model A) into the headline risk quantities by
// simulating the local epidemic process many times and reading risk off the
// ensemble (README.md "What we forecast"):
//
//   - R_local,i = R0 * s_i, with R0 an uncertain measles basic reproduction
//     number (~12-18). P(R_local > 1) = P(self-sustaining transmission given an
//     importation) is computed by integrating over R0.
//   - A stochastic branching process per UTLA, seeded by one importation, gives
//     the cluster-size predictive distribution and the large-outbreak tail
//     (offspring are negative-binomial to capture measles superspreading).
//
// Susceptible depletion: the process is a stochastic SIR, not a pure branching
// process — the effective reproduction number declines as the local susceptible
// pool burns down, R_eff = R_local,0 * S/S_0. This gives a mechanistic final-size
// distribution (extinction, or a self-limiting outbreak that stops when
// susceptibles run out) rather than growth to an arbitrary cap. The susceptible
// pool is s * population * EffectivePoolFraction; the default fraction of 1.0 is
// the HOMOGENEOUS-MIXING UPPER BOUND on outbreak size — real outbreaks are
// smaller (heterogeneity, clustering, public-health response), so absolute large
// sizes are over-estimates and the fraction is the knob to reflect that.
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
	// LargeOutbreakThreshold: a cluster of at least this many cases counts as a
	// "large outbreak" for ProbLargeOutbreak.
	LargeOutbreakThreshold int
	// EffectivePoolFraction: fraction of the nominal susceptible pool (s*N) that is
	// well-mixed and reachable by one introduction. 1.0 = homogeneous-mixing upper
	// bound; lower values reflect heterogeneity/clustering/intervention.
	EffectivePoolFraction float64
	// ReachableClusterSize: absolute ceiling on the well-mixed susceptible
	// community reachable from a single introduction (a school/neighbourhood
	// scale). The susceptible pool is min(s*N*fraction, ReachableClusterSize). This
	// is what keeps absolute outbreak sizes finite and community-scaled rather than
	// UTLA-wide. It is a COARSE assumption pending calibration to observed outbreak
	// sizes — so absolute cluster-size magnitudes are illustrative; P(R_local>1)
	// and P(large outbreak) do not depend on it.
	ReachableClusterSize float64
	// MaxOutbreakSize: hard safety bound on a single simulated outbreak.
	MaxOutbreakSize int
}

// DefaultTransmissionParams are the measles defaults used by the forecast.
func DefaultTransmissionParams() TransmissionParams {
	return TransmissionParams{
		R0Min: 12, R0Max: 18, Dispersion: 0.5, NumSims: 4000,
		LargeOutbreakThreshold: 100, EffectivePoolFraction: 1.0,
		// ~400 from calibration to observed 2025 UTLA totals (see calibration.go);
		// the forecast overrides this with the freshly-fitted value.
		ReachableClusterSize: 400, MaxOutbreakSize: 5_000_000,
	}
}

// effectivePool returns the reachable susceptible pool for one introduction.
func effectivePool(s, population, fraction, reachableCeiling float64) float64 {
	pool := s * population * fraction
	if reachableCeiling > 0 && pool > reachableCeiling {
		return reachableCeiling
	}
	return pool
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

// simulateCluster runs one stochastic-SIR outbreak seeded by a single
// importation, with susceptible depletion: the effective reproduction number
// each generation is rLocal0 * remaining/susceptiblePool, so the outbreak
// self-limits when susceptibles run out. Returns the final size (including the
// index case), bounded by the pool and by maxSize (a safety cap only).
func simulateCluster(rLocal0, dispersion, susceptiblePool float64, maxSize int, rng *rand.Rand) int {
	if susceptiblePool < 1 || rLocal0 <= 0 {
		return 1 // no reachable susceptibles / non-transmissible: just the index case
	}
	cumulative := 1.0 // index case
	infectious := 1
	for infectious > 0 && cumulative < susceptiblePool && int(cumulative) < maxSize {
		remaining := susceptiblePool - cumulative
		rEff := rLocal0 * (remaining / susceptiblePool)
		next := nextGeneration(infectious, rEff, dispersion, rng)
		if float64(next) > remaining {
			next = int(remaining) // cannot infect more than remain susceptible
		}
		cumulative += float64(next)
		infectious = next
	}
	return int(cumulative)
}

// UTLARisk is the per-UTLA risk readout for the committed map.
type UTLARisk struct {
	Susceptibility    float64 // s_i from sub-model A
	ProbRLocalGt1     float64 // P(R_local > 1)
	MedianRLocal      float64 // R0Mid * s (median R_local)
	ClusterP50        int     // final outbreak-size quantiles (given one importation)
	ClusterP90        int
	ClusterP99        int
	ProbLargeOutbreak float64 // P(final size >= LargeOutbreakThreshold)
	MeanClusterSize   float64
}

// ForecastUTLARisk runs the transmission ensemble for one UTLA at susceptibility
// s and population, with susceptible depletion, and reads the risk quantities off
// the simulated final outbreak sizes.
func ForecastUTLARisk(s, population float64, p TransmissionParams, rng *rand.Rand) UTLARisk {
	pool := effectivePool(s, population, p.EffectivePoolFraction, p.ReachableClusterSize)
	sizes := make([]int, p.NumSims)
	large, total := 0, 0.0
	for i := 0; i < p.NumSims; i++ {
		r0 := p.R0Min + rng.Float64()*(p.R0Max-p.R0Min)
		sz := simulateCluster(r0*s, p.Dispersion, pool, p.MaxOutbreakSize, rng)
		sizes[i] = sz
		total += float64(sz)
		if sz >= p.LargeOutbreakThreshold {
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
