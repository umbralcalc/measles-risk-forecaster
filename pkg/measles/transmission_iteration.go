// This file expresses sub-model B's branching process as a stochadex
// simulator.Iteration. Each Iterate advances the epidemic one generation, so a
// single simulation run is one realisation of a UTLA outbreak and an ensemble of
// seeds is the cluster-size distribution — the engine-native form of the model,
// and the per-UTLA building block the joint co-occurrence simulation composes.
//
// It shares nextGeneration() with the standalone kernel (transmission.go), so the
// engine-native model and the unit-tested kernel are guaranteed identical.
package measles

import (
	"math/rand/v2"

	"github.com/umbralcalc/stochadex/pkg/simulator"
)

// BranchingProcessIteration advances a per-UTLA measles outbreak by one
// generation.
//
//	State:  [infectious, cumulative]  (current generation size, total cases so far)
//	Params: r0               basic reproduction number
//	        susceptibility   effective susceptible fraction s (R_local,0 = r0 * s)
//	        dispersion       negative-binomial dispersion k (superspreading)
//	        susceptible_pool absolute reachable susceptibles S_0 (= s * N * fraction);
//	                         R_eff = R_local,0 * remaining/S_0, so the outbreak
//	                         self-limits by susceptible depletion (absorbing when
//	                         susceptibles run out)
//
// Seed the simulation with init_state_values [1, 1] (one importation).
type BranchingProcessIteration struct {
	rng *rand.Rand
}

func (b *BranchingProcessIteration) Configure(
	partitionIndex int,
	settings *simulator.Settings,
) {
	seed := settings.Iterations[partitionIndex].Seed
	b.rng = rand.New(rand.NewPCG(seed, seed))
}

func (b *BranchingProcessIteration) Iterate(
	params *simulator.Params,
	partitionIndex int,
	stateHistories []*simulator.StateHistory,
	timestepsHistory *simulator.CumulativeTimestepsHistory,
) []float64 {
	state := stateHistories[partitionIndex]
	infectious := int(state.Values.At(0, 0))
	cumulative := state.Values.At(0, 1)

	r0 := params.Map["r0"][0]
	s := params.Map["susceptibility"][0]
	dispersion := params.Map["dispersion"][0]
	pool := params.Map["susceptible_pool"][0]

	// Absorbing states: extinct (no infectious left) or susceptibles exhausted.
	remaining := pool - cumulative
	if infectious <= 0 || remaining <= 0 {
		return []float64{0, cumulative}
	}

	// Effective reproduction number declines with susceptible depletion.
	rEff := r0 * s * (remaining / pool)
	next := nextGeneration(infectious, rEff, dispersion, b.rng)
	if float64(next) > remaining {
		next = int(remaining)
	}
	cumulative += float64(next)
	return []float64{float64(next), cumulative}
}
