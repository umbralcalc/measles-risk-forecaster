// This file is seam 1 of the stochadex integration: sub-model B's branching
// process expressed as a stochadex simulator.Iteration. Each Iterate advances
// the epidemic one generation, so a single simulation run is one realisation of
// a UTLA outbreak and an ensemble of seeds is the cluster-size distribution. This
// is the form that compiles to WASM for the dashboard, letting a reader run
// per-UTLA outbreak ensembles client-side (PLAN.md "Dashboard & distribution").
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
//	Params: r0            basic reproduction number
//	        susceptibility effective susceptible fraction s (R_local = r0 * s)
//	        dispersion    negative-binomial dispersion k (superspreading)
//	        outbreak_cap  cumulative size at which the outbreak is called "large"
//	                      and the process is frozen (absorbing)
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
	cap := params.Map["outbreak_cap"][0]

	// Absorbing states: extinct (no infectious left) or already a large outbreak.
	if infectious <= 0 || cumulative >= cap {
		return []float64{0, cumulative}
	}

	next := nextGeneration(infectious, r0*s, dispersion, b.rng)
	cumulative += float64(next)
	if cumulative > cap {
		cumulative = cap
	}
	return []float64{float64(next), cumulative}
}
