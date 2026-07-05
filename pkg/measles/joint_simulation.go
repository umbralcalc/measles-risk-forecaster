// This file is the joint co-occurrence simulation: a genuinely multi-partition
// stochadex model in which a SHARED national importation latent feeds every
// UTLA's seed count, so outbreaks co-occur when a national importation surge
// hits. It is the one place multi-partition earns its keep (see the design
// discussion): the analytic importation multiplier gives correct per-UTLA
// marginals but discards the correlation structure; this simulation samples the
// JOINT distribution, so the national-total tail (many areas surging together)
// comes out correctly fat.
//
// Architecture (2 partitions, coupled by params_from_upstream):
//
//	partition 0  NationalImportationIteration  -> draws a shared national seed
//	             total M per scenario (from the wide importation-index band)
//	partition 1  JointOutbreakIteration        -> reads M, seeds each UTLA with
//	             Poisson(M * receptivity_i), then branches every UTLA one
//	             generation per step at its own R_local = R0 * s_i
//
// The shared M is what correlates the UTLAs: a high-M scenario seeds everyone at
// once. Partition 1 shares nextGeneration() with the standalone kernel and the
// single-UTLA BranchingProcessIteration, so all three stay identical.
package measles

import (
	"math"
	"math/rand/v2"

	"github.com/umbralcalc/stochadex/pkg/simulator"
	"gonum.org/v1/gonum/stat/distuv"
)

// NationalImportationIteration (partition 0) draws the shared national seed total
// M once per scenario and holds it across generations. M ~ log-uniform over the
// [seed_low, seed_high] band — the importation-pressure index's wide uncertainty.
//
//	State:  [M]
//	Params: seed_low, seed_high  (national seed-total band for the scenario)
type NationalImportationIteration struct {
	rng *rand.Rand
}

func (n *NationalImportationIteration) Configure(partitionIndex int, settings *simulator.Settings) {
	seed := settings.Iterations[partitionIndex].Seed
	n.rng = rand.New(rand.NewPCG(seed, seed))
}

func (n *NationalImportationIteration) Iterate(
	params *simulator.Params, partitionIndex int,
	stateHistories []*simulator.StateHistory,
	timestepsHistory *simulator.CumulativeTimestepsHistory,
) []float64 {
	if timestepsHistory.CurrentStepNumber > 1 {
		return []float64{stateHistories[partitionIndex].Values.At(0, 0)} // hold M
	}
	lo := params.Map["seed_low"][0]
	hi := params.Map["seed_high"][0]
	if lo <= 0 {
		lo = 1e-6
	}
	// Log-uniform draw across the band.
	m := math.Exp(math.Log(lo) + n.rng.Float64()*(math.Log(hi)-math.Log(lo)))
	return []float64{m}
}

// JointOutbreakIteration (partition 1) holds all UTLAs in one state vector and
// advances them jointly, seeded by the shared national total from partition 0.
//
//	State:  [infectious_1..N, cumulative_1..N]   (width 2N)
//	Params: susceptibility (N), receptivity (N), r0, dispersion, outbreak_cap,
//	        national_seed_total (from upstream partition 0)
type JointOutbreakIteration struct {
	rng *rand.Rand
}

func (j *JointOutbreakIteration) Configure(partitionIndex int, settings *simulator.Settings) {
	seed := settings.Iterations[partitionIndex].Seed
	j.rng = rand.New(rand.NewPCG(seed, seed))
}

func (j *JointOutbreakIteration) Iterate(
	params *simulator.Params, partitionIndex int,
	stateHistories []*simulator.StateHistory,
	timestepsHistory *simulator.CumulativeTimestepsHistory,
) []float64 {
	s := params.Map["susceptibility"]
	w := params.Map["receptivity"]
	pool := params.Map["susceptible_pool"]
	r0 := params.Map["r0"][0]
	dispersion := params.Map["dispersion"][0]
	m := params.Map["national_seed_total"][0]
	n := len(s)

	out := make([]float64, 2*n)
	pois := distuv.Poisson{Src: j.rng}
	if timestepsHistory.CurrentStepNumber == 1 {
		// Generation 0: seed each UTLA from the shared national total.
		for i := 0; i < n; i++ {
			pois.Lambda = m * w[i]
			seeds := pois.Rand()
			out[i] = seeds
			out[n+i] = seeds
		}
		return out
	}
	// Later generations: branch every UTLA with susceptible depletion.
	state := stateHistories[partitionIndex].Values
	for i := 0; i < n; i++ {
		infectious := int(state.At(0, i))
		cumulative := state.At(0, n+i)
		remaining := pool[i] - cumulative
		if infectious <= 0 || remaining <= 0 || pool[i] < 1 {
			out[i] = 0
			out[n+i] = cumulative
			continue
		}
		rEff := r0 * s[i] * (remaining / pool[i])
		next := nextGeneration(infectious, rEff, dispersion, j.rng)
		if float64(next) > remaining {
			next = int(remaining)
		}
		cumulative += float64(next)
		out[i] = float64(next)
		out[n+i] = cumulative
	}
	return out
}

// JointParams configures the joint ensemble.
type JointParams struct {
	R0, Dispersion    float64
	SeedLow, SeedHigh float64 // national seed-total band M ~ logUniform(low, high)
	MaxGenerations    int
}

// JointScenario is one sampled season: the drawn national seed total and the
// resulting per-UTLA cumulative case counts.
type JointScenario struct {
	NationalSeedTotal float64
	Cumulative        []float64
}

// RunJointEnsemble runs the 2-partition joint simulation nScenarios times and
// returns the sampled joint outcomes. susceptibility, receptivity and
// susceptiblePool align by UTLA (pool_i = s_i * N_i * fraction). If seedLow ==
// seedHigh the national total is fixed (no shared-latent correlation) — the
// independent baseline used to demonstrate co-occurrence.
func RunJointEnsemble(
	susceptibility, receptivity, susceptiblePool []float64,
	p JointParams, nScenarios int, seed uint64,
) []JointScenario {
	n := len(susceptibility)
	scenarios := make([]JointScenario, nScenarios)
	for k := 0; k < nScenarios; k++ {
		runSeed := seed + uint64(k)
		store := runOneJointScenario(susceptibility, receptivity, susceptiblePool, p, runSeed)
		joint := store.GetValues("joint")
		nat := store.GetValues("national")
		last := joint[len(joint)-1]
		cum := make([]float64, n)
		copy(cum, last[n:2*n])
		scenarios[k] = JointScenario{
			NationalSeedTotal: nat[len(nat)-1][0],
			Cumulative:        cum,
		}
	}
	return scenarios
}

func runOneJointScenario(
	susceptibility, receptivity, susceptiblePool []float64, p JointParams, runSeed uint64,
) *simulator.StateTimeStorage {
	n := len(susceptibility)
	jointInit := make([]float64, 2*n)

	settings := &simulator.Settings{
		Iterations: []simulator.IterationSettings{
			{
				Name: "national",
				Params: simulator.NewParams(map[string][]float64{
					"seed_low":  {p.SeedLow},
					"seed_high": {p.SeedHigh},
				}),
				InitStateValues:   []float64{0},
				StateWidth:        1,
				StateHistoryDepth: 2,
				Seed:              runSeed,
			},
			{
				Name: "joint",
				Params: simulator.NewParams(map[string][]float64{
					"susceptibility":      susceptibility,
					"receptivity":         receptivity,
					"susceptible_pool":    susceptiblePool,
					"r0":                  {p.R0},
					"dispersion":          {p.Dispersion},
					"national_seed_total": {0},
				}),
				ParamsFromUpstream: map[string]simulator.UpstreamConfig{
					"national_seed_total": {Upstream: 0},
				},
				InitStateValues:   jointInit,
				StateWidth:        2 * n,
				StateHistoryDepth: 2,
				Seed:              runSeed + 7919, // distinct stream
			},
		},
		InitTimeValue:         0,
		TimestepsHistoryDepth: 2,
	}
	settings.Init()

	store := simulator.NewStateTimeStorage()
	iterations := []simulator.Iteration{
		&NationalImportationIteration{},
		&JointOutbreakIteration{},
	}
	impl := &simulator.Implementations{
		Iterations:      iterations,
		OutputCondition: &simulator.EveryStepOutputCondition{},
		OutputFunction:  &simulator.StateTimeStorageOutputFunction{Store: store},
		TerminationCondition: &simulator.NumberOfStepsTerminationCondition{
			MaxNumberOfSteps: p.MaxGenerations,
		},
		TimestepFunction: &simulator.ConstantTimestepFunction{Stepsize: 1.0},
	}
	for i, it := range iterations {
		it.Configure(i, settings)
	}
	simulator.NewPartitionCoordinator(settings, impl).Run()
	return store
}
