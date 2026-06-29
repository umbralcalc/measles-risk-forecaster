package measles

import (
	"testing"

	"github.com/umbralcalc/stochadex/pkg/simulator"
)

// branchingImpls builds the stochadex Implementations for a single branching
// partition run for the given number of generations, writing to store.
func branchingImpls(store *simulator.StateTimeStorage, generations int) *simulator.Implementations {
	return &simulator.Implementations{
		Iterations:      []simulator.Iteration{&BranchingProcessIteration{}},
		OutputCondition: &simulator.EveryStepOutputCondition{},
		OutputFunction:  &simulator.StateTimeStorageOutputFunction{Store: store},
		TerminationCondition: &simulator.NumberOfStepsTerminationCondition{
			MaxNumberOfSteps: generations,
		},
		TimestepFunction: &simulator.ConstantTimestepFunction{Stepsize: 1.0},
	}
}

// TestBranchingIterationRunsUnderStochadex is the integration check: the model
// runs under the stochadex partition coordinator and the test harness, driven by
// transmission_settings.yaml — i.e. it is a well-formed stochadex Iteration.
func TestBranchingIterationRunsUnderStochadex(t *testing.T) {
	t.Run("coordinator", func(t *testing.T) {
		settings := simulator.LoadSettingsFromYaml("./transmission_settings.yaml")
		store := simulator.NewStateTimeStorage()
		impls := branchingImpls(store, 30)
		for i, it := range impls.Iterations { // configure the instances the coordinator runs
			it.Configure(i, settings)
		}
		simulator.NewPartitionCoordinator(settings, impls).Run()
	})

	t.Run("harness", func(t *testing.T) {
		settings := simulator.LoadSettingsFromYaml("./transmission_settings.yaml")
		store := simulator.NewStateTimeStorage()
		if err := simulator.RunWithHarnesses(settings, branchingImpls(store, 30)); err != nil {
			t.Errorf("harness failed: %v", err)
		}
	})
}

// runBranchingOnce runs the Iteration via the coordinator at a given
// susceptibility and seed, returning the final cumulative cluster size.
func runBranchingOnce(t *testing.T, susceptibility float64, seed uint64, cap float64) float64 {
	t.Helper()
	settings := simulator.LoadSettingsFromYaml("./transmission_settings.yaml")
	settings.Iterations[0].Params.Map["susceptibility"] = []float64{susceptibility}
	settings.Iterations[0].Params.Map["outbreak_cap"] = []float64{cap}
	settings.Iterations[0].Seed = seed

	store := simulator.NewStateTimeStorage()
	impls := branchingImpls(store, 50)
	for i, it := range impls.Iterations {
		it.Configure(i, settings)
	}
	simulator.NewPartitionCoordinator(settings, impls).Run()

	traj := store.GetValues("branching") // [][infectious, cumulative] per step
	last := traj[len(traj)-1]
	cumulative := last[1]

	// Sanity: cumulative is monotone non-decreasing along the trajectory.
	prev := 1.0
	for _, row := range traj {
		if row[1] < prev-1e-9 {
			t.Fatalf("cumulative decreased: %v < %v", row[1], prev)
		}
		prev = row[1]
	}
	return cumulative
}

// TestBranchingIterationRegime confirms the engine-native model reproduces the
// regime distinction the kernel tests established: supercritical susceptibility
// (R_local >> 1) drives outbreaks to the cap a substantial fraction of the time,
// while subcritical susceptibility (below the herd-immunity threshold) does not.
func TestBranchingIterationRegime(t *testing.T) {
	const cap = 2000.0
	const reps = 200

	largeHi := 0
	for seed := uint64(0); seed < reps; seed++ {
		if runBranchingOnce(t, 0.30, seed, cap) >= cap {
			largeHi++
		}
	}
	largeLo := 0
	for seed := uint64(0); seed < reps; seed++ {
		if runBranchingOnce(t, 0.04, seed, cap) >= cap {
			largeLo++
		}
	}
	fracHi := float64(largeHi) / reps
	fracLo := float64(largeLo) / reps
	t.Logf("P(large): supercritical s=0.30 -> %.2f ; subcritical s=0.04 -> %.2f",
		fracHi, fracLo)

	if fracHi < 0.3 {
		t.Errorf("supercritical P(large) = %.2f, expected substantial", fracHi)
	}
	if fracLo > 0.02 {
		t.Errorf("subcritical P(large) = %.2f, expected ~0", fracLo)
	}
}
