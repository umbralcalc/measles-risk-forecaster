package measles

import (
	"fmt"
	"strconv"
)

// LoadPopulation reads dat/population_utla.csv (ONS mid-year estimates via Nomis)
// into a map from ONS geography code to population. It is the denominator/offset
// for the susceptibility->cases observation model (seam 2).
func LoadPopulation(path string) (map[string]float64, error) {
	rows, col, err := readCSVWithHeader(path)
	if err != nil {
		return nil, err
	}
	out := make(map[string]float64, len(rows))
	for _, r := range rows {
		pop, err := strconv.ParseFloat(r[col["population"]], 64)
		if err != nil {
			return nil, fmt.Errorf("bad population %q: %w", r[col["population"]], err)
		}
		out[r[col["code"]]] = pop
	}
	return out, nil
}

// PopulationVector aligns a population map to a graph's node order, returning the
// per-node population and an error if any node is missing (no silent zeros — a
// zero population would break the log offset).
func PopulationVector(g *AdjacencyGraph, pop map[string]float64) ([]float64, error) {
	out := make([]float64, g.NumNodes())
	for i, code := range g.Codes {
		p, ok := pop[code]
		if !ok || p <= 0 {
			return nil, fmt.Errorf("missing/zero population for %s (%s)", g.Names[i], code)
		}
		out[i] = p
	}
	return out, nil
}
