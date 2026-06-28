// This file implements the spatial layer (gating check #4): an intrinsic
// conditional autoregressive (ICAR / Besag) prior over the UTLA adjacency graph,
// and the Gaussian smoother that pools strength across neighbours.
//
// Why it matters (PLAN.md sub-model A): UTLA coverage and case data are sparse
// and — for cases — heavily suppressed (see censored_poisson.go). A spatial prior
// that borrows strength from neighbouring UTLAs is what turns that sparse,
// censored signal into a usable susceptibility surface. A fully-suppressed UTLA
// carries no direct information (its observation precision is ~0); under this
// prior it is imputed from its neighbours rather than left at zero. That is the
// direct link between gating check #3 (censoring) and #4 (spatial pooling).
package measles

import (
	"encoding/csv"
	"fmt"
	"os"

	"gonum.org/v1/gonum/mat"
)

// AdjacencyGraph is an undirected graph of areal units (UTLAs), with a stable
// node ordering so it can index vectors of per-node quantities.
type AdjacencyGraph struct {
	Codes      []string       // node i's UTLA code, in matrix order
	Names      []string       // node i's UTLA name (parallel to Codes)
	Neighbours [][]int        // Neighbours[i] = indices adjacent to node i
	index      map[string]int // code -> node index
}

// NumNodes returns the number of areal units.
func (g *AdjacencyGraph) NumNodes() int { return len(g.Codes) }

// NumEdges returns the number of undirected edges.
func (g *AdjacencyGraph) NumEdges() int {
	total := 0
	for _, nbrs := range g.Neighbours {
		total += len(nbrs)
	}
	return total / 2
}

// Degree returns the number of neighbours of node i.
func (g *AdjacencyGraph) Degree(i int) int { return len(g.Neighbours[i]) }

// IndexOf returns the node index for a UTLA code (-1 if absent).
func (g *AdjacencyGraph) IndexOf(code string) int {
	if i, ok := g.index[code]; ok {
		return i
	}
	return -1
}

// LoadAdjacency reads the node list (code,name) and undirected edge list
// (code_a,code_b) produced by dat/build_adjacency.py and assembles the graph.
func LoadAdjacency(nodesCSV, edgesCSV string) (*AdjacencyGraph, error) {
	nodeRows, err := readCSV(nodesCSV)
	if err != nil {
		return nil, fmt.Errorf("reading nodes: %w", err)
	}
	g := &AdjacencyGraph{index: make(map[string]int)}
	for _, row := range nodeRows {
		if len(row) < 2 {
			return nil, fmt.Errorf("node row needs code,name: %v", row)
		}
		g.index[row[0]] = len(g.Codes)
		g.Codes = append(g.Codes, row[0])
		g.Names = append(g.Names, row[1])
	}
	g.Neighbours = make([][]int, len(g.Codes))

	edgeRows, err := readCSV(edgesCSV)
	if err != nil {
		return nil, fmt.Errorf("reading edges: %w", err)
	}
	seen := make(map[[2]int]bool)
	for _, row := range edgeRows {
		if len(row) < 2 {
			return nil, fmt.Errorf("edge row needs code_a,code_b: %v", row)
		}
		a, okA := g.index[row[0]]
		b, okB := g.index[row[1]]
		if !okA || !okB {
			return nil, fmt.Errorf("edge references unknown code: %v", row)
		}
		key := [2]int{a, b}
		if a > b {
			key = [2]int{b, a}
		}
		if a == b || seen[key] {
			continue
		}
		seen[key] = true
		g.Neighbours[a] = append(g.Neighbours[a], b)
		g.Neighbours[b] = append(g.Neighbours[b], a)
	}
	return g, nil
}

// readCSV reads a CSV with a header row, returning the data rows.
func readCSV(path string) ([][]string, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	rows, err := csv.NewReader(f).ReadAll()
	if err != nil {
		return nil, err
	}
	if len(rows) == 0 {
		return nil, fmt.Errorf("empty file: %s", path)
	}
	return rows[1:], nil // drop header
}

// ICARSmoother computes the Gaussian posterior mean of a latent field over an
// adjacency graph under an intrinsic CAR (Besag) prior, given per-node
// observations and their precisions.
//
// Prior:       x ~ ICAR(tau), precision tau*(D - W)  (D=degree diag, W=adjacency)
// Likelihood:  y_i ~ N(x_i, 1/obsPrec_i);  obsPrec_i = 0 means "no information"
//              (e.g. a fully-suppressed UTLA).
// Posterior mean solves the linear system
//
//	( tau*(D - W) + diag(obsPrec) ) x = diag(obsPrec) y .
//
// The ICAR precision is rank-deficient (it has a constant null space per
// connected component), so the obsPrec terms are what make the system solvable:
// each connected component must contain at least one node with obsPrec > 0,
// otherwise that component's level is unidentified and Factorize fails.
type ICARSmoother struct {
	Graph *AdjacencyGraph
	Tau   float64 // smoothing strength (larger = smoother / more pooling)
}

// PosteriorMean returns the smoothed field. y and obsPrec must have length
// NumNodes; obsPrec entries must be >= 0. It errors if any connected component
// is unanchored (all obsPrec == 0), which would leave that component's level
// undetermined.
func (s *ICARSmoother) PosteriorMean(y, obsPrec []float64) ([]float64, error) {
	n := s.Graph.NumNodes()
	if len(y) != n || len(obsPrec) != n {
		return nil, fmt.Errorf("y and obsPrec must have length %d", n)
	}
	for i, p := range obsPrec {
		if p < 0 {
			return nil, fmt.Errorf("obsPrec[%d] = %v < 0", i, p)
		}
	}

	// M = tau*(D - W) + diag(obsPrec), symmetric positive definite when anchored.
	m := mat.NewSymDense(n, nil)
	for i := 0; i < n; i++ {
		m.SetSym(i, i, s.Tau*float64(s.Graph.Degree(i))+obsPrec[i])
		for _, j := range s.Graph.Neighbours[i] {
			if j > i { // set each off-diagonal once
				m.SetSym(i, j, -s.Tau)
			}
		}
	}

	var chol mat.Cholesky
	if ok := chol.Factorize(m); !ok {
		return nil, fmt.Errorf(
			"posterior precision not positive definite: a connected component " +
				"has no anchoring observation (all obsPrec == 0)",
		)
	}

	b := mat.NewVecDense(n, nil)
	for i := 0; i < n; i++ {
		b.SetVec(i, obsPrec[i]*y[i])
	}
	var x mat.VecDense
	if err := chol.SolveVecTo(&x, b); err != nil {
		return nil, fmt.Errorf("solving posterior mean: %w", err)
	}
	return x.RawVector().Data, nil
}
