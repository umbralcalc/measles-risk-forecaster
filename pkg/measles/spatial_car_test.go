package measles

import (
	"math"
	"math/rand/v2"
	"testing"
)

const (
	nodesCSV = "../../dat/utla_nodes.csv"
	edgesCSV = "../../dat/utla_adjacency.csv"
)

func loadRealGraph(t *testing.T) *AdjacencyGraph {
	t.Helper()
	g, err := LoadAdjacency(nodesCSV, edgesCSV)
	if err != nil {
		t.Fatalf("LoadAdjacency: %v (did you run dat/fetch_geography.sh + "+
			"dat/build_adjacency.py?)", err)
	}
	return g
}

func nodeByName(g *AdjacencyGraph, name string) int {
	for i, n := range g.Names {
		if n == name {
			return i
		}
	}
	return -1
}

// TestRealAdjacencyGraph validates the graph the Python pipeline produced, from
// Go: the expected size, and geographically-known neighbour relations. If the
// boundary source or adjacency logic regresses, these break.
func TestRealAdjacencyGraph(t *testing.T) {
	g := loadRealGraph(t)

	if g.NumNodes() != 153 {
		t.Errorf("nodes = %d, want 153", g.NumNodes())
	}
	if g.NumEdges() != 360 {
		t.Errorf("edges = %d, want 360", g.NumEdges())
	}

	// Westminster's four neighbours (a well-understood case).
	wmin := nodeByName(g, "Westminster")
	if wmin < 0 {
		t.Fatal("Westminster not found")
	}
	got := map[string]bool{}
	for _, j := range g.Neighbours[wmin] {
		got[g.Names[j]] = true
	}
	for _, want := range []string{
		"Camden", "City of London", "Kensington and Chelsea", "Brent",
	} {
		if !got[want] {
			t.Errorf("Westminster missing neighbour %q", want)
		}
	}
	if len(got) != 4 {
		t.Errorf("Westminster has %d neighbours, want 4", len(got))
	}

	// Cornwall borders only Devon on land.
	corn := nodeByName(g, "Cornwall")
	if g.Degree(corn) != 1 || g.Names[g.Neighbours[corn][0]] != "Devon" {
		t.Errorf("Cornwall neighbours = %v, want [Devon]", neighbourNames(g, corn))
	}

	// Islands are isolated (and must be handled, not silently unpooled).
	for _, island := range []string{"Isle of Wight", "Isles of Scilly"} {
		if i := nodeByName(g, island); i < 0 || g.Degree(i) != 0 {
			t.Errorf("%s degree = %d, want 0 (isolated)", island, g.Degree(i))
		}
	}
}

func neighbourNames(g *AdjacencyGraph, i int) []string {
	out := make([]string, 0, g.Degree(i))
	for _, j := range g.Neighbours[i] {
		out = append(out, g.Names[j])
	}
	return out
}

// smoothFieldOverGraph builds a field that varies smoothly with respect to the
// adjacency structure, by diffusing random node values (repeated neighbour
// averaging). This is exactly the kind of spatially-clustered surface the ICAR
// prior is designed for — coverage pockets that bleed across UTLA borders.
func smoothFieldOverGraph(g *AdjacencyGraph, rounds int, rng *rand.Rand) []float64 {
	n := g.NumNodes()
	f := make([]float64, n)
	for i := range f {
		f[i] = rng.NormFloat64()
	}
	for r := 0; r < rounds; r++ {
		next := make([]float64, n)
		for i := 0; i < n; i++ {
			sum, count := f[i], 1.0
			for _, j := range g.Neighbours[i] {
				sum += f[j]
				count++
			}
			next[i] = sum / count
		}
		f = next
	}
	return f
}

func rmse(a, b []float64, idx []int) float64 {
	sum := 0.0
	for _, i := range idx {
		d := a[i] - b[i]
		sum += d * d
	}
	return math.Sqrt(sum / float64(len(idx)))
}

// TestICARSmootherPoolsNoise checks the headline behaviour: on the real UTLA
// graph, smoothing a noisy observation of a spatially-smooth truth recovers the
// truth better than the raw observations do. Pooling helps.
func TestICARSmootherPoolsNoise(t *testing.T) {
	g := loadRealGraph(t)
	rng := rand.New(rand.NewPCG(1, 2))

	truth := smoothFieldOverGraph(g, 12, rng)

	const noiseSD = 0.5
	n := g.NumNodes()
	obs := make([]float64, n)
	obsPrec := make([]float64, n)
	for i := range obs {
		obs[i] = truth[i] + noiseSD*rng.NormFloat64()
		obsPrec[i] = 1.0 / (noiseSD * noiseSD)
	}

	s := &ICARSmoother{Graph: g, Tau: 4.0}
	post, err := s.PosteriorMean(obs, obsPrec)
	if err != nil {
		t.Fatalf("PosteriorMean: %v", err)
	}

	all := make([]int, n)
	for i := range all {
		all[i] = i
	}
	rawErr := rmse(obs, truth, all)
	smoothErr := rmse(post, truth, all)
	t.Logf("RMSE raw=%.4f  smoothed=%.4f  (%.0f%% reduction)",
		rawErr, smoothErr, 100*(1-smoothErr/rawErr))

	if smoothErr >= rawErr {
		t.Errorf("smoothing did not help: raw=%.4f smoothed=%.4f", rawErr, smoothErr)
	}
}

// TestICARSmootherImputesCensored is the #3<->#4 link: nodes with no
// observation (obsPrec == 0, the fully-suppressed UTLA case) must be imputed
// from their neighbours, not left at zero. We assert the *exact* ICAR invariant:
// at an unanchored node the posterior mean equals the average of its neighbours'
// posterior means.
func TestICARSmootherImputesCensored(t *testing.T) {
	g := loadRealGraph(t)
	rng := rand.New(rand.NewPCG(7, 11))
	truth := smoothFieldOverGraph(g, 12, rng)
	n := g.NumNodes()

	// Censor ~30% of non-island nodes (give them no observation).
	obs := make([]float64, n)
	obsPrec := make([]float64, n)
	censored := make([]bool, n)
	nCensored := 0
	for i := 0; i < n; i++ {
		if g.Degree(i) > 0 && rng.Float64() < 0.30 {
			censored[i] = true
			obs[i] = 0 // the naive value we must NOT end up at
			obsPrec[i] = 0
			nCensored++
		} else {
			obs[i] = truth[i] + 0.3*rng.NormFloat64()
			obsPrec[i] = 1.0 / 0.09
		}
	}

	s := &ICARSmoother{Graph: g, Tau: 3.0}
	post, err := s.PosteriorMean(obs, obsPrec)
	if err != nil {
		t.Fatalf("PosteriorMean: %v", err)
	}

	imputeErr, zeroNaiveErr := 0.0, 0.0
	cidx := make([]int, 0, nCensored)
	for i := 0; i < n; i++ {
		if !censored[i] {
			continue
		}
		cidx = append(cidx, i)

		// Exact invariant: posterior at an unanchored node == neighbour mean.
		mean := 0.0
		for _, j := range g.Neighbours[i] {
			mean += post[j]
		}
		mean /= float64(g.Degree(i))
		if math.Abs(post[i]-mean) > 1e-7 {
			t.Errorf("node %d: posterior %.6f != neighbour mean %.6f",
				i, post[i], mean)
		}
		imputeErr += math.Abs(post[i] - truth[i])
		zeroNaiveErr += math.Abs(0 - truth[i])
	}
	imputeErr /= float64(nCensored)
	zeroNaiveErr /= float64(nCensored)
	t.Logf("censored=%d  mean|imputed-truth|=%.4f  mean|zero-truth|=%.4f",
		nCensored, imputeErr, zeroNaiveErr)

	// Imputation from neighbours must beat the naive zero-fill at censored nodes.
	if imputeErr >= zeroNaiveErr {
		t.Errorf("neighbour imputation (%.4f) no better than zero-fill (%.4f)",
			imputeErr, zeroNaiveErr)
	}
}

// TestUnanchoredComponentErrors confirms the smoother refuses an ill-posed
// system rather than returning garbage: an isolated island with no observation
// is an unanchored connected component.
func TestUnanchoredComponentErrors(t *testing.T) {
	g := loadRealGraph(t)
	n := g.NumNodes()
	obs := make([]float64, n)
	obsPrec := make([]float64, n)
	for i := range obsPrec {
		obsPrec[i] = 1.0
	}
	// Remove the anchor from an isolated island -> unanchored component.
	iow := nodeByName(g, "Isle of Wight")
	obsPrec[iow] = 0

	s := &ICARSmoother{Graph: g, Tau: 1.0}
	if _, err := s.PosteriorMean(obs, obsPrec); err == nil {
		t.Error("expected error for unanchored island component, got nil")
	}
}
