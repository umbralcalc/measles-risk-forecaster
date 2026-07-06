// This file scores the risk map out-of-sample: does the model's risk ordering
// put subsequently-affected UTLAs high, and does it beat the naive baselines the
// public already has — coverage-alone ("low coverage = high risk") and recent
// case counts ("where it's already happening")? Beating coverage-alone is the
// bar the plan sets. The scoring is rank-based (AUC / Mann-Whitney): no
// probability calibration required, so no extra fitting can leak in.
//
// The susceptibility surface is built from COVER coverage only (never from case
// data), so predicting year-t clusters from coverage as-of year t-1 has no
// case leakage by construction.
package measles

import "sort"

// ClusterLabels returns, aligned to the graph, whether each UTLA had a listed
// cluster (>=10 confirmed cases) in the report year — the scoring target.
func ClusterLabels(g *AdjacencyGraph, cases []UTLACumulativeCase, year int) []bool {
	obs := BuildCensoredCaseObservations(g, cases, year)
	labels := make([]bool, g.NumNodes())
	for i := range labels {
		labels[i] = !obs.Censored[i] // listed (>=10 cases) = a cluster occurred
	}
	return labels
}

// AUC is the area under the ROC curve (Mann-Whitney statistic): the probability
// that a randomly chosen positive (cluster) is scored higher than a randomly
// chosen negative. 0.5 is no skill; 1.0 is perfect ranking. Ties get average
// ranks. Returns 0.5 if either class is empty.
func AUC(scores []float64, labels []bool) float64 {
	n := len(scores)
	type sl struct {
		s float64
		y bool
	}
	arr := make([]sl, n)
	for i := range arr {
		arr[i] = sl{scores[i], labels[i]}
	}
	sort.Slice(arr, func(a, b int) bool { return arr[a].s < arr[b].s })

	ranks := make([]float64, n)
	for i := 0; i < n; {
		j := i
		for j < n && arr[j].s == arr[i].s {
			j++
		}
		avg := float64(i+j+1) / 2.0 // average of 1-based ranks (i+1 .. j)
		for k := i; k < j; k++ {
			ranks[k] = avg
		}
		i = j
	}

	sumPos, nPos := 0.0, 0
	for k := 0; k < n; k++ {
		if arr[k].y {
			sumPos += ranks[k]
			nPos++
		}
	}
	nNeg := n - nPos
	if nPos == 0 || nNeg == 0 {
		return 0.5
	}
	return (sumPos - float64(nPos)*float64(nPos+1)/2) / (float64(nPos) * float64(nNeg))
}
