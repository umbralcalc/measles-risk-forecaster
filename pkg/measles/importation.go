// This file builds the coarse national importation-pressure index (PLAN.md
// importation-feasibility verdict, 2026-07-03). England's measles importations
// track upstream European activity, so the index is European measles incidence
// relative to baseline — deliberately annual and wide-uncertainty, because the
// data supports nothing finer: the Europe->UK link is weakly identified (the UK
// import-count calibration is prose-only, quarterly, tiny). It is an importation-
// *pressure* signal, not a precise weekly forecast.
//
// The index feeds per-UTLA seeding: when European activity surges, every
// receptive, susceptible UTLA's expected outbreak pressure rises together — the
// co-occurrence mechanism. This is the analogue of the "shared national epidemic
// anomaly" latent factor named in the model design, done at the honest
// resolution rather than an over-precise weekly latent process.
package measles

import (
	"fmt"
	"math"
	"sort"
	"strconv"
)

// Importation-elasticity band: expected UK importations scale as
// (European incidence / baseline)^gamma, with gamma poorly identified. The wide
// [0.5, 1.5] band around 1 is the honest uncertainty, carried into the index.
const (
	importElasticityLow  = 0.5
	importElasticityMid  = 1.0
	importElasticityHigh = 1.5
)

// LoadEuropeanIncidence reads dat/european_measles_annual.csv (year -> cases).
func LoadEuropeanIncidence(path string) (map[int]float64, error) {
	rows, col, err := readCSVWithHeader(path)
	if err != nil {
		return nil, err
	}
	out := make(map[int]float64, len(rows))
	for _, r := range rows {
		y, err := strconv.Atoi(r[col["year"]])
		if err != nil {
			return nil, fmt.Errorf("bad year %q: %w", r[col["year"]], err)
		}
		c, err := strconv.ParseFloat(r[col["cases"]], 64)
		if err != nil {
			return nil, fmt.Errorf("bad cases %q: %w", r[col["cases"]], err)
		}
		out[y] = c
	}
	return out, nil
}

// ImportationPressure is the coarse national index for one year: European
// activity relative to a baseline geomean, with a wide band from the elasticity
// uncertainty. Relative* are multiples of "typical" (1.0 = baseline).
type ImportationPressure struct {
	Year            int
	EuropeanCases   float64
	BaselineGeomean float64
	RelativeCentral float64 // (E/Ebar)^1.0
	RelativeLow     float64 // (E/Ebar)^0.5  — wide-band low
	RelativeHigh    float64 // (E/Ebar)^1.5  — wide-band high
}

// RecentBaselineYears returns up to n years before targetYear present in the
// series, excluding the 2020-2022 pandemic-collapse years (measles surveillance
// and transmission were anomalously suppressed, so they would bias "typical" low).
func RecentBaselineYears(series map[int]float64, targetYear, n int) []int {
	years := make([]int, 0, len(series))
	for y := range series {
		if y < targetYear && (y < 2020 || y > 2022) {
			years = append(years, y)
		}
	}
	sort.Sort(sort.Reverse(sort.IntSlice(years)))
	if len(years) > n {
		years = years[:n]
	}
	return years
}

// geomean returns the geometric mean of the series over the given years
// (positive values only).
func geomean(series map[int]float64, years []int) float64 {
	sumLog, k := 0.0, 0
	for _, y := range years {
		if v, ok := series[y]; ok && v > 0 {
			sumLog += math.Log(v)
			k++
		}
	}
	if k == 0 {
		return math.NaN()
	}
	return math.Exp(sumLog / float64(k))
}

// ComputeImportationPressure builds the index for `year` against a baseline
// geomean over baselineYears. If `year` is absent from the series, the most
// recent available year <= year is used (annual data lags).
func ComputeImportationPressure(
	series map[int]float64, year int, baselineYears []int,
) (ImportationPressure, error) {
	e, ok := series[year]
	usedYear := year
	if !ok {
		// Fall back to the most recent available year (annual data lag).
		best := -1
		for y := range series {
			if y <= year && y > best {
				best = y
			}
		}
		if best < 0 {
			return ImportationPressure{}, fmt.Errorf("no European data at or before %d", year)
		}
		usedYear, e = best, series[best]
	}
	// Never let the year being scored sit inside its own baseline.
	base := make([]int, 0, len(baselineYears))
	for _, y := range baselineYears {
		if y != usedYear {
			base = append(base, y)
		}
	}
	ebar := geomean(series, base)
	if math.IsNaN(ebar) || ebar <= 0 {
		return ImportationPressure{}, fmt.Errorf("bad baseline geomean")
	}
	ratio := e / ebar
	return ImportationPressure{
		Year:            usedYear,
		EuropeanCases:   e,
		BaselineGeomean: ebar,
		RelativeCentral: math.Pow(ratio, importElasticityMid),
		RelativeLow:     math.Pow(ratio, importElasticityLow),
		RelativeHigh:    math.Pow(ratio, importElasticityHigh),
	}, nil
}

// ReceptivityWeights returns per-UTLA importation receptivity as population
// share (a coarse proxy: larger places receive proportionally more importations).
// Weights sum to 1. Flight-connectivity refinement is a documented future step.
func ReceptivityWeights(population []float64) []float64 {
	total := 0.0
	for _, p := range population {
		total += p
	}
	w := make([]float64, len(population))
	if total <= 0 {
		return w
	}
	for i, p := range population {
		w[i] = p / total
	}
	return w
}

// ImportationAdjustedPressure combines the national index, per-UTLA receptivity,
// and local amplification P(self-sustaining) into a relative expected-outbreak-
// pressure score per UTLA:
//
//	pressure_i = index * receptivity_i * pSelfSustaining_i.
//
// This is the co-occurrence mechanism: a higher national index raises every
// UTLA's pressure together, concentrated where receptivity (population) and
// amplification (susceptibility) are both high.
func ImportationAdjustedPressure(index float64, receptivity, pSelfSustaining []float64) []float64 {
	out := make([]float64, len(receptivity))
	for i := range receptivity {
		out[i] = index * receptivity[i] * pSelfSustaining[i]
	}
	return out
}
