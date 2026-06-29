// This file is the susceptibility layer's ingest + sub-model A core (gating
// check #5): load COVER MMR coverage (dat/cover_mmr_coverage.csv), turn it into
// an effective-susceptibility surface per UTLA, and reconcile it onto the
// boundary/adjacency set (#4) — handling the two documented join mismatches:
//
//  1. Code-vintage: COVER still emits the old county codes for areas that became
//     unitary in 2023 (North Yorkshire, Somerset). coverCodeCrosswalk remaps them.
//  2. Small-LA combination (a documented COVER disclosure caveat): some UTLAs
//     (Rutland, City of London, Isles of Scilly) have no COVER row at all. Those
//     get no direct observation and are imputed by the CAR spatial prior from
//     their neighbours — except islands with no land neighbour, which are
//     assigned their documented combining-parent's value (Isles of Scilly ->
//     Cornwall) since the prior cannot pool an isolated node.
package measles

import (
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"strconv"
)

// Default MMR vaccine efficacies (per-dose, measles component). One dose is
// ~93% effective, two doses ~97% — the standard public-health figures.
const (
	DefaultMMR1Efficacy = 0.93
	DefaultMMR2Efficacy = 0.97
)

// coverCodeCrosswalk maps COVER's (old) geography code to the current boundary
// code for areas reorganised after COVER's coding was set. Same physical area.
var coverCodeCrosswalk = map[string]string{
	"E10000023": "E06000065", // North Yorkshire (county -> unitary, 2023)
	"E10000027": "E06000066", // Somerset (county -> unitary, 2023)
}

// combinedIntoNeighbour records UTLAs that COVER suppresses by combination and
// which have no usable spatial neighbour to borrow from (isolated islands), so
// must take a documented parent's value instead of CAR imputation.
var combinedIntoNeighbour = map[string]string{
	"E06000053": "E06000052", // Isles of Scilly -> Cornwall (island, combined)
}

// CoverageRecord is one row of the long COVER coverage table.
type CoverageRecord struct {
	Code        string // ONS geography code (COVER vintage)
	Name        string
	Dose        string // "MMR1" | "MMR2"
	Stratum     string // "24m" | "5y"
	Year        int
	Date        string
	CoveragePct float64 // 0..100
}

// LoadCoverage reads dat/cover_mmr_coverage.csv.
func LoadCoverage(path string) ([]CoverageRecord, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	rows, err := csv.NewReader(f).ReadAll()
	if err != nil {
		return nil, err
	}
	if len(rows) < 2 {
		return nil, fmt.Errorf("no coverage data in %s", path)
	}
	col := map[string]int{}
	for i, h := range rows[0] {
		col[h] = i
	}
	out := make([]CoverageRecord, 0, len(rows)-1)
	for _, r := range rows[1:] {
		year, err := strconv.Atoi(r[col["year"]])
		if err != nil {
			return nil, fmt.Errorf("bad year %q: %w", r[col["year"]], err)
		}
		pct, err := strconv.ParseFloat(r[col["coverage_pct"]], 64)
		if err != nil {
			return nil, fmt.Errorf("bad coverage_pct %q: %w", r[col["coverage_pct"]], err)
		}
		out = append(out, CoverageRecord{
			Code:        r[col["geography_code"]],
			Name:        r[col["geography"]],
			Dose:        r[col["dose"]],
			Stratum:     r[col["stratum"]],
			Year:        year,
			Date:        r[col["date"]],
			CoveragePct: pct,
		})
	}
	return out, nil
}

// SusceptibilityFromCoverage converts one- and two-dose MMR coverage (as
// fractions in [0,1]) into the effective susceptible fraction, given per-dose
// efficacies. Children with two doses are protected with probability e2; those
// with only one dose with probability e1; the rest are susceptible:
//
//	immune = c2*e2 + (c1 - c2)*e1   (c1 >= c2 expected; clamped if not)
//	s      = 1 - immune
//
// This is a cohort-free snapshot susceptibility (no waning / natural immunity),
// deliberately simple for the surface; richer cohort accounting is sub-model A's
// later refinement.
func SusceptibilityFromCoverage(c1, c2, e1, e2 float64) float64 {
	if c2 > c1 {
		c2 = c1 // two-dose coverage cannot exceed one-dose; guard data noise
	}
	immune := c2*e2 + (c1-c2)*e1
	s := 1 - immune
	if s < 0 {
		s = 0
	}
	return s
}

// latestDoseCoverage returns, per geography code, the latest-year MMR1 and MMR2
// coverage fraction for the given stratum. Codes are crosswalked to current
// boundary codes.
func latestDoseCoverage(records []CoverageRecord, stratum string) map[string]struct{ c1, c2 float64 } {
	type yv struct {
		year int
		pct  float64
	}
	latest := map[string]map[string]yv{} // code -> dose -> latest year/value
	for _, r := range records {
		if r.Stratum != stratum {
			continue
		}
		code := r.Code
		if mapped, ok := coverCodeCrosswalk[code]; ok {
			code = mapped
		}
		if latest[code] == nil {
			latest[code] = map[string]yv{}
		}
		if cur, ok := latest[code][r.Dose]; !ok || r.Year > cur.year {
			latest[code][r.Dose] = yv{r.Year, r.CoveragePct}
		}
	}
	out := map[string]struct{ c1, c2 float64 }{}
	for code, doses := range latest {
		out[code] = struct{ c1, c2 float64 }{
			c1: doses["MMR1"].pct / 100,
			c2: doses["MMR2"].pct / 100,
		}
	}
	return out
}

// SusceptibilitySurface is the reconciled, spatially-smoothed susceptibility
// field over a UTLA adjacency graph.
type SusceptibilitySurface struct {
	Graph    *AdjacencyGraph
	Raw      []float64 // per-node observed susceptibility (NaN where unobserved)
	Smoothed []float64 // CAR posterior-mean surface (all nodes finite)
	Imputed  []bool    // node had no direct COVER observation (neighbour-imputed)
}

// BuildSusceptibilitySurface joins COVER coverage onto the graph, computes
// per-UTLA susceptibility for the given stratum, and CAR-smooths it. Nodes with
// no COVER row are left unobserved (obsPrec 0) so the spatial prior imputes them
// from neighbours; documented island combinations are assigned their parent's
// value first so no connected component is left unanchored.
func BuildSusceptibilitySurface(
	g *AdjacencyGraph,
	records []CoverageRecord,
	stratum string,
	e1, e2 float64,
	tau, obsPrecision float64,
) (*SusceptibilitySurface, error) {
	cov := latestDoseCoverage(records, stratum)
	n := g.NumNodes()
	raw := make([]float64, n)
	obsPrec := make([]float64, n)
	imputed := make([]bool, n)

	for i, code := range g.Codes {
		src := code
		if parent, ok := combinedIntoNeighbour[code]; ok {
			src = parent // island combination: take the documented parent's data
		}
		if c, ok := cov[src]; ok && c.c1 > 0 {
			raw[i] = SusceptibilityFromCoverage(c.c1, c.c2, e1, e2)
			obsPrec[i] = obsPrecision
			imputed[i] = combinedIntoNeighbour[code] != "" // combined counts as imputed
		} else {
			raw[i] = math.NaN()
			obsPrec[i] = 0
			imputed[i] = true
		}
	}

	y := make([]float64, n)
	for i := range raw {
		if obsPrec[i] > 0 {
			y[i] = raw[i]
		}
	}
	s := &ICARSmoother{Graph: g, Tau: tau}
	smoothed, err := s.PosteriorMean(y, obsPrec)
	if err != nil {
		return nil, fmt.Errorf("smoothing susceptibility surface: %w", err)
	}
	return &SusceptibilitySurface{
		Graph: g, Raw: raw, Smoothed: smoothed, Imputed: imputed,
	}, nil
}
