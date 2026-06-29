// This file ingests the case layers (gating check #1 resolved) into Go:
//   - Region x week counts (the dynamics layer; feeds the nowcast, sub-model C).
//   - UTLA cumulative counts from the report, where suppression is row omission:
//     a UTLA absent from the report had <10 cases, i.e. an interval-censored
//     count on {0,...,9}. BuildCensoredCaseObservations reconstructs that
//     censored set against the boundary UTLAs (#4) and emits the observation
//     vector the censored Poisson likelihood (#3) consumes — Suppressed (NaN)
//     for every UTLA the report omitted.
package measles

import (
	"encoding/csv"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// RegionWeekCase is one Region x week measles count.
type RegionWeekCase struct {
	RegionCode          string
	Region              string
	Year                int
	Epiweek             int
	Date                string
	Cases               float64
	InReportingDelay    bool
}

// LoadRegionWeekCases reads dat/measles_cases_region_week.csv.
func LoadRegionWeekCases(path string) ([]RegionWeekCase, error) {
	rows, col, err := readCSVWithHeader(path)
	if err != nil {
		return nil, err
	}
	out := make([]RegionWeekCase, 0, len(rows))
	for _, r := range rows {
		year, _ := strconv.Atoi(r[col["year"]])
		week, _ := strconv.Atoi(r[col["epiweek"]])
		cases, err := strconv.ParseFloat(r[col["cases"]], 64)
		if err != nil {
			return nil, fmt.Errorf("bad cases %q: %w", r[col["cases"]], err)
		}
		out = append(out, RegionWeekCase{
			RegionCode:       r[col["geography_code"]],
			Region:           r[col["region"]],
			Year:             year,
			Epiweek:          week,
			Date:             r[col["date"]],
			Cases:            cases,
			InReportingDelay: strings.EqualFold(r[col["in_reporting_delay_period"]], "true"),
		})
	}
	return out, nil
}

// UTLACumulativeCase is one listed UTLA total from a report year (UTLAs with
// <10 cases are NOT present — that absence is the censoring).
type UTLACumulativeCase struct {
	ReportYear int
	UTLA       string
	Region     string
	TotalCases int
}

// LoadUTLACases reads dat/measles_cases_utla.csv.
func LoadUTLACases(path string) ([]UTLACumulativeCase, error) {
	rows, col, err := readCSVWithHeader(path)
	if err != nil {
		return nil, err
	}
	out := make([]UTLACumulativeCase, 0, len(rows))
	for _, r := range rows {
		yr, _ := strconv.Atoi(r[col["report_year"]])
		total, err := strconv.Atoi(r[col["total_cases"]])
		if err != nil {
			return nil, fmt.Errorf("bad total_cases %q: %w", r[col["total_cases"]], err)
		}
		out = append(out, UTLACumulativeCase{
			ReportYear: yr,
			UTLA:       r[col["utla"]],
			Region:     r[col["region"]],
			TotalCases: total,
		})
	}
	return out, nil
}

// utlaNameAliases maps report UTLA names to boundary (CTYUA24NM) names where
// they differ in spelling. Extend as the reports surface new spellings.
var utlaNameAliases = map[string]string{
	"Bristol":            "Bristol, City of",
	"Herefordshire":      "Herefordshire, County of",
	"Kingston upon Hull":  "Kingston upon Hull, City of",
	"City of Kingston upon Hull": "Kingston upon Hull, City of",
}

func normaliseUTLA(name string) string {
	name = strings.Join(strings.Fields(name), " ")
	if alias, ok := utlaNameAliases[name]; ok {
		return alias
	}
	return name
}

// CensoredCaseObservations holds the per-node case observation vector for a
// report year, ready for the censored Poisson likelihood: listed UTLAs carry
// their observed count; every other (row-omitted) UTLA carries Suppressed (NaN).
type CensoredCaseObservations struct {
	Graph      *AdjacencyGraph
	ReportYear int
	Counts     []float64 // observed count, or Suppressed (NaN) if censored
	Censored   []bool    // true where the UTLA was row-omitted (<10 cases)
	Unmatched  []string  // report UTLA names that did not match any boundary node
}

// BuildCensoredCaseObservations reconstructs the censored case vector for one
// report year on the boundary/adjacency graph. The report lists only UTLAs with
// >=10 cases; every boundary UTLA not in that list is interval-censored [0,9].
func BuildCensoredCaseObservations(
	g *AdjacencyGraph, cases []UTLACumulativeCase, reportYear int,
) *CensoredCaseObservations {
	byName := make(map[string]int, g.NumNodes())
	for i, n := range g.Names {
		byName[n] = i
	}

	counts := make([]float64, g.NumNodes())
	censored := make([]bool, g.NumNodes())
	for i := range counts {
		counts[i] = Suppressed // default: censored unless the report lists it
		censored[i] = true
	}

	var unmatched []string
	for _, c := range cases {
		if c.ReportYear != reportYear {
			continue
		}
		idx, ok := byName[normaliseUTLA(c.UTLA)]
		if !ok {
			unmatched = append(unmatched, c.UTLA)
			continue
		}
		counts[idx] = float64(c.TotalCases)
		censored[idx] = false
	}
	return &CensoredCaseObservations{
		Graph: g, ReportYear: reportYear,
		Counts: counts, Censored: censored, Unmatched: unmatched,
	}
}

// readCSVWithHeader reads a CSV and returns the data rows plus a header->index map.
func readCSVWithHeader(path string) ([][]string, map[string]int, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	defer f.Close()
	rows, err := csv.NewReader(f).ReadAll()
	if err != nil {
		return nil, nil, err
	}
	if len(rows) < 2 {
		return nil, nil, fmt.Errorf("no data rows in %s", path)
	}
	col := make(map[string]int, len(rows[0]))
	for i, h := range rows[0] {
		col[h] = i
	}
	return rows[1:], col, nil
}
