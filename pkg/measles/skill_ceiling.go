// This file computes the skill-ceiling decomposition (cmd/skill-ceiling): how
// much of where-and-when measles actually happens is predictable from the
// susceptibility surface, versus irreducible. This is the project's load-bearing
// honest figure (README.md "the predictable-vs-irreducible decomposition").
//
// It has two axes, kept separate because they answer different questions and the
// feasibility study (see repo notes) showed one is easily mistaken for the other:
//
//   - Spatial: does the susceptibility surface predict WHERE cases concentrate,
//     beyond the mechanical population effect? The residual is the irreducible
//     spark-placement component the England & Wales literature established for
//     high-coverage settings.
//   - Temporal: how persistent/predictable is national activity week-to-week?
//     CRITICAL caveat carried in the output: this momentum is dominated by
//     transmission (an outbreak this week implies cases next week), NOT by
//     importation-rate predictability. It is an upper bound, not evidence that
//     imports are forecastable.
package measles

import "math"

// SpatialSkill is the spatial predictable-vs-irreducible decomposition, fitted to
// the censored UTLA case counts with population as a known exposure offset.
type SpatialSkill struct {
	LogLikNull     float64 // population-offset only (susceptibility unused, beta=0)
	LogLikFull     float64 // + susceptibility
	Alpha, Beta    float64 // full-model coefficients used
	LRStatistic    float64 // 2*(ll_full - ll_null), chi-square (1 df): is susceptibility informative?
	McFaddenR2     float64 // 1 - ll_full/ll_null: susceptibility skill beyond population
	IrreducibleFrac float64 // 1 - McFaddenR2: the share left to spark placement etc.
	NListed        int     // UTLAs with observed counts (>=10)
	NCensored      int     // UTLAs interval-censored [0,9]
}

// SpatialSkillDecomposition quantifies how much susceptibility explains about
// case concentration beyond population. Pass the full-model (alpha, beta) from
// FitCaseLink; the null model (beta=0) is fitted here. Inputs align by UTLA;
// censoredCounts carries observed counts or Suppressed for row-omitted UTLAs.
func SpatialSkillDecomposition(
	susceptibility, population, censoredCounts []float64,
	threshold int,
	alphaFull, betaFull float64,
) SpatialSkill {
	z := standardisedLogSusceptibility(susceptibility)
	offset := logPopOffsets(population)

	// Null model: fit alpha with susceptibility switched off (beta = 0).
	llNull := func(a float64) float64 {
		mu := make([]float64, len(offset))
		for i := range offset {
			mu[i] = math.Exp(a + offset[i])
		}
		return CensoredPoissonLogLike(censoredCounts, mu, threshold)
	}
	alphaNull := goldenSectionMax(llNull, -6, 12)
	logLikNull := llNull(alphaNull)

	// Full model log-likelihood at the supplied (alpha, beta).
	muFull := caseLinkMeans(alphaFull, betaFull, z, offset)
	logLikFull := CensoredPoissonLogLike(censoredCounts, muFull, threshold)

	listed, censored := 0, 0
	for _, c := range censoredCounts {
		if IsSuppressed(c) {
			censored++
		} else {
			listed++
		}
	}

	mcf := 0.0
	if logLikNull != 0 {
		mcf = 1 - logLikFull/logLikNull
	}
	return SpatialSkill{
		LogLikNull: logLikNull, LogLikFull: logLikFull,
		Alpha: alphaFull, Beta: betaFull,
		LRStatistic: 2 * (logLikFull - logLikNull),
		McFaddenR2:  mcf, IrreducibleFrac: 1 - mcf,
		NListed: listed, NCensored: censored,
	}
}

// TemporalSkill is the national activity persistence decomposition. See the file
// header: MomentumR2 measures transmission persistence, NOT importation skill.
type TemporalSkill struct {
	NWeeks       int
	Lag1Autocorr float64 // persistence of log-incidence at 1 week
	MomentumR2   float64 // expanding-window AR(1) 1-step-ahead skill vs climatology
	Caveat       string
}

// TemporalSkillDecomposition measures short-horizon predictability of the
// national weekly series (onset-ordered counts, right-truncated tail removed by
// the caller). Works on log(x+1).
func TemporalSkillDecomposition(weekly []float64) TemporalSkill {
	n := len(weekly)
	lx := make([]float64, n)
	for i, v := range weekly {
		lx[i] = math.Log(v + 1)
	}
	mean := 0.0
	for _, v := range lx {
		mean += v
	}
	mean /= float64(n)
	varTot := 0.0
	for _, v := range lx {
		varTot += (v - mean) * (v - mean)
	}

	// Lag-1 autocorrelation.
	acf1 := 0.0
	if varTot > 0 {
		num := 0.0
		for i := 1; i < n; i++ {
			num += (lx[i] - mean) * (lx[i-1] - mean)
		}
		acf1 = num / varTot
	}

	// Expanding-window AR(1) one-step-ahead skill vs a climatology baseline.
	sseAR, sseClim := 0.0, 0.0
	start := 20
	if start >= n {
		start = n / 2
	}
	for t := start; t < n; t++ {
		hist := lx[:t]
		xa, ya := hist[:len(hist)-1], hist[1:]
		mh := seriesMean(hist)
		mx, my := seriesMean(xa), seriesMean(ya)
		sxx, sxy := 0.0, 0.0
		for i := range xa {
			sxx += (xa[i] - mx) * (xa[i] - mx)
			sxy += (xa[i] - mx) * (ya[i] - my)
		}
		b := 0.0
		if sxx > 0 {
			b = sxy / sxx
		}
		a := my - b*mx
		predAR := a + b*lx[t-1]
		sseAR += (lx[t] - predAR) * (lx[t] - predAR)
		sseClim += (lx[t] - mh) * (lx[t] - mh)
	}
	momentumR2 := 0.0
	if sseClim > 0 {
		momentumR2 = 1 - sseAR/sseClim
	}

	return TemporalSkill{
		NWeeks: n, Lag1Autocorr: acf1, MomentumR2: momentumR2,
		Caveat: "MomentumR2 reflects transmission persistence (ongoing outbreaks), " +
			"not importation-rate predictability; it is an upper bound on temporal skill.",
	}
}

func seriesMean(x []float64) float64 {
	s := 0.0
	for _, v := range x {
		s += v
	}
	return s / float64(len(x))
}

// goldenSectionMax maximises a unimodal f on [lo, hi] by golden-section search.
func goldenSectionMax(f func(float64) float64, lo, hi float64) float64 {
	const invPhi = 0.6180339887498949
	a, b := lo, hi
	c := b - invPhi*(b-a)
	d := a + invPhi*(b-a)
	fc, fd := f(c), f(d)
	for i := 0; i < 200 && (b-a) > 1e-9; i++ {
		if fc > fd {
			b, d, fd = d, c, fc
			c = b - invPhi*(b-a)
			fc = f(c)
		} else {
			a, c, fc = c, d, fd
			d = a + invPhi*(b-a)
			fd = f(d)
		}
	}
	return 0.5 * (a + b)
}
