// Package align implements pairwise sequence alignment by dynamic
// programming: Needleman-Wunsch (global), Smith-Waterman (local) and
// semi-global/overlap alignment, all with affine gap penalties following
// Gotoh (1982).
//
// Gap model: a gap of length k costs GapOpen + k*GapExtend, where both
// penalties are given as positive magnitudes. With GapOpen=10 and
// GapExtend=1 this matches the NCBI BLAST convention of gap costs (10,1),
// where a gap of length 1 costs 11.
package align

import (
	"errors"
	"fmt"
	"strings"

	"github.com/LlvallVeiz/sequence-alignment/internal/matrix"
)

// Algorithm selects the alignment mode.
type Algorithm int

const (
	// Global aligns the full length of both sequences (Needleman-Wunsch).
	Global Algorithm = iota
	// Local finds the best-scoring substring pair (Smith-Waterman).
	Local
	// SemiGlobal aligns with free end gaps on both sequences
	// (overlap alignment): leading and trailing gaps cost nothing.
	SemiGlobal
)

// String returns the short name used on the command line.
func (a Algorithm) String() string {
	switch a {
	case Global:
		return "global"
	case Local:
		return "local"
	case SemiGlobal:
		return "semi-global"
	default:
		return "unknown"
	}
}

// FullName returns the eponymous name of the algorithm.
func (a Algorithm) FullName() string {
	switch a {
	case Global:
		return "Needleman-Wunsch (global)"
	case Local:
		return "Smith-Waterman (local)"
	case SemiGlobal:
		return "Semi-global (overlap)"
	default:
		return "unknown"
	}
}

// ParseAlgorithm parses "nw", "global", "sw", "local", "semi", "overlap"
// (case-insensitive).
func ParseAlgorithm(s string) (Algorithm, error) {
	switch strings.ToLower(s) {
	case "nw", "global", "needleman-wunsch", "needlemanwunsch":
		return Global, nil
	case "sw", "local", "smith-waterman", "smithwaterman":
		return Local, nil
	case "semi", "semi-global", "semiglobal", "overlap", "fitting":
		return SemiGlobal, nil
	}
	return Global, fmt.Errorf("unknown algorithm %q (want nw/global, sw/local or semi/overlap)", s)
}

// Params holds the scoring configuration for an alignment.
type Params struct {
	// Matrix scores aligned residue pairs.
	Matrix *matrix.Matrix
	// GapOpen is the positive penalty for opening a gap.
	GapOpen int
	// GapExtend is the positive penalty for each residue in a gap.
	GapExtend int
}

// Result is one pairwise alignment.
type Result struct {
	// AlignedA and AlignedB have equal length; '-' marks a gap.
	AlignedA string
	AlignedB string
	// Score is the alignment score under Params.
	Score int
}

// Align aligns a against b. Inputs are upper-cased and whitespace is
// stripped; characters outside the matrix alphabet are an error.
func Align(a, b string, algo Algorithm, p Params) (Result, error) {
	if p.Matrix == nil {
		return Result{}, errors.New("align: nil scoring matrix")
	}
	if p.GapOpen < 0 || p.GapExtend < 0 {
		return Result{}, fmt.Errorf("align: gap penalties must be non-negative (got open=%d extend=%d)", p.GapOpen, p.GapExtend)
	}
	a = matrix.Normalize(a)
	b = matrix.Normalize(b)
	if !p.Matrix.IsProtein {
		// Accept RNA input for DNA scoring: U is read as T.
		a = strings.ReplaceAll(a, "U", "T")
		b = strings.ReplaceAll(b, "U", "T")
	}
	if c, pos := p.Matrix.FirstInvalid(a); pos >= 0 {
		return Result{}, fmt.Errorf("align: sequence A: character %q at position %d not in %s alphabet", c, pos+1, p.Matrix.Name)
	}
	if c, pos := p.Matrix.FirstInvalid(b); pos >= 0 {
		return Result{}, fmt.Errorf("align: sequence B: character %q at position %d not in %s alphabet", c, pos+1, p.Matrix.Name)
	}
	m, n := len(a), len(b)

	// Degenerate cases with an empty sequence.
	if m == 0 && n == 0 {
		return Result{Score: 0}, nil
	}
	if m == 0 || n == 0 {
		switch algo {
		case Local:
			return Result{Score: 0}, nil
		default:
			l := m + n
			gap := gapCost(l, p.GapOpen, p.GapExtend)
			if m == 0 {
				return Result{AlignedA: strings.Repeat("-", n), AlignedB: b, Score: -gap}, nil
			}
			return Result{AlignedA: a, AlignedB: strings.Repeat("-", m), Score: -gap}, nil
		}
	}

	eng := newEngine(a, b, p)
	switch algo {
	case Global:
		return eng.global(), nil
	case Local:
		return eng.local(), nil
	case SemiGlobal:
		return eng.semiGlobal(), nil
	}
	return Result{}, fmt.Errorf("align: unknown algorithm %v", algo)
}

func gapCost(k, open, extend int) int { return open + k*extend }

// negInf is smaller than any reachable alignment score.
const negInf = -1 << 28

const (
	stM byte = 'M' // ends with an aligned pair (diagonal step)
	stX byte = 'X' // ends with a gap in A (consumes a residue of B)
	stY byte = 'Y' // ends with a gap in B (consumes a residue of A)
)

// engine holds the three Gotoh DP matrices and their traceback pointers.
type engine struct {
	a, b string
	m, n int
	p    Params
	// M, X, Y are (m+1) x (n+1) score matrices; tM, tX, tY record, for
	// each cell, which matrix its value came from.
	M, X, Y    [][]int
	tM, tX, tY [][]byte
	sub        [][]int // sub[i][j] = score(a[i], b[j]) for i,j >= 1
}

func newEngine(a, b string, p Params) *engine {
	m, n := len(a), len(b)
	e := &engine{a: a, b: b, m: m, n: n, p: p}
	e.M = newIntMatrix(m+1, n+1)
	e.X = newIntMatrix(m+1, n+1)
	e.Y = newIntMatrix(m+1, n+1)
	e.tM = newByteMatrix(m+1, n+1)
	e.tX = newByteMatrix(m+1, n+1)
	e.tY = newByteMatrix(m+1, n+1)
	e.sub = newIntMatrix(m+1, n+1)
	for i := 1; i <= m; i++ {
		for j := 1; j <= n; j++ {
			s, _ := p.Matrix.Score(a[i-1], b[j-1])
			e.sub[i][j] = s
		}
	}
	return e
}

func newIntMatrix(rows, cols int) [][]int {
	mx := make([][]int, rows)
	for i := range mx {
		mx[i] = make([]int, cols)
	}
	return mx
}

func newByteMatrix(rows, cols int) [][]byte {
	mx := make([][]byte, rows)
	for i := range mx {
		mx[i] = make([]byte, cols)
	}
	return mx
}

func max3(x, y, z int) (int, byte) {
	// Tie-break: prefer M, then X, then Y (substitution over gaps).
	if x >= y && x >= z {
		return x, stM
	}
	if y >= z {
		return y, stX
	}
	return z, stY
}

// fill runs the DP. local floors every cell at zero; semi frees end gaps.
func (e *engine) fill(algo Algorithm) {
	d, ex := e.p.GapOpen, e.p.GapExtend
	m, n := e.m, e.n

	e.M[0][0] = 0
	e.X[0][0] = negInf
	e.Y[0][0] = negInf

	switch algo {
	case Local:
		// A local alignment may start anywhere: every boundary cell
		// costs 0, so a leading gap is never charged.
		for i := 1; i <= m; i++ {
			e.M[i][0] = 0
			e.Y[i][0] = 0
			e.X[i][0] = 0
			e.tM[i][0] = stM
			e.tY[i][0] = stM
		}
		for j := 1; j <= n; j++ {
			e.M[0][j] = 0
			e.X[0][j] = 0
			e.Y[0][j] = 0
			e.tM[0][j] = stM
			e.tX[0][j] = stM
		}
	case SemiGlobal:
		// Leading gaps are free: first row/column of M cost 0, and a
		// leading gap run in either sequence costs 0 as well.
		for i := 1; i <= m; i++ {
			e.M[i][0] = 0
			e.Y[i][0] = 0
			e.X[i][0] = negInf
			e.tM[i][0] = stM
			if i == 1 {
				e.tY[i][0] = stM
			} else {
				e.tY[i][0] = stY
			}
		}
		for j := 1; j <= n; j++ {
			e.M[0][j] = 0
			e.X[0][j] = 0
			e.Y[0][j] = negInf
			e.tM[0][j] = stM
			if j == 1 {
				e.tX[0][j] = stM
			} else {
				e.tX[0][j] = stX
			}
		}
	default:
		for i := 1; i <= m; i++ {
			e.M[i][0] = -(d + i*ex)
			e.Y[i][0] = -(d + i*ex)
			e.X[i][0] = negInf
			e.tM[i][0] = stM
			if i == 1 {
				e.tY[i][0] = stM
			} else {
				e.tY[i][0] = stY
			}
		}
		for j := 1; j <= n; j++ {
			e.M[0][j] = -(d + j*ex)
			e.X[0][j] = -(d + j*ex)
			e.Y[0][j] = negInf
			e.tM[0][j] = stM
			if j == 1 {
				e.tX[0][j] = stM
			} else {
				e.tX[0][j] = stX
			}
		}
	}

	for i := 1; i <= m; i++ {
		for j := 1; j <= n; j++ {
			// M[i][j]: diagonal step.
			best, from := max3(e.M[i-1][j-1], e.X[i-1][j-1], e.Y[i-1][j-1])
			val := best + e.sub[i][j]
			if algo == Local && val < 0 {
				val = 0
				from = 0
			}
			e.M[i][j] = val
			e.tM[i][j] = from
			// X[i][j]: gap in A (consume b[j]).
			open := e.M[i][j-1] - (d + ex)
			ext := e.X[i][j-1] - ex
			if ext >= open { // tie: extend the open gap rather than opening a new one
				e.X[i][j] = ext
				e.tX[i][j] = stX
			} else {
				e.X[i][j] = open
				e.tX[i][j] = stM
			}
			if algo == Local && e.X[i][j] < 0 {
				e.X[i][j] = 0
				e.tX[i][j] = 0
			}
			// Y[i][j]: gap in B (consume a[i]).
			open = e.M[i-1][j] - (d + ex)
			ext = e.Y[i-1][j] - ex
			if ext >= open {
				e.Y[i][j] = ext
				e.tY[i][j] = stY
			} else {
				e.Y[i][j] = open
				e.tY[i][j] = stM
			}
			if algo == Local && e.Y[i][j] < 0 {
				e.Y[i][j] = 0
				e.tY[i][j] = 0
			}
		}
	}
}

// traceback walks from (i, j) in the given state back to the start,
// returning the alignment built in reverse. It stops early when stopAtZero
// is set and the current cell scores 0 (local alignment).
func (e *engine) traceback(i, j int, state byte, stopAtZero bool) (string, string) {
	var ra, rb []byte
	for i > 0 || j > 0 {
		switch state {
		case stM:
			if stopAtZero && e.M[i][j] == 0 {
				return reverseString(ra), reverseString(rb)
			}
			if i == 0 { // leading gap in A
				ra = append(ra, '-')
				rb = append(rb, e.b[j-1])
				j--
				continue
			}
			if j == 0 { // leading gap in B
				ra = append(ra, e.a[i-1])
				rb = append(rb, '-')
				i--
				continue
			}
			ra = append(ra, e.a[i-1])
			rb = append(rb, e.b[j-1])
			state = e.tM[i][j]
			i--
			j--
		case stX:
			if stopAtZero && e.X[i][j] == 0 {
				return reverseString(ra), reverseString(rb)
			}
			ra = append(ra, '-')
			rb = append(rb, e.b[j-1])
			state = e.tX[i][j]
			j--
		case stY:
			if stopAtZero && e.Y[i][j] == 0 {
				return reverseString(ra), reverseString(rb)
			}
			ra = append(ra, e.a[i-1])
			rb = append(rb, '-')
			state = e.tY[i][j]
			i--
		}
	}
	return reverseString(ra), reverseString(rb)
}

func reverseString(b []byte) string {
	for i, j := 0, len(b)-1; i < j; i, j = i+1, j-1 {
		b[i], b[j] = b[j], b[i]
	}
	return string(b)
}

func (e *engine) global() Result {
	e.fill(Global)
	score, state := max3(e.M[e.m][e.n], e.X[e.m][e.n], e.Y[e.m][e.n])
	aa, bb := e.traceback(e.m, e.n, state, false)
	return Result{AlignedA: aa, AlignedB: bb, Score: score}
}

func (e *engine) local() Result {
	e.fill(Local)
	best, bi, bj, state := 0, 0, 0, byte(0)
	for i := 1; i <= e.m; i++ {
		for j := 1; j <= e.n; j++ {
			for _, cell := range []struct {
				v int
				s byte
			}{{e.M[i][j], stM}, {e.X[i][j], stX}, {e.Y[i][j], stY}} {
				if cell.v > best {
					best, bi, bj, state = cell.v, i, j, cell.s
				}
			}
		}
	}
	if best == 0 {
		return Result{Score: 0}
	}
	aa, bb := e.traceback(bi, bj, state, true)
	return Result{AlignedA: aa, AlignedB: bb, Score: best}
}

func (e *engine) semiGlobal() Result {
	e.fill(SemiGlobal)
	// Best end cell on the last row or last column (trailing gaps free).
	// Scan the last row left-to-right, then the last column top-to-bottom;
	// strict > keeps the leftmost/topmost optimum.
	best, bi, bj := e.M[e.m][0], e.m, 0
	for j := 1; j <= e.n; j++ {
		if e.M[e.m][j] > best {
			best, bi, bj = e.M[e.m][j], e.m, j
		}
	}
	for i := 0; i <= e.m; i++ {
		if e.M[i][e.n] > best {
			best, bi, bj = e.M[i][e.n], i, e.n
		}
	}
	aa, bb := e.traceback(bi, bj, stM, false)
	// Append the free trailing overhang of whichever sequence ran past
	// the end cell.
	if bj < e.n {
		aa += strings.Repeat("-", e.n-bj)
		bb += e.b[bj:]
	}
	if bi < e.m {
		aa += e.a[bi:]
		bb += strings.Repeat("-", e.m-bi)
	}
	return Result{AlignedA: aa, AlignedB: bb, Score: best}
}

// ScoreAlignment scores two already-aligned strings (equal length, '-'
// for gaps) under the given matrix and affine gap penalties. It is the
// inverse of Align: the score of Align's output re-computed here must
// match, except for local and semi-global alignments, whose free/zero
// regions this function scores literally.
func ScoreAlignment(alignedA, alignedB string, m *matrix.Matrix, gapOpen, gapExtend int) (int, error) {
	if len(alignedA) != len(alignedB) {
		return 0, fmt.Errorf("score: aligned sequences have different lengths (%d vs %d)", len(alignedA), len(alignedB))
	}
	if gapOpen < 0 || gapExtend < 0 {
		return 0, fmt.Errorf("score: gap penalties must be non-negative (got open=%d extend=%d)", gapOpen, gapExtend)
	}
	score := 0
	inGapA, inGapB := false, false
	for k := 0; k < len(alignedA); k++ {
		ca, cb := alignedA[k], alignedB[k]
		switch {
		case ca == '-' && cb == '-':
			return 0, fmt.Errorf("score: column %d aligns gap to gap", k+1)
		case ca == '-':
			if inGapA {
				score -= gapExtend
			} else {
				score -= gapOpen + gapExtend
				inGapA = true
			}
			inGapB = false
		case cb == '-':
			if inGapB {
				score -= gapExtend
			} else {
				score -= gapOpen + gapExtend
				inGapB = true
			}
			inGapA = false
		default:
			s, err := m.MustScore(ca, cb)
			if err != nil {
				return 0, fmt.Errorf("score: column %d: %w", k+1, err)
			}
			score += s
			inGapA, inGapB = false, false
		}
	}
	return score, nil
}
