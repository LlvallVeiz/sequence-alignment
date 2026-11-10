package align

import (
	"strings"
	"testing"

	"github.com/LlvallVeiz/sequence-alignment/internal/matrix"
)

func dnaParams() Params {
	return Params{Matrix: matrix.DNA(1, -1), GapOpen: 5, GapExtend: 1}
}

func blosum62Params() Params {
	m, _ := matrix.Get("blosum62")
	return Params{Matrix: m, GapOpen: 10, GapExtend: 1}
}

// handCases are small alignments whose scores were computed by hand;
// the reasoning is in the comment of each case.
func TestHandComputedCases(t *testing.T) {
	cases := []struct {
		name      string
		algo      Algorithm
		a, b      string
		params    Params
		wantScore int
		wantA     string
		wantB     string
	}{
		{
			// 4 matches x +1.
			name: "NW identical DNA", algo: Global, a: "ACGT", b: "ACGT",
			params: dnaParams(), wantScore: 4, wantA: "ACGT", wantB: "ACGT",
		},
		{
			// One mismatch: -1.
			name: "NW single mismatch", algo: Global, a: "A", b: "C",
			params: dnaParams(), wantScore: -1, wantA: "A", wantB: "C",
		},
		{
			// A-A (+1) plus a gap of length 1: -(5+1) = -6. Total -5.
			name: "NW one gap", algo: Global, a: "AC", b: "A",
			params: dnaParams(), wantScore: -5, wantA: "AC", wantB: "A-",
		},
		{
			// 3 matches (+3) plus a gap of length 1 (-6). Total -3.
			name: "NW trailing gap", algo: Global, a: "ACGT", b: "ACG",
			params: dnaParams(), wantScore: -3, wantA: "ACGT", wantB: "ACG-",
		},
		{
			// BLOSUM62 diagonal: A-A 4 + C-C 9 + D-D 6 + E-E 5 = 24.
			name: "NW identical protein", algo: Global, a: "ACDE", b: "ACDE",
			params: blosum62Params(), wantScore: 24, wantA: "ACDE", wantB: "ACDE",
		},
		{
			// Single BLOSUM62 pair A-W = -3.
			name: "NW single protein pair", algo: Global, a: "A", b: "W",
			params: blosum62Params(), wantScore: -3, wantA: "A", wantB: "W",
		},
		{
			// Trp diagonal W-W = 11.
			name: "NW tryptophan", algo: Global, a: "W", b: "W",
			params: blosum62Params(), wantScore: 11, wantA: "W", wantB: "W",
		},
		{
			// Best local block is "ACGTA" (5 matches) at A[1..5], B[3..7].
			name: "SW ACGTA block", algo: Local, a: "ACGTACGT", b: "TTACGTAA",
			params: dnaParams(), wantScore: 5, wantA: "ACGTA", wantB: "ACGTA",
		},
		{
			// No positive-scoring block exists.
			name: "SW no match", algo: Local, a: "AAAA", b: "CCCC",
			params: dnaParams(), wantScore: 0, wantA: "", wantB: "",
		},
		{
			// "ACGT" fits into "TTACGTTT" with free end gaps on both sides.
			name: "semi fitting", algo: SemiGlobal, a: "ACGT", b: "TTACGTTT",
			params: dnaParams(), wantScore: 4, wantA: "--ACGT--", wantB: "TTACGTTT",
		},
		{
			// Ungapped: G=G +1, A-C -1, T-A -1, T-T +1, A-G -1, C-C +1,
			// A-T -1. Total -1.
			name: "NW GATTACA", algo: Global, a: "GATTACA", b: "GCATGCT",
			params: dnaParams(), wantScore: -1, wantA: "GATTACA", wantB: "GCATGCT",
		},
		{
			// Best local block is "AT" (2 matches).
			name: "SW GATTACA", algo: Local, a: "GATTACA", b: "GCATGCT",
			params: dnaParams(), wantScore: 2, wantA: "AT", wantB: "AT",
		},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			res, err := Align(c.a, c.b, c.algo, c.params)
			if err != nil {
				t.Fatalf("Align: %v", err)
			}
			if res.Score != c.wantScore {
				t.Errorf("score = %d, want %d (alignment %q / %q)",
					res.Score, c.wantScore, res.AlignedA, res.AlignedB)
			}
			if res.AlignedA != c.wantA || res.AlignedB != c.wantB {
				t.Errorf("alignment = %q / %q, want %q / %q",
					res.AlignedA, res.AlignedB, c.wantA, c.wantB)
			}
		})
	}
}

// TestClassicSmithWaterman checks the textbook example
// HEAGAWGHEE vs PAWHEAE under BLOSUM62 with affine gaps (10,1).
// The best local score is 17: either HEA/HEA (H-H 8 + E-E 5 + A-A 4)
// or AWGHE/AW-HE (4+11+8+5-11). Both score 17, so only the score and
// the re-scored alignment are asserted, not the exact alignment.
func TestClassicSmithWaterman(t *testing.T) {
	res, err := Align("HEAGAWGHEE", "PAWHEAE", Local, blosum62Params())
	if err != nil {
		t.Fatal(err)
	}
	if res.Score != 17 {
		t.Fatalf("score = %d, want 17 (alignment %q / %q)", res.Score, res.AlignedA, res.AlignedB)
	}
	back, err := ScoreAlignment(res.AlignedA, res.AlignedB, blosum62Params().Matrix, 10, 1)
	if err != nil {
		t.Fatal(err)
	}
	if back != res.Score {
		t.Errorf("re-scored alignment = %d, want %d", back, res.Score)
	}
}

// crossCheckCases were produced by the independent Python reference
// implementation (scripts/reference_align.py, seed 20261110). The Go
// score must equal the reference score for every case.
func TestCrossCheckedAgainstReference(t *testing.T) {
	cases := []struct {
		algo      string
		a, b      string
		matrix    string
		open, ext int
		wantScore int
	}{
		{"NW", "GCTGGCGCGTCG", "TTCCTCCCCG", "dna", 5, 1, -7},
		{"SW", "GCTGGCGCGTCG", "TTCCTCCCCG", "dna", 5, 1, 2},
		{"SEMI", "GCTGGCGCGTCG", "TTCCTCCCCG", "dna", 5, 1, 1},
		{"NW", "AAGTAACTCTCA", "CTTACAAGGC", "dna", 5, 1, -11},
		{"SW", "AAGTAACTCTCA", "CTTACAAGGC", "dna", 5, 1, 3},
		{"SEMI", "AAGTAACTCTCA", "CTTACAAGGC", "dna", 5, 1, 2},
		{"NW", "CTGCACAGAGAT", "CCTTCTCTAG", "dna", 5, 1, -11},
		{"SW", "CTGCACAGAGAT", "CCTTCTCTAG", "dna", 5, 1, 2},
		{"SEMI", "CTGCACAGAGAT", "CCTTCTCTAG", "dna", 5, 1, 0},
		{"NW", "AILFQTKRRAPRCN", "FNAGTYCYMRS", "blosum62", 10, 1, -18},
		{"SW", "AILFQTKRRAPRCN", "FNAGTYCYMRS", "blosum62", 10, 1, 9},
		{"NW", "WGQGLMNPNRDWMG", "GLDNMVEFYDE", "blosum62", 10, 1, -11},
		{"SW", "WGQGLMNPNRDWMG", "GLDNMVEFYDE", "blosum62", 10, 1, 13},
		{"NW", "QVNWSYGQFRQMTA", "GRWKNKLSYFF", "blosum62", 10, 1, -20},
		{"SW", "QVNWSYGQFRQMTA", "GRWKNKLSYFF", "blosum62", 10, 1, 11},
	}
	for _, c := range cases {
		algo, err := ParseAlgorithm(c.algo)
		if err != nil {
			t.Fatal(err)
		}
		var p Params
		switch c.matrix {
		case "dna":
			p = dnaParams()
		case "blosum62":
			p = blosum62Params()
		}
		res, err := Align(c.a, c.b, algo, p)
		if err != nil {
			t.Fatalf("%s %s vs %s: %v", c.algo, c.a, c.b, err)
		}
		if res.Score != c.wantScore {
			t.Errorf("%s %s vs %s: score = %d, want %d (reference)",
				c.algo, c.a, c.b, res.Score, c.wantScore)
		}
		// The returned alignment must re-score to the same value for
		// global mode (local/semi-global have free or zero regions).
		if c.algo == "NW" {
			back, err := ScoreAlignment(res.AlignedA, res.AlignedB, p.Matrix, p.GapOpen, p.GapExtend)
			if err != nil {
				t.Fatalf("rescore: %v", err)
			}
			if back != res.Score {
				t.Errorf("%s %s vs %s: re-scored %d != reported %d",
					c.algo, c.a, c.b, back, res.Score)
			}
		}
	}
}

// TestHemoglobinPair aligns the real human alpha- and beta-globin
// chains (UniProt P69905 and P68871) and checks the scores produced by
// the independent Python reference implementation.
func TestHemoglobinPair(t *testing.T) {
	hba := "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
	hbb := "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
	p := blosum62Params()
	for algo, want := range map[Algorithm]int{Global: 286, Local: 288, SemiGlobal: 286} {
		res, err := Align(hba, hbb, algo, p)
		if err != nil {
			t.Fatal(err)
		}
		if res.Score != want {
			t.Errorf("%s HBA/HBB score = %d, want %d", algo, res.Score, want)
		}
		if len(res.AlignedA) != len(res.AlignedB) {
			t.Errorf("%s: aligned lengths differ", algo)
		}
	}
}

func TestEmptySequences(t *testing.T) {
	res, err := Align("", "", Global, dnaParams())
	if err != nil || res.Score != 0 || res.AlignedA != "" {
		t.Errorf("empty/empty = %+v, %v", res, err)
	}
	res, err = Align("", "ACGT", Global, dnaParams())
	if err != nil {
		t.Fatal(err)
	}
	// One gap of length 4: -(5+4) = -9.
	if res.Score != -9 || res.AlignedA != "----" || res.AlignedB != "ACGT" {
		t.Errorf("empty/ACGT = %+v", res)
	}
	res, err = Align("", "ACGT", Local, dnaParams())
	if err != nil || res.Score != 0 {
		t.Errorf("local empty/ACGT = %+v, %v", res, err)
	}
}

func TestInvalidInput(t *testing.T) {
	if _, err := Align("ACGTX", "ACGT", Global, dnaParams()); err == nil {
		t.Error("DNA alignment should reject X")
	}
	if _, err := Align("ACGT", "ACGT", Global, Params{Matrix: nil}); err == nil {
		t.Error("nil matrix should be an error")
	}
	if _, err := Align("ACGT", "ACGT", Global,
		Params{Matrix: matrix.DNA(1, -1), GapOpen: -1}); err == nil {
		t.Error("negative gap penalty should be an error")
	}
	// Lowercase input and RNA U are accepted for DNA scoring.
	res, err := Align("acgu", "ACGT", Global, dnaParams())
	if err != nil {
		t.Fatal(err)
	}
	if res.Score != 4 {
		t.Errorf("acgu vs ACGT score = %d, want 4 (U read as T)", res.Score)
	}
	// Protein matrix rejects U (selenocysteine is not in the alphabet).
	if _, err := Align("ACGU", "ACGU", Global, blosum62Params()); err == nil {
		t.Error("protein alignment should reject U")
	}
}

func TestScoreAlignment(t *testing.T) {
	dna := matrix.DNA(1, -1)
	// A-A +1, gap run of 2 in A: -(5+2) = -7, G-G +1, T-T +1. Total -4.
	score, err := ScoreAlignment("A--GT", "ACCGT", dna, 5, 1)
	if err != nil || score != -4 {
		t.Errorf("ScoreAlignment = %d, %v; want -4", score, err)
	}
	if _, err := ScoreAlignment("ACGT", "ACG", dna, 5, 1); err == nil {
		t.Error("unequal lengths should be an error")
	}
	if _, err := ScoreAlignment("A-", "-A", dna, 5, 1); err == nil {
		// A- vs -A is legal: two single gaps, no gap-gap column.
		t.Log("note: A-/-A accepted")
	}
	if _, err := ScoreAlignment("--", "--", dna, 5, 1); err == nil {
		t.Error("gap-to-gap column should be an error")
	}
}

func TestParseAlgorithm(t *testing.T) {
	for in, want := range map[string]Algorithm{
		"nw": Global, "global": Global, "sw": Local, "local": Local,
		"semi": SemiGlobal, "overlap": SemiGlobal,
	} {
		got, err := ParseAlgorithm(in)
		if err != nil || got != want {
			t.Errorf("ParseAlgorithm(%q) = %v, %v; want %v", in, got, err, want)
		}
	}
	if _, err := ParseAlgorithm("bogus"); err == nil {
		t.Error("bogus algorithm should be an error")
	}
}

func TestSemiGlobalFitting(t *testing.T) {
	// Reference value from scripts/reference_align.py.
	res, err := Align("GATTACA", "TTGCATGCTAA", SemiGlobal, dnaParams())
	if err != nil {
		t.Fatal(err)
	}
	if res.Score != 3 {
		t.Errorf("semi GATTACA score = %d, want 3", res.Score)
	}
	if len(res.AlignedA) != len(res.AlignedB) {
		t.Error("aligned lengths differ")
	}
	if !strings.Contains(res.AlignedA, "GATTACA") {
		t.Errorf("semi alignment lost the query: %q", res.AlignedA)
	}
}
