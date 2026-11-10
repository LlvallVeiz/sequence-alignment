package matrix

import "testing"

// The expected values below are the NCBI-distributed tables (see data.go
// for provenance). Each entry was cross-checked against at least one
// independent source while building the project; see docs/RESEARCH.md.

func mustGet(t *testing.T, name string) *Matrix {
	t.Helper()
	m, ok := Get(name)
	if !ok {
		t.Fatalf("Get(%q) not found", name)
	}
	return m
}

func TestBLOSUM62KnownEntries(t *testing.T) {
	m := mustGet(t, "BLOSUM62")
	// Diagonal (also listed in the MELO/DeepWiki BLOSUM62 documentation).
	diag := map[byte]int{
		'A': 4, 'R': 5, 'N': 6, 'D': 6, 'C': 9, 'Q': 5, 'E': 5, 'G': 6,
		'H': 8, 'I': 4, 'L': 4, 'K': 5, 'M': 5, 'F': 6, 'P': 7, 'S': 4,
		'T': 5, 'W': 11, 'Y': 7, 'V': 4,
	}
	for a, want := range diag {
		if got, _ := m.Score(a, a); got != want {
			t.Errorf("BLOSUM62 %c-%c = %d, want %d", a, a, got, want)
		}
	}
	// Off-diagonal spot checks from the NCBI file (blosum62.iij).
	cases := []struct {
		a, b byte
		want int
	}{
		{'A', 'W', -3}, {'F', 'Y', 3}, {'H', 'D', -1}, {'I', 'L', 2},
		{'D', 'E', 2}, {'K', 'R', 2}, {'S', 'T', 1}, {'M', 'L', 2},
		{'P', 'T', -1}, {'G', 'A', 0}, {'W', 'Y', 2}, {'C', 'S', -1},
	}
	for _, c := range cases {
		if got, _ := m.Score(c.a, c.b); got != c.want {
			t.Errorf("BLOSUM62 %c-%c = %d, want %d", c.a, c.b, got, c.want)
		}
	}
}

func TestBLOSUM45KnownEntries(t *testing.T) {
	m := mustGet(t, "BLOSUM45")
	// Full 24x24 table verified 0 diffs vs CalState LA blosum45.txt.
	cases := []struct {
		a, b byte
		want int
	}{
		{'A', 'A', 5}, {'C', 'C', 12}, {'W', 'W', 15}, {'H', 'H', 10},
		{'D', 'D', 7}, {'P', 'P', 9}, {'B', 'D', 5}, {'Z', 'E', 4},
		{'X', 'X', -1}, {'*', '*', 1}, {'*', 'A', -5}, {'F', 'Y', 3},
	}
	for _, c := range cases {
		if got, _ := m.Score(c.a, c.b); got != c.want {
			t.Errorf("BLOSUM45 %c-%c = %d, want %d", c.a, c.b, got, c.want)
		}
	}
}

func TestBLOSUM80KnownEntries(t *testing.T) {
	m := mustGet(t, "BLOSUM80")
	// NCBI file blosum80_3.iij (1/3 bit units, entropy 0.9868).
	// NOTE: the `blosum` PyPI package ships a different BLOSUM80 variant
	// (A-A=5, W-W=11); this project uses the NCBI file (A-A=7, W-W=16).
	cases := []struct {
		a, b byte
		want int
	}{
		{'A', 'A', 7}, {'R', 'R', 9}, {'N', 'N', 9}, {'D', 'D', 10},
		{'C', 'C', 13}, {'Q', 'Q', 9}, {'E', 'E', 8}, {'G', 'G', 9},
		{'H', 'H', 12}, {'I', 'I', 7}, {'L', 'L', 6}, {'K', 'K', 8},
		{'M', 'M', 9}, {'F', 'F', 10}, {'P', 'P', 12}, {'S', 'S', 7},
		{'T', 'T', 8}, {'W', 'W', 16}, {'Y', 'Y', 11}, {'V', 'V', 7},
		{'A', 'R', -3}, {'D', 'E', 2}, {'F', 'Y', 4}, {'K', 'R', 3},
	}
	for _, c := range cases {
		if got, _ := m.Score(c.a, c.b); got != c.want {
			t.Errorf("BLOSUM80 %c-%c = %d, want %d", c.a, c.b, got, c.want)
		}
	}
}

func TestPAM250KnownEntries(t *testing.T) {
	m := mustGet(t, "PAM250")
	// Standard Dayhoff PAM250 log-odds table.
	cases := []struct {
		a, b byte
		want int
	}{
		{'A', 'A', 2}, {'C', 'C', 12}, {'W', 'W', 17}, {'S', 'S', 2},
		{'F', 'F', 9}, {'Y', 'Y', 10}, {'R', 'R', 6}, {'G', 'G', 5},
		{'D', 'E', 3}, {'I', 'L', 2}, {'K', 'R', 3}, {'M', 'L', 4},
		{'S', 'T', 1}, {'N', 'D', 2}, {'Q', 'E', 2}, {'V', 'I', 4},
		{'A', 'W', -6}, {'C', 'W', -8}, {'P', 'F', -5},
	}
	for _, c := range cases {
		if got, _ := m.Score(c.a, c.b); got != c.want {
			t.Errorf("PAM250 %c-%c = %d, want %d", c.a, c.b, got, c.want)
		}
	}
}

// TestProteinMatricesSymmetricAndComplete checks every built-in protein
// matrix is square, symmetric and fully populated over its alphabet.
func TestProteinMatricesSymmetricAndComplete(t *testing.T) {
	for _, name := range []string{"BLOSUM45", "BLOSUM62", "BLOSUM80", "PAM250"} {
		m := mustGet(t, name)
		n := len(m.Alphabet)
		if len(m.Scores) != n {
			t.Fatalf("%s: %d rows for %d-letter alphabet", name, len(m.Scores), n)
		}
		for i := 0; i < n; i++ {
			if len(m.Scores[i]) != n {
				t.Fatalf("%s: row %d has %d columns, want %d", name, i, len(m.Scores[i]), n)
			}
			for j := 0; j < n; j++ {
				if m.Scores[i][j] != m.Scores[j][i] {
					t.Errorf("%s not symmetric at %c-%c: %d vs %d",
						name, m.Alphabet[i], m.Alphabet[j], m.Scores[i][j], m.Scores[j][i])
				}
			}
		}
	}
}

func TestGetUnknown(t *testing.T) {
	if _, ok := Get("BLOSUM99"); ok {
		t.Error("Get(BLOSUM99) should fail")
	}
	if _, ok := Get(""); ok {
		t.Error("Get(\"\") should fail")
	}
}

func TestDNA(t *testing.T) {
	m := DNA(2, -3)
	if got, _ := m.Score('A', 'A'); got != 2 {
		t.Errorf("DNA A-A = %d, want 2", got)
	}
	if got, _ := m.Score('A', 'C'); got != -3 {
		t.Errorf("DNA A-C = %d, want -3", got)
	}
	// Ambiguous base N scores 0 against anything, including itself.
	for _, b := range "ACGTN" {
		if got, _ := m.Score('N', byte(b)); got != 0 {
			t.Errorf("DNA N-%c = %d, want 0", b, got)
		}
	}
	if _, ok := m.Score('A', 'X'); ok {
		t.Error("DNA matrix should reject X")
	}
	if !IsDNA("acgtnACGTN") {
		t.Error("IsDNA should accept ACGTN in any case")
	}
	if IsDNA("ACGTU") || IsDNA("ACGT-") || IsDNA("") {
		t.Error("IsDNA should reject U, '-' and empty input")
	}
}

func TestNormalize(t *testing.T) {
	if got := Normalize("ac gt\nAC"); got != "ACGTAC" {
		t.Errorf("Normalize = %q, want %q", got, "ACGTAC")
	}
}
