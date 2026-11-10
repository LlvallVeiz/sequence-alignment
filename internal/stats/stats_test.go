package stats

import (
	"testing"

	"github.com/LlvallVeiz/sequence-alignment/internal/matrix"
)

// Hand-computed: alignment of "ACGT" vs "A-GT" has 3 matches, 0
// mismatches, 1 gap char in 1 opening, length 4.
func TestComputeDNA(t *testing.T) {
	s := Compute("ACGT", "A-GT", matrix.DNA(1, -1))
	if s.Length != 4 || s.Matches != 3 || s.Mismatches != 0 {
		t.Errorf("stats = %+v", s)
	}
	if s.GapChars != 1 || s.GapOpens != 1 {
		t.Errorf("gaps = %d chars / %d openings, want 1/1", s.GapChars, s.GapOpens)
	}
	if s.IdentityPct != 75.0 {
		t.Errorf("identity = %v, want 75", s.IdentityPct)
	}
	// DNA matrix: similarity equals identity.
	if s.Similar != s.Matches || s.SimilarityPct != s.IdentityPct {
		t.Errorf("DNA similarity should equal identity: %+v", s)
	}
}

// Hand-computed with BLOSUM62: "AWGHE" vs "AW-HE" has 4 matches,
// 0 mismatches, 1 gap char in 1 opening; all 4 matched pairs score > 0
// (A-A 4, W-W 11, H-H 8, E-E 5), so similarity is 4/5 = 80%.
func TestComputeProtein(t *testing.T) {
	m, _ := matrix.Get("blosum62")
	s := Compute("AWGHE", "AW-HE", m)
	if s.Length != 5 || s.Matches != 4 || s.Mismatches != 0 {
		t.Errorf("stats = %+v", s)
	}
	if s.Similar != 4 || s.SimilarityPct != 80.0 {
		t.Errorf("similarity = %d (%.2f%%), want 4 (80%%)", s.Similar, s.SimilarityPct)
	}
	if s.IdentityPct != 80.0 {
		t.Errorf("identity = %v, want 80", s.IdentityPct)
	}
}

// A similar-but-not-identical pair (D-E scores +2 in BLOSUM62) counts
// toward similarity but not identity.
func TestSimilarNotIdentical(t *testing.T) {
	m, _ := matrix.Get("blosum62")
	s := Compute("DE", "ED", m)
	if s.Matches != 0 || s.Mismatches != 2 {
		t.Errorf("stats = %+v", s)
	}
	if s.Similar != 2 {
		t.Errorf("D-E and E-D both score +2, similarity should be 2, got %d", s.Similar)
	}
	if s.IdentityPct != 0 || s.SimilarityPct != 100.0 {
		t.Errorf("identity %.2f%% similarity %.2f%%, want 0%% / 100%%", s.IdentityPct, s.SimilarityPct)
	}
}

// Two separate gap runs count as two openings.
func TestGapRuns(t *testing.T) {
	s := Compute("A--GT--A", "ACCGTTCA", matrix.DNA(1, -1))
	if s.GapChars != 4 || s.GapOpens != 2 {
		t.Errorf("gaps = %d chars / %d openings, want 4/2", s.GapChars, s.GapOpens)
	}
	if s.Length != 8 || s.Matches != 4 {
		t.Errorf("stats = %+v", s)
	}
}

func TestComputeUnequalLengths(t *testing.T) {
	s := Compute("ACGT", "ACG", matrix.DNA(1, -1))
	if s.Length != 0 {
		t.Errorf("unequal lengths should give zero stats, got %+v", s)
	}
}
