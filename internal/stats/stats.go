// Package stats computes descriptive statistics of a pairwise alignment:
// score, identity, similarity, gaps.
package stats

import (
	"github.com/LlvallVeiz/sequence-alignment/internal/matrix"
)

// Stats describes one alignment.
type Stats struct {
	// Score is the alignment score (0 when no matrix was supplied).
	Score int
	// Length is the number of alignment columns.
	Length int
	// Matches is the number of columns with identical residues.
	Matches int
	// Mismatches is the number of columns with two different residues.
	Mismatches int
	// Similar is the number of columns whose substitution score is
	// positive (includes identities). For a DNA matrix this equals Matches.
	Similar int
	// GapChars is the total number of '-' characters.
	GapChars int
	// GapOpens is the number of gap runs.
	GapOpens int
	// IdentityPct is 100*Matches/Length.
	IdentityPct float64
	// SimilarityPct is 100*Similar/Length.
	SimilarityPct float64
}

// Compute derives statistics from two aligned strings of equal length.
// m may be nil, in which case Score stays 0 and Similar equals Matches.
func Compute(alignedA, alignedB string, m *matrix.Matrix) Stats {
	var s Stats
	if len(alignedA) != len(alignedB) {
		return s // zero value; callers validate lengths first
	}
	s.Length = len(alignedA)
	inGapA, inGapB := false, false
	for k := 0; k < len(alignedA); k++ {
		ca, cb := alignedA[k], alignedB[k]
		switch {
		case ca == '-' && cb == '-':
			// Not a valid alignment column; ignore it.
		case ca == '-':
			s.GapChars++
			if !inGapA {
				s.GapOpens++
				inGapA = true
			}
			inGapB = false
		case cb == '-':
			s.GapChars++
			if !inGapB {
				s.GapOpens++
				inGapB = true
			}
			inGapA = false
		default:
			inGapA, inGapB = false, false
			if ca == cb {
				s.Matches++
			} else {
				s.Mismatches++
			}
			if m != nil {
				if v, ok := m.Score(ca, cb); ok && v > 0 {
					s.Similar++
				}
			} else if ca == cb {
				s.Similar++
			}
		}
	}
	if s.Length > 0 {
		s.IdentityPct = 100 * float64(s.Matches) / float64(s.Length)
		s.SimilarityPct = 100 * float64(s.Similar) / float64(s.Length)
	}
	return s
}
