// Package matrix provides amino-acid substitution matrices and a simple
// DNA match/mismatch matrix for pairwise sequence alignment.
//
// The built-in protein matrices (BLOSUM45, BLOSUM62, BLOSUM80, PAM250)
// are the canonical NCBI-distributed tables; see data.go for provenance
// and the cross-checks performed while building them.
package matrix

import (
	"fmt"
	"strings"
)

// Matrix is a substitution matrix: Scores[i][j] is the score for aligning
// Alphabet[i] with Alphabet[j]. Protein matrices are symmetric.
type Matrix struct {
	// Name is the short identifier, e.g. "BLOSUM62".
	Name string
	// Description is a one-line human-readable description.
	Description string
	// Alphabet lists the characters this matrix scores, in Scores order.
	Alphabet string
	// Scores is a len(Alphabet) x len(Alphabet) table.
	Scores [][]int
	// IsProtein reports whether the matrix scores amino acids.
	IsProtein bool
}

// index returns the position of c in the alphabet, or -1.
func (m *Matrix) index(c byte) int {
	return strings.IndexByte(m.Alphabet, c)
}

// Score returns the substitution score for aligning a with b.
// The second return value is false when either character is not in the
// matrix alphabet.
func (m *Matrix) Score(a, b byte) (int, bool) {
	ia, ib := m.index(a), m.index(b)
	if ia < 0 || ib < 0 {
		return 0, false
	}
	return m.Scores[ia][ib], true
}

// MustScore is Score but returns an error for characters outside the
// alphabet instead of a boolean.
func (m *Matrix) MustScore(a, b byte) (int, error) {
	s, ok := m.Score(a, b)
	if !ok {
		return 0, fmt.Errorf("matrix %s: character %q not in alphabet %q", m.Name, a, m.Alphabet)
	}
	return s, nil
}

// Valid reports whether every character of s is in the matrix alphabet.
func (m *Matrix) Valid(s string) bool {
	for i := 0; i < len(s); i++ {
		if m.index(s[i]) < 0 {
			return false
		}
	}
	return true
}

// FirstInvalid returns the first character of s that is not in the matrix
// alphabet, and its position, or 0 and -1 when s is fully valid.
func (m *Matrix) FirstInvalid(s string) (byte, int) {
	for i := 0; i < len(s); i++ {
		if m.index(s[i]) < 0 {
			return s[i], i
		}
	}
	return 0, -1
}

// fromTable builds a Matrix from a square array in proteinAlphabet order.
func fromTable(name, desc string, table [24][24]int) *Matrix {
	scores := make([][]int, 24)
	for i := range table {
		row := make([]int, 24)
		copy(row, table[i][:])
		scores[i] = row
	}
	return &Matrix{
		Name:        name,
		Description: desc,
		Alphabet:    proteinAlphabet,
		Scores:      scores,
		IsProtein:   true,
	}
}

// builtins maps lower-cased lookup names to constructors.
var builtins = map[string]func() *Matrix{
	"blosum45": func() *Matrix {
		return fromTable("BLOSUM45",
			"Henikoff & Henikoff (1992), 45% identity clustering; for distant homologs",
			blosum45Data)
	},
	"blosum62": func() *Matrix {
		return fromTable("BLOSUM62",
			"Henikoff & Henikoff (1992), 62% identity clustering; general-purpose default",
			blosum62Data)
	},
	"blosum80": func() *Matrix {
		return fromTable("BLOSUM80",
			"Henikoff & Henikoff (1992), 80% identity clustering; for closely related sequences",
			blosum80Data)
	},
	"pam250": func() *Matrix {
		return fromTable("PAM250",
			"Dayhoff, Schwartz & Orcutt (1978); 250 accepted point mutations per 100 residues",
			pam250Data)
	},
}

// Get returns a built-in matrix by name (case-insensitive, e.g. "BLOSUM62"
// or "pam250"). The special names "dna" and "identity" return DNA scoring
// matrices; see DNA.
func Get(name string) (*Matrix, bool) {
	if c, ok := builtins[strings.ToLower(name)]; ok {
		return c(), true
	}
	switch strings.ToLower(name) {
	case "dna":
		return DNA(1, -1), true
	case "identity":
		return DNA(1, 0), true
	}
	return nil, false
}

// List returns the built-in matrices in a stable order, DNA matrices last.
func List() []*Matrix {
	names := []string{"BLOSUM45", "BLOSUM62", "BLOSUM80", "PAM250"}
	out := make([]*Matrix, 0, len(names)+2)
	for _, n := range names {
		m, _ := Get(n)
		out = append(out, m)
	}
	out = append(out, DNA(1, -1), DNA(1, 0))
	return out
}

// DNA builds a nucleotide scoring matrix over the alphabet ACGTN.
// match is the score for identical bases, mismatch for different bases;
// the ambiguous base N scores 0 against anything (including N).
func DNA(match, mismatch int) *Matrix {
	alpha := "ACGTN"
	scores := make([][]int, len(alpha))
	for i, a := range alpha {
		row := make([]int, len(alpha))
		for j, b := range alpha {
			switch {
			case a == 'N' || b == 'N':
				row[j] = 0
			case a == b:
				row[j] = match
			default:
				row[j] = mismatch
			}
		}
		scores[i] = row
	}
	return &Matrix{
		Name:        fmt.Sprintf("DNA(+%d/%d)", match, mismatch),
		Description: fmt.Sprintf("nucleotide match=%+d mismatch=%+d, N scores 0", match, mismatch),
		Alphabet:    alpha,
		Scores:      scores,
		IsProtein:   false,
	}
}

// IsDNA reports whether s looks like a DNA sequence (only A, C, G, T, N,
// case-insensitive). It is used to pick a sensible default matrix.
func IsDNA(s string) bool {
	if s == "" {
		return false
	}
	for i := 0; i < len(s); i++ {
		switch s[i] {
		case 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n':
		default:
			return false
		}
	}
	return true
}

// Normalize upper-cases s and strips whitespace. It does not validate the
// alphabet; use Valid afterwards.
func Normalize(s string) string {
	var b strings.Builder
	b.Grow(len(s))
	for i := 0; i < len(s); i++ {
		c := s[i]
		switch {
		case c == ' ' || c == '\t' || c == '\n' || c == '\r':
		case 'a' <= c && c <= 'z':
			b.WriteByte(c - 'a' + 'A')
		default:
			b.WriteByte(c)
		}
	}
	return b.String()
}
