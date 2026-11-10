// Command seqalign aligns pairs of biological sequences.
//
// Usage:
//
//	seqalign align [flags] SEQ1 SEQ2     pairwise alignment
//	seqalign batch -in pairs.fasta ...    align many pairs from one FASTA file
//	seqalign matrices [-name BLOSUM62]    list (and print) scoring matrices
//	seqalign score [flags] ALN1 ALN2      score two already-aligned sequences
//
// A sequence argument is either a raw sequence (ACGT..., MVLSPAD...)
// or "@path/to/file.fasta", which reads the first record of a FASTA file.
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/LlvallVeiz/sequence-alignment/internal/align"
	"github.com/LlvallVeiz/sequence-alignment/internal/fasta"
	"github.com/LlvallVeiz/sequence-alignment/internal/matrix"
	"github.com/LlvallVeiz/sequence-alignment/internal/stats"
)

const version = "1.0.0"

func main() {
	if len(os.Args) < 2 {
		usage(os.Stdout)
		os.Exit(2)
	}
	var err error
	switch os.Args[1] {
	case "align":
		err = cmdAlign(os.Args[2:])
	case "batch":
		err = cmdBatch(os.Args[2:])
	case "matrices":
		err = cmdMatrices(os.Args[2:])
	case "score":
		err = cmdScore(os.Args[2:])
	case "-version", "--version", "version":
		fmt.Printf("seqalign %s\n", version)
		return
	case "-h", "--help", "help":
		usage(os.Stdout)
		return
	default:
		fmt.Fprintf(os.Stderr, "seqalign: unknown subcommand %q\n\n", os.Args[1])
		usage(os.Stderr)
		os.Exit(2)
	}
	if err != nil {
		fmt.Fprintf(os.Stderr, "seqalign: %v\n", err)
		os.Exit(1)
	}
}

func usage(w io.Writer) {
	fmt.Fprint(w, `seqalign — pairwise sequence alignment (Needleman-Wunsch, Smith-Waterman, semi-global)

Usage:
  seqalign align [flags] SEQ1 SEQ2        align two sequences
  seqalign batch -in pairs.fasta [flags]   align consecutive FASTA record pairs
  seqalign matrices [-name NAME]           list scoring matrices (or print one)
  seqalign score [flags] ALN1 ALN2         score two pre-aligned sequences
  seqalign -version                         print version

Sequence arguments are raw sequences, or @file.fasta to read the first
record of a FASTA file. Flags come before the sequence arguments.
Run 'seqalign <cmd> -h' for the flags of a command.
`)
}

// resolveSeq turns a CLI sequence argument into a sequence string and a
// display ID. "@file" reads the first FASTA record of file; a raw
// sequence gets the fallback ID ("seq1"/"seq2").
func resolveSeq(arg, fallback string) (seq, id string, err error) {
	if strings.HasPrefix(arg, "@") {
		path := arg[1:]
		rec, err := fasta.First(path)
		if err != nil {
			return "", "", err
		}
		return rec.Sequence, rec.ID, nil
	}
	return matrix.Normalize(arg), fallback, nil
}

// pickMatrix resolves the -matrix flag, defaulting to DNA when both
// sequences look like DNA and BLOSUM62 otherwise.
func pickMatrix(name, a, b string) (*matrix.Matrix, error) {
	if name != "" && name != "auto" {
		m, ok := matrix.Get(name)
		if !ok {
			return nil, fmt.Errorf("unknown matrix %q (see 'seqalign matrices')", name)
		}
		return m, nil
	}
	if matrix.IsDNA(a) && matrix.IsDNA(b) {
		return matrix.DNA(1, -1), nil
	}
	m, _ := matrix.Get("blosum62")
	return m, nil
}

// defaultGaps returns the conventional gap penalties for a matrix:
// (10,1) for proteins (NCBI recommendation for BLOSUM62), (5,1) for DNA.
func defaultGaps(m *matrix.Matrix) (open, extend int) {
	if m.IsProtein {
		return 10, 1
	}
	return 5, 1
}

// alignParams bundles the flags shared by align/batch/score.
type alignParams struct {
	algo      align.Algorithm
	matrix    *matrix.Matrix
	gapOpen   int
	gapExtend int
}

// addAlignFlags registers the shared flags on fs.
func addAlignFlags(fs *flag.FlagSet) {
	fs.String("algo", "nw", "alignment algorithm: nw (global), sw (local), semi (overlap)")
	fs.String("matrix", "auto", "scoring matrix: auto, dna, blosum45, blosum62, blosum80, pam250")
	fs.Int("gap-open", 0, "gap opening penalty (default 10 protein / 5 DNA)")
	fs.Int("gap-extend", 0, "gap extension penalty per residue (default 1)")
	fs.Int("match", 1, "DNA match score (only with a DNA matrix)")
	fs.Int("mismatch", -1, "DNA mismatch score (only with a DNA matrix)")
}

type alignOutput struct {
	Algorithm string      `json:"algorithm"`
	Matrix    string      `json:"matrix"`
	GapOpen   int         `json:"gap_open"`
	GapExtend int         `json:"gap_extend"`
	Score     int         `json:"score"`
	Stats     stats.Stats `json:"stats"`
	AlignedA  string      `json:"aligned_a"`
	AlignedB  string      `json:"aligned_b"`
	IDa       string      `json:"id_a"`
	IDb       string      `json:"id_b"`
}

func cmdAlign(args []string) error {
	fs := flag.NewFlagSet("align", flag.ContinueOnError)
	addAlignFlags(fs)
	format := fs.String("format", "text", "output format: text or json")
	outPath := fs.String("out", "", "write output to this file instead of stdout")
	fs.Usage = func() {
		fmt.Fprint(fs.Output(), "Usage: seqalign align [flags] SEQ1 SEQ2\n\nFlags:\n")
		fs.PrintDefaults()
	}
	// Parse flags first, then resolve the sequences, then derive the
	// scoring parameters (the matrix auto-pick needs the sequences).
	if err := fs.Parse(args); err != nil {
		return err
	}
	rest := fs.Args()
	if len(rest) != 2 {
		return fmt.Errorf("align needs exactly two sequences (got %d)", len(rest))
	}
	seqA, idA, err := resolveSeq(rest[0], "seq1")
	if err != nil {
		return err
	}
	seqB, idB, err := resolveSeq(rest[1], "seq2")
	if err != nil {
		return err
	}
	p, err := deriveParams(fs, seqA, seqB)
	if err != nil {
		return err
	}
	res, err := align.Align(seqA, seqB, p.algo, align.Params{Matrix: p.matrix, GapOpen: p.gapOpen, GapExtend: p.gapExtend})
	if err != nil {
		return err
	}
	st := stats.Compute(res.AlignedA, res.AlignedB, p.matrix)
	st.Score = res.Score

	var w io.Writer = os.Stdout
	if *outPath != "" {
		f, err := os.Create(*outPath)
		if err != nil {
			return err
		}
		defer f.Close()
		w = f
	}
	if *format == "json" {
		enc := json.NewEncoder(w)
		enc.SetIndent("", "  ")
		return enc.Encode(alignOutput{
			Algorithm: p.algo.String(), Matrix: p.matrix.Name,
			GapOpen: p.gapOpen, GapExtend: p.gapExtend,
			Score: res.Score, Stats: st,
			AlignedA: res.AlignedA, AlignedB: res.AlignedB,
			IDa: idA, IDb: idB,
		})
	}
	if *format != "text" {
		return fmt.Errorf("unknown format %q (want text or json)", *format)
	}
	writeTextAlignment(w, idA, idB, p, res, st)
	return nil
}

// deriveParams reads the already-parsed flag set and resolves the matrix
// and gap penalties for the given sequences.
func deriveParams(fs *flag.FlagSet, a, b string) (alignParams, error) {
	var p alignParams
	algo, err := align.ParseAlgorithm(fs.Lookup("algo").Value.String())
	if err != nil {
		return p, err
	}
	p.algo = algo
	m, err := pickMatrix(fs.Lookup("matrix").Value.String(), a, b)
	if err != nil {
		return p, err
	}
	set := map[string]bool{}
	fs.Visit(func(f *flag.Flag) { set[f.Name] = true })
	if !m.IsProtein && (set["match"] || set["mismatch"]) {
		m = matrix.DNA(fs.Lookup("match").Value.(flag.Getter).Get().(int),
			fs.Lookup("mismatch").Value.(flag.Getter).Get().(int))
	}
	p.matrix = m
	defOpen, defExtend := defaultGaps(m)
	p.gapOpen, p.gapExtend = defOpen, defExtend
	if set["gap-open"] {
		p.gapOpen = fs.Lookup("gap-open").Value.(flag.Getter).Get().(int)
	}
	if set["gap-extend"] {
		p.gapExtend = fs.Lookup("gap-extend").Value.(flag.Getter).Get().(int)
	}
	return p, nil
}

// writeTextAlignment prints a classic pairwise alignment with a match line.
func writeTextAlignment(w io.Writer, idA, idB string, p alignParams, res align.Result, st stats.Stats) {
	fmt.Fprintf(w, "%s alignment (%s, gap open %d, extend %d)\n",
		p.algo.FullName(), p.matrix.Name, p.gapOpen, p.gapExtend)
	fmt.Fprintf(w, "Score: %d\n", res.Score)
	fmt.Fprintf(w, "Identity: %.2f%% (%d/%d)  Similarity: %.2f%% (%d/%d)  Gaps: %d chars in %d opening(s)  Length: %d\n\n",
		st.IdentityPct, st.Matches, st.Length, st.SimilarityPct, st.Similar, st.Length,
		st.GapChars, st.GapOpens, st.Length)
	const width = 60
	posA, posB := 0, 0
	for off := 0; off < len(res.AlignedA); off += width {
		end := off + width
		if end > len(res.AlignedA) {
			end = len(res.AlignedA)
		}
		segA, segB := res.AlignedA[off:end], res.AlignedB[off:end]
		var mid strings.Builder
		for k := 0; k < len(segA); k++ {
			switch {
			case segA[k] == '-' || segB[k] == '-':
				mid.WriteByte(' ')
			case segA[k] == segB[k]:
				mid.WriteByte('|')
			default:
				if v, ok := p.matrix.Score(segA[k], segB[k]); ok && v > 0 {
					mid.WriteByte(':')
				} else {
					mid.WriteByte(' ')
				}
			}
		}
		startA, startB := posA+1, posB+1
		for k := 0; k < len(segA); k++ {
			if segA[k] != '-' {
				posA++
			}
			if segB[k] != '-' {
				posB++
			}
		}
		endA, endB := posA, posB
		if len(segA) == 0 {
			endA, endB = posA, posB
		}
		fmt.Fprintf(w, "%-12s %4d  %s  %d\n", idA, startA, segA, endA)
		fmt.Fprintf(w, "%-12s %4s  %s\n", "", "", mid.String())
		fmt.Fprintf(w, "%-12s %4d  %s  %d\n\n", idB, startB, segB, endB)
	}
}

type batchRow struct {
	Pair          int     `json:"pair"`
	IDa           string  `json:"id_a"`
	IDb           string  `json:"id_b"`
	Algorithm     string  `json:"algorithm"`
	Matrix        string  `json:"matrix"`
	Score         int     `json:"score"`
	IdentityPct   float64 `json:"identity_pct"`
	SimilarityPct float64 `json:"similarity_pct"`
	Matches       int     `json:"matches"`
	Mismatches    int     `json:"mismatches"`
	GapChars      int     `json:"gap_chars"`
	GapOpens      int     `json:"gap_opens"`
	Length        int     `json:"length"`
	AlignedA      string  `json:"aligned_a,omitempty"`
	AlignedB      string  `json:"aligned_b,omitempty"`
}

func cmdBatch(args []string) error {
	fs := flag.NewFlagSet("batch", flag.ContinueOnError)
	addAlignFlags(fs)
	inPath := fs.String("in", "", "input FASTA file with consecutive record pairs (required)")
	outPath := fs.String("out", "", "output file (default: stdout)")
	format := fs.String("format", "tsv", "output format: tsv or json")
	alnPath := fs.String("alignments", "", "also write full alignments to this file")
	fs.Usage = func() {
		fmt.Fprint(fs.Output(), "Usage: seqalign batch -in pairs.fasta [-out result.tsv] [flags]\n\nFlags:\n")
		fs.PrintDefaults()
	}
	if err := fs.Parse(args); err != nil {
		return err
	}
	if *inPath == "" {
		return fmt.Errorf("batch needs -in pairs.fasta")
	}
	records, err := fasta.ParseFile(*inPath)
	if err != nil {
		return err
	}
	pairs, err := fasta.Pairs(records)
	if err != nil {
		return err
	}
	// Scoring parameters are resolved per pair: with -matrix auto each
	// pair gets its sensible default (DNA for DNA pairs, BLOSUM62 for
	// protein pairs); an explicit -matrix applies to every pair.
	rows := make([]batchRow, 0, len(pairs))
	var alnOut io.Writer
	var alnFile *os.File
	if *alnPath != "" {
		alnFile, err = os.Create(*alnPath)
		if err != nil {
			return err
		}
		defer alnFile.Close()
		alnOut = alnFile
	}
	for i, pr := range pairs {
		p, err := deriveParams(fs, pr[0].Sequence, pr[1].Sequence)
		if err != nil {
			return fmt.Errorf("pair %d (%s vs %s): %w", i+1, pr[0].ID, pr[1].ID, err)
		}
		res, err := align.Align(pr[0].Sequence, pr[1].Sequence, p.algo,
			align.Params{Matrix: p.matrix, GapOpen: p.gapOpen, GapExtend: p.gapExtend})
		if err != nil {
			return fmt.Errorf("pair %d (%s vs %s): %w", i+1, pr[0].ID, pr[1].ID, err)
		}
		st := stats.Compute(res.AlignedA, res.AlignedB, p.matrix)
		st.Score = res.Score
		rows = append(rows, batchRow{
			Pair: i + 1, IDa: pr[0].ID, IDb: pr[1].ID,
			Algorithm: p.algo.String(), Matrix: p.matrix.Name,
			Score: res.Score, IdentityPct: st.IdentityPct, SimilarityPct: st.SimilarityPct,
			Matches: st.Matches, Mismatches: st.Mismatches,
			GapChars: st.GapChars, GapOpens: st.GapOpens, Length: st.Length,
			AlignedA: res.AlignedA, AlignedB: res.AlignedB,
		})
		if alnOut != nil {
			fmt.Fprintf(alnOut, "=== Pair %d: %s vs %s ===\n", i+1, pr[0].ID, pr[1].ID)
			writeTextAlignment(alnOut, pr[0].ID, pr[1].ID, p, res, st)
		}
	}

	var w io.Writer = os.Stdout
	if *outPath != "" {
		f, err := os.Create(*outPath)
		if err != nil {
			return err
		}
		defer f.Close()
		w = f
	}
	switch *format {
	case "tsv":
		fmt.Fprintln(w, "pair\tid_a\tid_b\talgorithm\tmatrix\tscore\tidentity_pct\tsimilarity_pct\tmatches\tmismatches\tgap_chars\tgap_opens\tlength")
		for _, r := range rows {
			fmt.Fprintf(w, "%d\t%s\t%s\t%s\t%s\t%d\t%.2f\t%.2f\t%d\t%d\t%d\t%d\t%d\n",
				r.Pair, r.IDa, r.IDb, r.Algorithm, r.Matrix, r.Score,
				r.IdentityPct, r.SimilarityPct, r.Matches, r.Mismatches,
				r.GapChars, r.GapOpens, r.Length)
		}
	case "json":
		enc := json.NewEncoder(w)
		enc.SetIndent("", "  ")
		return enc.Encode(rows)
	default:
		return fmt.Errorf("unknown format %q (want tsv or json)", *format)
	}
	return nil
}

func cmdMatrices(args []string) error {
	fs := flag.NewFlagSet("matrices", flag.ContinueOnError)
	name := fs.String("name", "", "print this matrix in full (e.g. BLOSUM62)")
	fs.Usage = func() {
		fmt.Fprint(fs.Output(), "Usage: seqalign matrices [-name BLOSUM62]\n\nFlags:\n")
		fs.PrintDefaults()
	}
	if err := fs.Parse(args); err != nil {
		return err
	}
	if *name != "" {
		m, ok := matrix.Get(*name)
		if !ok {
			return fmt.Errorf("unknown matrix %q", *name)
		}
		printMatrix(os.Stdout, m)
		return nil
	}
	for _, m := range matrix.List() {
		fmt.Printf("%-12s %s\n", m.Name, m.Description)
	}
	return nil
}

func printMatrix(w io.Writer, m *matrix.Matrix) {
	fmt.Fprintf(w, "# %s — %s\n", m.Name, m.Description)
	fmt.Fprintf(w, "%4s", "")
	for i := 0; i < len(m.Alphabet); i++ {
		fmt.Fprintf(w, "%4c", m.Alphabet[i])
	}
	fmt.Fprintln(w)
	for i := 0; i < len(m.Alphabet); i++ {
		fmt.Fprintf(w, "%4c", m.Alphabet[i])
		for j := 0; j < len(m.Alphabet); j++ {
			fmt.Fprintf(w, "%4d", m.Scores[i][j])
		}
		fmt.Fprintln(w)
	}
}

func cmdScore(args []string) error {
	fs := flag.NewFlagSet("score", flag.ContinueOnError)
	addAlignFlags(fs)
	format := fs.String("format", "text", "output format: text or json")
	fs.Usage = func() {
		fmt.Fprint(fs.Output(), "Usage: seqalign score [flags] ALIGNED1 ALIGNED2\n\nFlags:\n")
		fs.PrintDefaults()
	}
	if err := fs.Parse(args); err != nil {
		return err
	}
	rest := fs.Args()
	if len(rest) != 2 {
		return fmt.Errorf("score needs exactly two pre-aligned sequences (got %d)", len(rest))
	}
	alnA, idA, err := resolveSeq(rest[0], "seq1")
	if err != nil {
		return err
	}
	alnB, idB, err := resolveSeq(rest[1], "seq2")
	if err != nil {
		return err
	}
	// Strip gaps for the matrix auto-pick, then keep the aligned strings.
	ungappedA := strings.ReplaceAll(alnA, "-", "")
	ungappedB := strings.ReplaceAll(alnB, "-", "")
	p, err := deriveParams(fs, ungappedA, ungappedB)
	if err != nil {
		return err
	}
	score, err := align.ScoreAlignment(alnA, alnB, p.matrix, p.gapOpen, p.gapExtend)
	if err != nil {
		return err
	}
	st := stats.Compute(alnA, alnB, p.matrix)
	st.Score = score
	if *format == "json" {
		enc := json.NewEncoder(os.Stdout)
		enc.SetIndent("", "  ")
		return enc.Encode(alignOutput{
			Algorithm: "scored", Matrix: p.matrix.Name,
			GapOpen: p.gapOpen, GapExtend: p.gapExtend,
			Score: score, Stats: st,
			AlignedA: alnA, AlignedB: alnB, IDa: idA, IDb: idB,
		})
	}
	if *format != "text" {
		return fmt.Errorf("unknown format %q (want text or json)", *format)
	}
	fmt.Printf("Matrix: %s (gap open %d, extend %d)\n", p.matrix.Name, p.gapOpen, p.gapExtend)
	fmt.Printf("Score: %d\n", score)
	fmt.Printf("Identity: %.2f%% (%d/%d)  Similarity: %.2f%% (%d/%d)  Gaps: %d chars in %d opening(s)  Length: %d\n",
		st.IdentityPct, st.Matches, st.Length, st.SimilarityPct, st.Similar, st.Length,
		st.GapChars, st.GapOpens, st.Length)
	return nil
}
