// Package fasta parses FASTA files: one or more records, each a '>'
// header line followed by one or more sequence lines.
package fasta

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"
)

// Record is one FASTA entry.
type Record struct {
	// ID is the first whitespace-separated token of the header.
	ID string
	// Description is the full header text after '>'.
	Description string
	// Sequence is the concatenated sequence lines, upper-cased,
	// whitespace removed.
	Sequence string
}

// Parse reads all records from r. Blank lines are ignored; sequence lines
// before the first header are an error.
func Parse(r io.Reader) ([]Record, error) {
	var records []Record
	var cur *Record
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 64*1024), 4*1024*1024)
	lineNo := 0
	for sc.Scan() {
		lineNo++
		line := strings.TrimSpace(sc.Text())
		if line == "" {
			continue
		}
		if strings.HasPrefix(line, ">") {
			if cur != nil {
				records = append(records, *cur)
			}
			header := strings.TrimSpace(line[1:])
			if header == "" {
				return nil, fmt.Errorf("fasta: empty header at line %d", lineNo)
			}
			id := header
			if i := strings.IndexAny(header, " \t"); i >= 0 {
				id = header[:i]
			}
			cur = &Record{ID: id, Description: header}
			continue
		}
		if cur == nil {
			return nil, fmt.Errorf("fasta: sequence data before first header at line %d", lineNo)
		}
		cur.Sequence += strings.ToUpper(strings.Join(strings.Fields(line), ""))
	}
	if err := sc.Err(); err != nil {
		return nil, fmt.Errorf("fasta: reading input: %w", err)
	}
	if cur != nil {
		records = append(records, *cur)
	}
	if len(records) == 0 {
		return nil, fmt.Errorf("fasta: no records found")
	}
	return records, nil
}

// ParseFile reads all records from the file at path.
func ParseFile(path string) ([]Record, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("fasta: %w", err)
	}
	defer f.Close()
	return Parse(f)
}

// ParseString parses FASTA text held in a string.
func ParseString(s string) ([]Record, error) {
	return Parse(strings.NewReader(s))
}

// First returns the first record of a FASTA file, or an error.
func First(path string) (Record, error) {
	recs, err := ParseFile(path)
	if err != nil {
		return Record{}, err
	}
	return recs[0], nil
}

// Pairs validates that records come in consecutive pairs (1-2, 3-4, ...)
// for batch alignment and returns them as pairs.
func Pairs(records []Record) ([][2]Record, error) {
	if len(records)%2 != 0 {
		return nil, fmt.Errorf("fasta: batch input needs an even number of records (pairs), got %d", len(records))
	}
	pairs := make([][2]Record, 0, len(records)/2)
	for i := 0; i < len(records); i += 2 {
		pairs = append(pairs, [2]Record{records[i], records[i+1]})
	}
	return pairs, nil
}
