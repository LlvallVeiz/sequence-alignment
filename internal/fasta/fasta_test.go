package fasta

import (
	"strings"
	"testing"
)

const multi = `>seq1 first record
ACGTACGT
ACGT
>seq2 second record | with description
acgt
>seq3
NNNN
`

func TestParseMultiRecord(t *testing.T) {
	recs, err := ParseString(multi)
	if err != nil {
		t.Fatal(err)
	}
	if len(recs) != 3 {
		t.Fatalf("got %d records, want 3", len(recs))
	}
	if recs[0].ID != "seq1" || recs[0].Sequence != "ACGTACGTACGT" {
		t.Errorf("record 1 = %+v", recs[0])
	}
	if recs[1].ID != "seq2" || recs[1].Description != "seq2 second record | with description" {
		t.Errorf("record 2 = %+v", recs[1])
	}
	if recs[1].Sequence != "ACGT" {
		t.Errorf("record 2 sequence = %q, want ACGT (lowercase input upper-cased)", recs[1].Sequence)
	}
	if recs[2].Sequence != "NNNN" {
		t.Errorf("record 3 sequence = %q, want NNNN", recs[2].Sequence)
	}
}

func TestParseErrors(t *testing.T) {
	for _, bad := range []string{
		"",
		"ACGT\n",         // sequence before any header
		">\nACGT\n",      // empty header
		">only-header\n", // header with no sequence is fine actually
	} {
		_, err := ParseString(bad)
		if bad == ">only-header\n" {
			// A header with an empty sequence is a valid (if useless) record.
			if err != nil {
				t.Errorf("header-only record should parse, got %v", err)
			}
			continue
		}
		if err == nil {
			t.Errorf("input %q should fail to parse", bad)
		}
	}
}

func TestParseCRLFAndBlanks(t *testing.T) {
	recs, err := ParseString(">a\r\nACGT\r\n\r\n> b\r\nTGCA\r\n")
	if err != nil {
		t.Fatal(err)
	}
	if len(recs) != 2 || recs[0].Sequence != "ACGT" || recs[1].Sequence != "TGCA" {
		t.Errorf("CRLF parse = %+v", recs)
	}
}

func TestPairs(t *testing.T) {
	recs, _ := ParseString(multi)
	pairs, err := Pairs(recs)
	if err == nil {
		t.Errorf("odd record count should fail, got %v", pairs)
	}
	recs4, _ := ParseString(multi + ">seq4\nACGT\n")
	pairs, err = Pairs(recs4)
	if err != nil {
		t.Fatal(err)
	}
	if len(pairs) != 2 || pairs[0][0].ID != "seq1" || pairs[0][1].ID != "seq2" {
		t.Errorf("pairs = %+v", pairs)
	}
}

func TestParseFileMissing(t *testing.T) {
	if _, err := ParseFile("/nonexistent/nope.fasta"); err == nil {
		t.Error("missing file should be an error")
	}
}

func TestRoundTripDescription(t *testing.T) {
	recs, err := ParseString(">HBA_HUMAN Hemoglobin subunit alpha\nMVLSPADKTNVK\n")
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(recs[0].Description, "Hemoglobin") {
		t.Errorf("description lost: %q", recs[0].Description)
	}
}
