package main

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

var seqalignBin string

func TestMain(m *testing.M) {
	dir, err := os.MkdirTemp("", "seqalign-test-*")
	if err != nil {
		panic(err)
	}
	defer os.RemoveAll(dir)
	seqalignBin = filepath.Join(dir, "seqalign")
	build := exec.Command("go", "build", "-o", seqalignBin, ".")
	build.Dir = "."
	if out, err := build.CombinedOutput(); err != nil {
		panic("go build for tests failed: " + err.Error() + "\n" + string(out))
	}
	os.Exit(m.Run())
}

func runCLI(t *testing.T, args ...string) (string, string, int) {
	t.Helper()
	cmd := exec.Command(seqalignBin, args...)
	var stdout, stderr strings.Builder
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	err := cmd.Run()
	code := 0
	if err != nil {
		if ee, ok := err.(*exec.ExitError); ok {
			code = ee.ExitCode()
		} else {
			t.Fatalf("running seqalign %v: %v", args, err)
		}
	}
	return stdout.String(), stderr.String(), code
}

func TestCLIAlignDNA(t *testing.T) {
	out, _, code := runCLI(t, "align", "GATTACA", "GCATGCT")
	if code != 0 {
		t.Fatalf("exit %d, output:\n%s", code, out)
	}
	if !strings.Contains(out, "Score: -1") {
		t.Errorf("expected 'Score: -1' in output:\n%s", out)
	}
	if !strings.Contains(out, "GATTACA") || !strings.Contains(out, "GCATGCT") {
		t.Errorf("alignment lines missing from output:\n%s", out)
	}
}

func TestCLIAlignLocalJSON(t *testing.T) {
	out, _, code := runCLI(t, "align", "-algo", "sw", "-format", "json", "-matrix", "blosum62", "HEAGAWGHEE", "PAWHEAE")
	if code != 0 {
		t.Fatalf("exit %d, output:\n%s", code, out)
	}
	if !strings.Contains(out, `"score": 17`) {
		t.Errorf("expected JSON score 17 in output:\n%s", out)
	}
}

func TestCLIAlignFASTAFile(t *testing.T) {
	dir := t.TempDir()
	fa := filepath.Join(dir, "hba.fasta")
	if err := os.WriteFile(fa, []byte(">HBA_HUMAN\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	fb := filepath.Join(dir, "hbb.fasta")
	if err := os.WriteFile(fb, []byte(">HBB_HUMAN\nMVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	out, _, code := runCLI(t, "align", "-matrix", "blosum62", "@"+fa, "@"+fb)
	if code != 0 {
		t.Fatalf("exit %d, output:\n%s", code, out)
	}
	if !strings.Contains(out, "Score: 286") {
		t.Errorf("expected 'Score: 286' for HBA/HBB in output:\n%s", out)
	}
}

func TestCLIMatrices(t *testing.T) {
	out, _, code := runCLI(t, "matrices")
	if code != 0 {
		t.Fatalf("exit %d", code)
	}
	for _, name := range []string{"BLOSUM45", "BLOSUM62", "BLOSUM80", "PAM250"} {
		if !strings.Contains(out, name) {
			t.Errorf("matrices list missing %s:\n%s", name, out)
		}
	}
	out, _, code = runCLI(t, "matrices", "-name", "BLOSUM62")
	if code != 0 {
		t.Fatalf("exit %d", code)
	}
	if !strings.Contains(out, "W") || !strings.Contains(out, "11") {
		t.Errorf("BLOSUM62 printout looks wrong:\n%s", out)
	}
}

func TestCLIScore(t *testing.T) {
	out, _, code := runCLI(t, "score", "A--GT", "ACCGT")
	if code != 0 {
		t.Fatalf("exit %d, output:\n%s", code, out)
	}
	if !strings.Contains(out, "Score: -4") {
		t.Errorf("expected 'Score: -4' in output:\n%s", out)
	}
}

func TestCLIBatch(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "pairs.fasta")
	pairs := ">dna1\nGATTACA\n>dna2\nGCATGCT\n>prot1\nHEAGAWGHEE\n>prot2\nPAWHEAE\n"
	if err := os.WriteFile(in, []byte(pairs), 0o644); err != nil {
		t.Fatal(err)
	}
	outPath := filepath.Join(dir, "out.tsv")
	_, _, code := runCLI(t, "batch", "-in", in, "-out", outPath, "-matrix", "blosum62", "-algo", "sw")
	if code != 0 {
		t.Fatalf("exit %d", code)
	}
	data, err := os.ReadFile(outPath)
	if err != nil {
		t.Fatal(err)
	}
	lines := strings.Split(strings.TrimSpace(string(data)), "\n")
	if len(lines) != 3 { // header + 2 pairs
		t.Fatalf("expected 3 TSV lines, got %d:\n%s", len(lines), string(data))
	}
	if !strings.HasPrefix(lines[0], "pair\tid_a\tid_b") {
		t.Errorf("bad TSV header: %q", lines[0])
	}
	if !strings.Contains(lines[1], "dna1") || !strings.Contains(lines[2], "prot1") {
		t.Errorf("TSV rows missing pair ids:\n%s", string(data))
	}
}

func TestCLIErrors(t *testing.T) {
	if _, _, code := runCLI(t, "align", "ACGT"); code == 0 {
		t.Error("align with one sequence should fail")
	}
	if _, _, code := runCLI(t, "align", "ACGT", "ACGT", "-matrix", "BLOSUM99"); code == 0 {
		t.Error("unknown matrix should fail")
	}
	if _, _, code := runCLI(t, "align", "-matrix", "dna", "ACGT", "ACGTX"); code == 0 {
		t.Error("invalid DNA character should fail")
	}
	if _, _, code := runCLI(t, "nonsense"); code == 0 {
		t.Error("unknown subcommand should fail")
	}
}
