#!/usr/bin/env python3
"""Independent reference implementation of the seqalign DP algorithms.

Purpose: cross-check the Go implementation (internal/align) with a
separately written implementation in another language. It deliberately
uses plain Python lists and the same conventions as the Go code:

  * gap of length k costs open + k*extend (both penalties positive),
  * three Gotoh matrices M (pair), X (gap in A), Y (gap in B),
  * global: Needleman-Wunsch; local: Smith-Waterman with a zero floor;
    semi-global: free end gaps on both sequences.

It is NOT used by the Go program; it only produces expected values for
the Go tests (see internal/align/align_test.go) and for the examples.

Usage:
    python3 scripts/reference_align.py NW  ACGT ACGT --match 1 --mismatch -1
    python3 scripts/reference_align.py SW  HEAGAWGHEE PAWHEAE --matrix BLOSUM62
    python3 scripts/reference_align.py SEMI ACGT TTACGTTT
    python3 scripts/reference_align.py --selftest   # runs the hand-checked cases
"""

import argparse
import sys

NEG = float("-inf")

# Small built-in matrices for the reference (DNA + BLOSUM62 diagonal spot
# values are NOT re-typed here: protein matrices are loaded from the
# Biopython data files when available, falling back to the Go data file).
PROTEIN_ORDER = "ARNDCQEGHILKMFPSTWYVBZX*"


def load_protein_matrix(name):
    """Load a protein matrix from Biopython's NCBI data files if present."""
    try:
        from Bio.Align import substitution_matrices
        m = substitution_matrices.load(name)
        return {(a, b): int(m[a, b]) for a in PROTEIN_ORDER for b in PROTEIN_ORDER}
    except Exception:
        pass
    # Fallback: parse internal/matrix/data.go (same NCBI tables).
    import re
    from pathlib import Path
    text = (Path(__file__).parent.parent / "internal" / "matrix" / "data.go").read_text()
    var = name.lower() + "Data"
    block = re.search(re.escape(var) + r" = \[24\]\[24\]int\{(.*?)\n\}", text, re.S).group(1)
    rows = re.findall(r"\{([^}]*)\}", block)
    mat = {}
    for a, row in zip(PROTEIN_ORDER, rows):
        vals = [int(x) for x in row.split(",") if x.strip()]
        for b, v in zip(PROTEIN_ORDER, vals):
            mat[(a, b)] = v
    return mat


def dna_matrix(match=1, mismatch=-1):
    order = "ACGTN"
    mat = {}
    for a in order:
        for b in order:
            if a == "N" or b == "N":
                mat[(a, b)] = 0
            elif a == b:
                mat[(a, b)] = match
            else:
                mat[(a, b)] = mismatch
    return mat


def align(a, b, algo, mat, gap_open, gap_extend):
    """Return (score, aligned_a, aligned_b)."""
    m, n = len(a), len(b)
    if m == 0 and n == 0:
        return 0, "", ""
    if m == 0 or n == 0:
        if algo == "SW":
            return 0, "", ""
        L = m + n
        if m == 0:
            return -(gap_open + L * gap_extend), "-" * n, b
        return -(gap_open + L * gap_extend), a, "-" * m

    def new_mat():
        return [[0.0] * (n + 1) for _ in range(m + 1)]

    M, X, Y = new_mat(), new_mat(), new_mat()
    tM = [[None] * (n + 1) for _ in range(m + 1)]
    tX = [[None] * (n + 1) for _ in range(m + 1)]
    tY = [[None] * (n + 1) for _ in range(m + 1)]
    d, e = gap_open, gap_extend

    if algo == "SEMI":
        for i in range(1, m + 1):
            M[i][0] = 0
            Y[i][0] = 0
            X[i][0] = NEG
            tM[i][0] = "M"
            tY[i][0] = "M" if i == 1 else "Y"
        for j in range(1, n + 1):
            M[0][j] = 0
            X[0][j] = 0
            Y[0][j] = NEG
            tM[0][j] = "M"
            tX[0][j] = "M" if j == 1 else "X"
    elif algo == "SW":
        pass  # all boundary cells stay 0
    else:  # NW
        for i in range(1, m + 1):
            M[i][0] = -(d + i * e)
            Y[i][0] = -(d + i * e)
            X[i][0] = NEG
            tM[i][0] = "M"
            tY[i][0] = "M" if i == 1 else "Y"
        for j in range(1, n + 1):
            M[0][j] = -(d + j * e)
            X[0][j] = -(d + j * e)
            Y[0][j] = NEG
            tM[0][j] = "M"
            tX[0][j] = "M" if j == 1 else "X"

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = mat[(a[i - 1], b[j - 1])]
            best = max((M[i - 1][j - 1], "M"), (X[i - 1][j - 1], "X"),
                       (Y[i - 1][j - 1], "Y"), key=lambda t: t[0])
            # tie-break M > X > Y: max() keeps the first maximum in order.
            val = best[0] + s
            if algo == "SW" and val < 0:
                val, best = 0, (0, None)
            M[i][j], tM[i][j] = val, best[1]
            o, x = M[i][j - 1] - (d + e), X[i][j - 1] - e
            X[i][j], tX[i][j] = (x, "X") if x >= o else (o, "M")
            if algo == "SW" and X[i][j] < 0:
                X[i][j], tX[i][j] = 0, None
            o, y = M[i - 1][j] - (d + e), Y[i - 1][j] - e
            Y[i][j], tY[i][j] = (y, "Y") if y >= o else (o, "M")
            if algo == "SW" and Y[i][j] < 0:
                Y[i][j], tY[i][j] = 0, None

    def traceback(i, j, state):
        ra, rb = [], []
        while i > 0 or j > 0:
            if state == "M":
                if algo == "SW" and M[i][j] == 0:
                    break
                if i == 0:
                    ra.append("-")
                    rb.append(b[j - 1])
                    j -= 1
                    continue
                if j == 0:
                    ra.append(a[i - 1])
                    rb.append("-")
                    i -= 1
                    continue
                ra.append(a[i - 1])
                rb.append(b[j - 1])
                state = tM[i][j]
                i -= 1
                j -= 1
            elif state == "X":
                if algo == "SW" and X[i][j] == 0:
                    break
                ra.append("-")
                rb.append(b[j - 1])
                state = tX[i][j]
                j -= 1
            else:
                if algo == "SW" and Y[i][j] == 0:
                    break
                ra.append(a[i - 1])
                rb.append("-")
                state = tY[i][j]
                i -= 1
        return "".join(reversed(ra)), "".join(reversed(rb))

    if algo == "NW":
        best = max((M[m][n], "M"), (X[m][n], "X"), (Y[m][n], "Y"),
                   key=lambda t: t[0])
        aa, bb = traceback(m, n, best[1])
        return int(best[0]), aa, bb
    if algo == "SW":
        best, bi, bj, bs = 0, 0, 0, None
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                for v, s in ((M[i][j], "M"), (X[i][j], "X"), (Y[i][j], "Y")):
                    if v > best:
                        best, bi, bj, bs = v, i, j, s
        if best == 0:
            return 0, "", ""
        aa, bb = traceback(bi, bj, bs)
        return int(best), aa, bb
    # SEMI
    best, bi, bj = M[m][0], m, 0
    for j in range(1, n + 1):
        if M[m][j] > best:
            best, bi, bj = M[m][j], m, j
    for i in range(m + 1):
        if M[i][n] > best:
            best, bi, bj = M[i][n], i, n
    aa, bb = traceback(bi, bj, "M")
    if bj < n:
        aa += "-" * (n - bj)
        bb += b[bj:]
    if bi < m:
        aa += a[bi:]
        bb += "-" * (m - bi)
    return int(best), aa, bb


def selftest():
    dna = dna_matrix(1, -1)
    cases = [
        # (algo, a, b, matrix, open, extend, expected_score, note)
        ("NW", "ACGT", "ACGT", dna, 5, 1, 4, "identical DNA"),
        ("NW", "A", "C", dna, 5, 1, -1, "single mismatch"),
        ("NW", "AC", "A", dna, 5, 1, -5, "1 match + gap of 1: 1-(5+1)"),
        ("NW", "ACGT", "ACG", dna, 5, 1, -3, "3 matches + gap of 1: 3-(5+1)"),
        ("SW", "ACGTACGT", "TTACGTAA", dna, 5, 1, 5, "local ACGTA (5 matches)"),
        ("SW", "AAAA", "CCCC", dna, 5, 1, 0, "no positive local alignment"),
        ("SEMI", "ACGT", "TTACGTTT", dna, 5, 1, 4, "fitting, free end gaps"),
        ("SEMI", "ACGT", "ACGT", dna, 5, 1, 4, "identical, semi-global"),
        ("NW", "", "", dna, 5, 1, 0, "both empty"),
        ("NW", "", "ACGT", dna, 5, 1, -9, "empty vs 4 bases: -(5+4)"),
    ]
    b62 = load_protein_matrix("BLOSUM62")
    cases += [
        ("NW", "ACDE", "ACDE", b62, 10, 1, 24, "diagonal 4+9+6+5"),
        ("NW", "A", "W", b62, 10, 1, -3, "single BLOSUM62 pair"),
        ("NW", "W", "W", b62, 10, 1, 11, "Trp diagonal"),
        ("NW", "HEAGAWGHEE", "PAWHEAE", b62, 10, 1, 2, "classic example, NW mode"),
        ("SW", "HEAGAWGHEE", "PAWHEAE", b62, 10, 1, 17, "classic SW example: HEA/HEA = 8+5+4"),
        ("NW", "GATTACA", "GCATGCT", dna, 5, 1, None, "textbook-style DNA pair"),
        ("SW", "GATTACA", "GCATGCT", dna, 5, 1, None, "textbook-style DNA pair, local"),
        ("SEMI", "GATTACA", "TTGCATGCTAA", dna, 5, 1, None, "DNA fitting"),
    ]
    failed = 0
    for algo, a, b, mat, op, ex, want, note in cases:
        score, aa, bb = align(a, b, algo, mat, op, ex)
        status = "ok" if want is None or score == want else "FAIL"
        if status == "FAIL":
            failed += 1
        print(f"[{status}] {algo} {a!r} vs {b!r} -> score {score} (want {want}) [{note}]")
        print(f"       {aa}\n       {bb}")
    return 1 if failed else 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("algo", nargs="?", choices=["NW", "SW", "SEMI"],
                    help="algorithm")
    ap.add_argument("a", nargs="?", help="first sequence")
    ap.add_argument("b", nargs="?", help="second sequence")
    ap.add_argument("--matrix", default="DNA", help="DNA or BLOSUM45/62/80, PAM250")
    ap.add_argument("--match", type=int, default=1)
    ap.add_argument("--mismatch", type=int, default=-1)
    ap.add_argument("--gap-open", type=int, default=5)
    ap.add_argument("--gap-extend", type=int, default=1)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(argv)
    if args.selftest:
        return selftest()
    if not args.algo or not args.a or not args.b:
        ap.error("need ALGO A B (or --selftest)")
    if args.matrix.upper() == "DNA":
        mat = dna_matrix(args.match, args.mismatch)
        op, ex = args.gap_open, args.gap_extend
        # Accept RNA input: U is read as T.
        args.a = args.a.upper().replace("U", "T")
        args.b = args.b.upper().replace("U", "T")
    else:
        mat = load_protein_matrix(args.matrix.upper())
        op, ex = (args.gap_open, args.gap_extend)
        if args.gap_open == 5 and args.gap_extend == 1:
            op, ex = 10, 1  # protein default, mirroring the Go CLI
    score, aa, bb = align(args.a.upper(), args.b.upper(), args.algo, mat, op, ex)
    print(f"score: {score}")
    print(aa)
    print(bb)
    return 0


if __name__ == "__main__":
    sys.exit(main())
