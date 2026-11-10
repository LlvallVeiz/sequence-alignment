#!/usr/bin/env python3
"""Regenerate internal/matrix/data.go from the NCBI matrix files.

The source of truth is the NCBI BLAST matrix files as shipped by
Biopython (Bio/Align/substitution_matrices/data/{BLOSUM45,BLOSUM62,
BLOSUM80,PAM250}). The script only re-emits the tables; the provenance
and verification notes in the header are kept here so a regeneration
does not silently drop them.

Requires: biopython  (pip install biopython)
Usage:    python3 scripts/generate_matrix_data.py
"""

from pathlib import Path

from Bio.Align import substitution_matrices

ORDER24 = "ARNDCQEGHILKMFPSTWYVBZX*"
NAMES = ["BLOSUM45", "BLOSUM62", "BLOSUM80", "PAM250"]

HEADER = '''package matrix

// Data tables below are the canonical NCBI-distributed substitution
// matrices, reproduced verbatim (integer values, 24-letter alphabet
// A R N D C Q E G H I L K M F P S T W Y V B Z X *).
//
// Provenance and verification (2026-11-10): values were taken from the
// NCBI BLAST matrix files as shipped by Biopython
// (Bio/Align/substitution_matrices/data/{BLOSUM45,BLOSUM62,BLOSUM80,PAM250})
// and cross-checked entry-by-entry against the `scoring-matrices` Python
// package (NCBI BLAST+ data): 0 differences on the 20x20 core of all four
// matrices. Spot checks against independent web sources:
//   - BLOSUM45: 0/576 differences vs CalState LA blosum45.txt
//     (https://www.calstatela.edu/sites/default/files/blosum45.txt)
//   - BLOSUM62: diagonal matches the CMU/DeepWiki listing exactly
//     (A=4 R=5 N=6 D=6 C=9 Q=5 E=5 G=6 H=8 I=4 L=4 K=5 M=5 F=6
//     P=7 S=4 T=5 W=11 Y=7 V=4); the CMU course page's T column
//     (T-T=4, P-T=1, ...) disagrees with the NCBI file (T-T=5,
//     P-T=-1, ...) and was treated as a transcription error there.
//   - PAM250: diagonal (A=2 R=6 N=2 D=4 C=12 Q=4 E=4 G=5 H=6 I=5
//     L=6 K=5 M=6 F=9 P=6 S=2 T=3 W=17 Y=10 V=4) matches the
//     standard Dayhoff table; the GPM/Pearson 1990 listing differs in
//     7 entries (C-C=4, S-S=3, F-A=-4, F-N=-4, P-N=-1, P-G=-1,
//     T-F=-2) and was NOT used.
//   - BLOSUM80: the `blosum` PyPI package ships a different BLOSUM80
//     variant (A-A=5, W-W=11); this project uses the NCBI file
//     (1/3 bit units, A-A=7, W-W=16, entropy 0.9868).
// Regenerate with: python3 scripts/generate_matrix_data.py

// proteinAlphabet is the NCBI ordering shared by all built-in protein matrices.
const proteinAlphabet = "ARNDCQEGHILKMFPSTWYVBZX*"
'''


def main() -> None:
    out = [HEADER]
    for name in NAMES:
        m = substitution_matrices.load(name)
        var = name.lower() + "Data"
        out.append(f"\n// {var} holds {name} in proteinAlphabet order.\n")
        out.append(f"var {var} = [24][24]int{{\n")
        for a in ORDER24:
            row = ", ".join(f"{int(m[a, b]):3d}" for b in ORDER24)
            out.append(f"\t{{{row}}}, // {a}\n")
        out.append("}\n")
    root = Path(__file__).resolve().parent.parent
    target = root / "internal" / "matrix" / "data.go"
    target.write_text("".join(out))
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
