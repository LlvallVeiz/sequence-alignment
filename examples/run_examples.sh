#!/bin/sh
# Regenerate the example outputs in examples/output/ from real runs.
# Usage: sh examples/run_examples.sh   (run from the repository root)
set -e
cd "$(dirname "$0")/.."

go build -o /tmp/seqalign ./cmd/seqalign

OUT=examples/output
mkdir -p "$OUT"

# 1. DNA global + local alignment (Needleman-Wunsch / Smith-Waterman,
#    DNA +1/-1, gap open 5, extend 1).
/tmp/seqalign align GATTACA GCATGCT > "$OUT/dna_nw.txt"
/tmp/seqalign align -algo sw GATTACA GCATGCT > "$OUT/dna_sw.txt"

# 2. Hemoglobin alpha vs beta (BLOSUM62, gap 10/1), global + local.
#    The @file form reads the first record of a FASTA file.
/tmp/seqalign align -matrix blosum62 @examples/hba_human.fasta @examples/hbb_human.fasta \
  > "$OUT/protein_hba_hbb_blosum62.txt"
/tmp/seqalign align -algo sw -matrix blosum62 @examples/hba_human.fasta @examples/hbb_human.fasta \
  > "$OUT/protein_hba_hbb_sw.txt"

# 3. Human vs chimpanzee hemoglobin alpha (BLOSUM62, gap 10/1).
/tmp/seqalign align -matrix blosum62 \
  MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR \
  VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR \
  > "$OUT/protein_hba_human_chimp.txt"

# 4. Batch mode over the mixed DNA/protein pair file (auto matrix per pair).
/tmp/seqalign batch -in examples/batch_pairs.fasta -out "$OUT/batch.tsv" -alignments "$OUT/batch_alignments.txt"
/tmp/seqalign batch -in examples/batch_pairs.fasta -out "$OUT/batch.json" -format json

# 5. Matrix catalogue + scoring a given alignment.
/tmp/seqalign matrices > "$OUT/matrices.txt"
/tmp/seqalign score A--GT ACCGT > "$OUT/score.txt"

echo "Example outputs written to $OUT/"
