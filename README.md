# sequence-alignment

Pairwise sequence alignment di Go: **Needleman–Wunsch** (global),
**Smith–Waterman** (local), dan **semi-global/overlap**, semuanya dengan
gap affine (Gotoh 1982). Matriks skor protein BLOSUM45, BLOSUM62,
BLOSUM80, dan PAM250 memakai tabel NCBI yang terverifikasi silang —
bukan hafalan (lihat [`docs/RESEARCH.md`](docs/RESEARCH.md) untuk
bukti verifikasinya). CLI `seqalign` bisa menyejajarkan dua sekuens,
memproses banyak pasangan dari satu file FASTA, dan menilai alignment
yang sudah jadi.

## Teori singkat

Alignment optimal dicari dengan pemrograman dinamis atas tiga matriks
Gotoh — M (berakhir dengan pasangan residu), X (berakhir dengan gap di
sekuens A), Y (berakhir dengan gap di sekuens B) — sehingga gap
sepanjang k berharga `open + k·extend` (bukan `k·d` linear). Global
menyejajarkan seluruh panjang kedua sekuens; lokal mencari blok
substring terbaik (sel dilantai di nol, traceback dari sel maksimum);
semi-global membebaskan gap di ujung kedua sekuens (overlap/fitting).
Skor pasangan residu diambil dari matriks substitusi: BLOSUM dari
frekuensi substitusi yang diamati di blok konservasi (Henikoff &
Henikoff 1992), PAM dari model evolusi Markov (Dayhoff dkk. 1978).

## Instalasi

Butuh Go 1.22+.

```sh
git clone https://github.com/LlvallVeiz/sequence-alignment.git
cd sequence-alignment
go build -o seqalign ./cmd/seqalign
```

Atau cukup `go run ./cmd/seqalign ...` tanpa build.

## Penggunaan

```
seqalign align [flags] SEQ1 SEQ2       sejajarkan dua sekuens
seqalign batch -in pairs.fasta ...      sejajarkan pasangan FASTA berurutan
seqalign matrices [-name BLOSUM62]      daftar (dan cetak) matriks skor
seqalign score [flags] ALN1 ALN2        nilai dua sekuens yang sudah disejajarkan
```

Argumen sekuens berupa sekuens mentah (`GATTACA`, `MVLSPAD...`) atau
`@file.fasta` (membaca record pertama file FASTA). Flag diletakkan
**sebelum** argumen sekuens.

| Flag | Default | Keterangan |
|---|---|---|
| `-algo` | `nw` | `nw` (global), `sw` (local), `semi` (overlap) |
| `-matrix` | `auto` | `auto`, `dna`, `blosum45`, `blosum62`, `blosum80`, `pam250` |
| `-gap-open` | 10 protein / 5 DNA | penalti buka gap |
| `-gap-extend` | 1 | penalti per residu gap |
| `-match` / `-mismatch` | +1 / −1 | skor DNA (hanya bila matriks DNA) |
| `-format` | `text` | `text` atau `json` (`align`, `score`); `tsv` atau `json` (`batch`) |
| `-out` | stdout | tulis ke file |
| `-alignments` | — | (`batch`) tulis juga alignment lengkap ke file |

Matriks `auto`: DNA bila kedua sekuens hanya berisi ACGTN, BLOSUM62
bila bukan. Input RNA (U) dibaca sebagai T untuk matriks DNA.

## Contoh output (hasil run beneran)

DNA, global (`seqalign align GATTACA GCATGCT`):

```
Needleman-Wunsch (global) alignment (DNA(+1/-1), gap open 5, extend 1)
Score: -1
Identity: 42.86% (3/7)  Similarity: 42.86% (3/7)  Gaps: 0 chars in 0 opening(s)  Length: 7

GATTACA         1  GATTACA  7
                   |  | |
GCATGCT         1  GCATGCT  7
```

DNA, lokal (`seqalign align -algo sw GATTACA GCATGCT`): skor **2**,
blok terbaik `AT`/`AT` (identitas 100%, panjang 2).

Hemoglobin manusia α vs β, global (BLOSUM62, gap 10/1) — sekuens asli
UniProt P69905 (142 aa) vs P68871 (147 aa):

```
Needleman-Wunsch (global) alignment (BLOSUM62, gap open 10, extend 1)
Score: 286
Identity: 43.62% (65/149)  Similarity: 60.40% (90/149)  Gaps: 9 chars in 4 opening(s)  Length: 149

HBA_HUMAN       1  MV-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF-DLS-----HGS  53
                   || |:| :|: | | ||||  :  | | ||| |: : :| |: :|  | |||      |:
HBB_HUMAN       1  MVHLTPEEKSAVTALWGKV--NVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGN  58
```

Mode lokal untuk pasangan yang sama: skor **288**, identitas 43,45%
(63/145) — lokal sedikit lebih tinggi karena membuang ujung N-terminal
yang buruk (MVL vs MVH).

Hemoglobin α manusia vs simpanse (P69907 residu 2–142): skor **717**,
identitas **99,30% (141/142)** — hanya Met awal yang tidak ada pada
catatan rekombinan simpanse.

Batch (`seqalign batch -in examples/batch_pairs.fasta -out out.tsv`):

```
pair  id_a         id_b         algorithm  matrix       score  identity_pct  ...
1     dna_pair1_a  dna_pair1_b  global     DNA(+1/-1)   -1     42.86
2     dna_pair2_a  dna_pair2_b  global     DNA(+1/-1)   -8     0.00
3     prot_pair1_a prot_pair1_b global     BLOSUM62     2      30.00
4     prot_pair2_a prot_pair2_b global     BLOSUM62     61     51.52
```

Pasangan 2 menunjukkan perbedaan global vs lokal: secara global
`ACGTACGT` vs `TTACGTAA` bernilai −8 (0 match tanpa gap), sementara
mode lokal menemukan blok `ACGTA` bernilai 5. Semua output contoh di
[`examples/output/`](examples/output/) dihasilkan ulang oleh
`sh examples/run_examples.sh`.

## Struktur proyek

```
sequence-alignment/
├── cmd/seqalign/            CLI: align, batch, matrices, score
├── internal/
│   ├── matrix/              matriks BLOSUM/PAM (data NCBI) + matriks DNA
│   │   └── data.go          tabel 24×24, dibuat dari file NCBI (lihat scripts/)
│   ├── align/               inti DP Gotoh: global, lokal, semi-global
│   ├── fasta/               parser FASTA multi-record
│   └── stats/               identitas, similaritas, gap
├── docs/RESEARCH.md         riset, verifikasi matriks & algoritma, referensi
├── examples/                FASTA contoh + output hasil run beneran
│   ├── dna_pair.fasta
│   ├── protein_hba_hbb.fasta
│   ├── protein_hba_human_chimp.fasta
│   ├── batch_pairs.fasta
│   └── output/
├── scripts/
│   ├── reference_align.py   implementasi referensi Python (cross-check)
│   └── generate_matrix_data.py  regenerasi data.go dari file NCBI
├── Makefile
└── LICENSE (MIT)
```

## Verifikasi

```sh
go test ./...        # 35 tes, semua lolos
go vet ./...
```

Yang diperiksa tes:

* **Nilai matriks**: entri kunci keempat matriks (mis. BLOSUM62 W–W=11,
  T–T=5; BLOSUM80 A–A=7, W–W=16; PAM250 C–C=12, W–W=17), simetri, dan
  kelengkapan tabel 24×24.
* **Hitung tangan**: skor kecil yang diverifikasi manual di komentar
  tes (mis. NW `AC` vs `A` = 1−(5+1) = −5; `ACDE` vs `ACDE` BLOSUM62
  = 4+9+6+5 = 24; SW `HEAGAWGHEE` vs `PAWHEAE` = 17).
* **Cross-check independen**: implementasi referensi Python
  (`scripts/reference_align.py`) dibandingkan pada 15 kasus acak
  (seed 20261110) untuk ketiga algoritma — semua skor cocok; pasangan
  hemoglobin α/β juga dicocokkan (NW 286, SW 288).
* **CLI end-to-end**: `align`/`batch`/`matrices`/`score` dijalankan
  sebagai binary dan outputnya diperiksa.

Verifikasi matriks menemukan bahwa beberapa sumber internet memakai
varian matriks yang berbeda (BLOSUM80 paket PyPI `blosum` berbeda 324
entri; daftar PAM250 GPM berbeda 7 entri; kolom T halaman CMU untuk
BLOSUM62 salah salin). Proyek ini memakai file NCBI dan mencatat
perbedaannya di `docs/RESEARCH.md` bagian 4.

## Keputusan desain

* **Gap affine Gotoh** (tiga matriks + pointer traceback), bukan gap
  linear: sesuai praktik BLAST/NCBI, tetap O(mn).
* **Konvensi gap `open + k·extend`** = konvensi NCBI `(a,b)`:
  gap panjang 1 dengan (10,1) berharga −11.
* **Default mengikuti NCBI**: protein (10,1), DNA (5,1); matriks
  auto (DNA bila ACGTN, BLOSUM62 bila bukan).
* **Matriks 24 huruf** (termasuk B/Z/X/\*) seperti file NCBI, supaya
  sekuens nyata yang mengandung X bisa disejajarkan.
* **Tie-breaking deterministik** (substitusi didahulukan atas gap;
  perpanjang gap didahulukan atas buka gap baru); tes mengunci skor,
  bukan bentuk alignment saat ada seri.
* **Inti DP terpisah dari CLI** (`internal/align` tidak bergantung
  pada flag atau I/O), sehingga bisa dipakai sebagai library.

## Keterbatasan

* Hanya pairwise alignment — bukan multiple alignment dan bukan
  database search (tidak ada heuristik seed-and-extend, tidak ada
  E-value/Karlin-Altschul).
* DP penuh O(mn) tanpa optimasi memori (tanpa Hirschberg): dua
  sekuens 10 kb butuh ~100 juta sel — lambat dan boros memori.
* Semi-global di sini varian *overlap* (gap ujung gratis di kedua
  sekuens); varian *fitting* satu arah tidak disediakan.
* Similaritas memakai ambang sederhana (skor substitusi > 0).
* U (selenosistein) dan O (pirolisin) ditolak pada matriks protein.

## Lisensi

MIT — lihat [LICENSE](LICENSE).
