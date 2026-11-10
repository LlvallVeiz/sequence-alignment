# Research Notes — sequence-alignment

Riset latar belakang untuk `seqalign`: algoritma pairwise alignment,
matriks skor substitusi, dan pemilihan parameter. Semua klaim di bawah
ini merujuk ke sumber yang tercantum; nilai matriks yang dipakai program
bukan hasil hafalan — cara verifikasinya didokumentasikan di bagian 4.

## 1. Masalah yang diselesaikan

Dua sekuens biologis (DNA atau protein) disejajarkan untuk mengukur
kemiripan mereka: residu yang sama/sejenis diberi skor positif,
substitusi yang jarang diberi skor negatif, dan sisipan/delesi (gap)
dikenai penalti. Alignment optimal adalah alignment dengan skor total
tertinggi. Ini adalah operasi dasar bioinformatika: dipakai untuk
mencari homologi, mengukur divergensi evolusioner, dan memvalidasi
rakitan sekuens.

## 2. Algoritma

### 2.1 Needleman–Wunsch — global alignment

Needleman & Wunsch (1970) memperkenalkan pemrograman dinamis untuk
alignment global: seluruh panjang kedua sekuens disejajarkan, dari
ujung ke ujung. Rekurensi dasar (tanpa gap affine):

```
F(i,j) = max( F(i-1,j-1) + s(a_i, b_j),
              F(i-1,j) - d,
              F(i,j-1) - d )
```

dengan `F(0,0) = 0` dan baris/kolom pertama diisi penalti gap kumulatif.
Kompleksitas waktu dan ruang O(mn). Implementasi di sini memakai gap
affine (lihat 2.4), sehingga rekurensinya memakai tiga matriks Gotoh,
bukan satu.

### 2.2 Smith–Waterman — local alignment

Smith & Waterman (1981) memodifikasi Needleman–Wunsch untuk alignment
lokal: alignment terbaik antara *substring* kedua sekuens. Perbedaannya
hanya dua: setiap sel dilantai di nol (`max(0, ...)`), dan traceback
dimulai dari sel bernilai maksimum di seluruh matriks, bukan dari sudut
kanan-bawah. Alignment lokal adalah yang dipakai BLAST dan pencari
motif: ia menemukan daerah mirip tanpa memaksa ujung-ujung yang tidak
terkait ikut disejajarkan.

### 2.3 Semi-global (overlap / fitting)

Varian yang gap di ujung alignment tidak dikenai penalti — berguna saat
satu sekuens (mis. read pendek) harus "pas" di dalam sekuens lain yang
lebih panjang. Implementasi di sini memakai varian *overlap*: gap di
awal dan akhir **kedua** sekuens gratis. Inisialisasi baris/kolom
pertama dengan 0 dan skor akhir diambil dari maksimum sel di baris
terakhir atau kolom terakhir (Durbin et al. 1998, Bab 2).

### 2.4 Gap affine (Gotoh 1982)

Gap dengan panjang k dikenai biaya `open + k·extend`, bukan `k·d`
linear. Ini mencerminkan biologi: membuka sebuah indel (satu peristiwa
mutasi) mahal, memperpanjangnya murah. Gotoh (1982) menunjukkan gap
affine dapat dihitung tetap dalam O(mn) dengan tiga matriks: M (berakhir
dengan pasangan residu), X (berakhir dengan gap di sekuens A), Y
(berakhir dengan gap di sekuens B). Konvensi biaya di sini —
`open + k·extend` — sama dengan konvensi NCBI BLAST `(a, b)`: gap
panjang 1 dengan biaya (10,1) bernilai −11, sesuai halaman matriks
substitusi NCBI ("a gap of k residues receives a total score of
-(a+bk)").

## 3. Matriks skor substitusi

### 3.1 BLOSUM — Henikoff & Henikoff (1992)

BLOSUM (BLOcks SUbstitution Matrix) dibangun dari frekuensi substitusi
yang *diamati langsung* di blok-blok ungapped dari database BLOCKS.
Angka di nama matriks adalah ambang identitas untuk clustering:
BLOSUM62 mengelompokkan sekuens yang ≥62% identik sehingga tiap
kelompok dihitung sekali. Entropi relatif BLOSUM62 adalah 0,6979 bit
(dari header file NCBI `blosum62.iij`); BLOSUM80 0,9868 bit
(`blosum80_3.iij`, satuan 1/3 bit). Aturan praktis (NCBI): BLOSUM45
untuk homologi jauh, BLOSUM62 default umum, BLOSUM80 untuk sekuens
sangat mirip. Biaya gap yang direkomendasikan NCBI untuk BLOSUM62
adalah (10,1) — dipakai sebagai default program ini.

### 3.2 PAM — Dayhoff, Schwartz & Orcutt (1978)

PAM (Point Accepted Mutation) dibangun dari model evolusi: matriks
probabilitas mutasi PAM-1 (1% residu berubah) dipangkatkan menjadi
PAM-250 (250 mutasi diterima per 100 residu), lalu diubah ke skor
log-odds. Berbeda dengan BLOSUM yang empiris, PAM adalah model
Markov ekstrapolasi — karena itu PAM250 cocok untuk hubungan jauh,
sementara PAM rendah (PAM30/70) untuk hubungan dekat. Skala file NCBI
yang dipakai: `ln(2)/3`, skor harapan −0,844, entropi 0,354 bit
(header file NCBI `PAM250`).

### 3.3 DNA

Untuk DNA program memakai skor match/mismatch sederhana (default
+1/−1, dapat diubah lewat flag `--match`/`--mismatch`), dengan basa
ambigu N bernilai 0 terhadap apa pun. Input RNA (U) dibaca sebagai T.
Default gap DNA (5,1).

## 4. Verifikasi nilai matriks (dilakukan 2026-11-10)

Matriks di `internal/matrix/data.go` diambil dari file NCBI yang
didistribusikan bersama Biopython
(`Bio/Align/substitution_matrices/data/`), lalu diperiksa silang
dengan tiga sumber independen. **Hasilnya: dua varian matriks yang
beredar di internet ternyata berbeda — proyek ini memakai varian
file NCBI dan mendokumentasikan perbedaannya, bukan menutupinya.**

| Matriks | Sumber pembanding | Hasil |
|---|---|---|
| BLOSUM45 | `scoring-matrices` (PyPI, data NCBI BLAST+) | 0 beda pada inti 20×20 |
| BLOSUM45 | CalState LA `blosum45.txt` (file lengkap 24×24) | **0/576 beda** |
| BLOSUM62 | `scoring-matrices` | 0 beda pada inti 20×20 |
| BLOSUM62 | Daftar diagonal DeepWiki/MELO + halaman CMU | diagonal cocok persis (A=4 … W=11 …); kolom T halaman CMU (T-T=4, P-T=1) **berbeda** dari file NCBI (T-T=5, P-T=−1) — halaman CMU dinilai salah salin |
| BLOSUM80 | `scoring-matrices` | 0 beda pada inti 20×20 |
| BLOSUM80 | Paket PyPI `blosum` | **berbeda 324 entri** — paket itu memakai varian lain (A-A=5, W-W=11); proyek memakai file NCBI (A-A=7, W-W=16, entropi 0,9868) |
| PAM250 | `scoring-matrices` | 0 beda pada inti 20×20 |
| PAM250 | Daftar GPM/Pearson 1990 (`thegpm.org`) | **berbeda 7 entri** (C-C=4, S-S=3, F-A=−4, F-N=−4, P-N=−1, P-G=−1, T-F=−2 vs NCBI C-C=12, S-S=2, …); varian GPM tidak dipakai |

Regenerasi matriks: `python3 scripts/generate_matrix_data.py`
(memerlukan Biopython; hasilnya harus identik dengan `data.go`).

Catatan Biostrings (R) yang relevan: definisi matriks BLOSUM/PAM
"tidak unik" — file NCBI FTP dan yang dibundel dalam BLAST
stand-alone bisa berbeda. Karena itu proyek ini mengunci satu sumber
(file NCBI) dan mencatat varian lainnya di atas.

## 5. Verifikasi algoritma (dilakukan 2026-11-10)

* **Hitung tangan** (lihat `internal/align/align_test.go`): skor
  kecil yang diverifikasi manual di komentar tes, mis. NW
  `AC` vs `A` = 1 − (5+1) = −5; NW `ACDE` vs `ACDE` BLOSUM62 =
  4+9+6+5 = 24; SW `HEAGAWGHEE` vs `PAWHEAE` = 17
  (HEA/HEA = 8+5+4, setara dengan AWGHE/AW-HE = 4+11+8+5−11).
* **Implementasi referensi independen** (`scripts/reference_align.py`,
  Python, ditulis terpisah): skor Go dan Python dibandingkan pada 15
  kasus acak (seed 20261110) untuk ketiga algoritma — semua cocok.
  Kasus klasik `GATTACA` vs `GCATGCT` (DNA +1/−1, gap 5/1):
  NW −1 (tanpa gap), SW 2 (`AT`/`AT`), semi-global 3.
* **Contoh nyata**: hemoglobin manusia α (P69905, 142 aa) vs β
  (P68871, 147 aa), NW BLOSUM62 gap (10,1): skor **286**, identitas
  **43,62% (65/149)**, similaritas 60,40% — sesuai literatur bahwa
  rantai α dan β globin memiliki ~44% identitas. α manusia vs α
  simpanse (P69907 residu 2–142): skor **717**, identitas **99,30%
  (141/142)** — hanya Met awal yang tidak ada pada catatan
  rekombinan simpanse.

## 6. Keputusan desain

* **Gap affine Gotoh, tiga matriks + traceback pointer.** Dipilih
  atas implementasi linear-gap yang lebih sederhana karena (a) sesuai
  praktik BLAST/NCBI, (b) tetap O(mn). Harga: memori 3×(m+1)(n+1)
  integer + 3 matriks byte — dapat diterima untuk CLI pairwise.
* **Tie-breaking deterministik dan terdokumentasi**: di matriks M,
  substitusi didahulukan atas gap; di matriks gap, memperpanjang gap
  didahulukan atas membuka gap baru (menggabungkan gap yang
  berdekatan). Alignment dengan skor sama tetap bisa berbeda bentuk;
  tes hanya mengunci skor (dan bentuk untuk kasus unik).
* **Matriks 24 huruf (termasuk B/Z/X/\*)** seperti file NCBI, bukan
  hanya 20: sekuens nyata mengandung X/B/Z dan harus bisa disejajarkan.
* **Paket terpisah** (`matrix`, `align`, `fasta`, `stats`): inti DP
  tidak bergantung pada CLI sehingga bisa dipakai sebagai library.
* **Default parameter mengikuti NCBI**: protein (10,1), DNA (5,1);
  matriks auto (DNA bila kedua sekuens hanya ACGTN, BLOSUM62 bila
  bukan). Pengguna tetap bisa mengubah semuanya lewat flag.

## 7. Keterbatasan

* Hanya pairwise alignment (bukan multiple alignment, bukan
  database search seperti BLAST — tidak ada heuristik seed-and-extend
  dan tidak ada statistik E-value/Karlin-Altschul).
* DP O(mn) penuh tanpa optimasi memori (tanpa Hirschberg/Myers-Miller):
  dua sekuens 10 kb butuh ~100 juta sel — lambat dan boros memori.
  Untuk sekuens sebesar itu, gunakan alat khusus (mis. minimap2).
* Semi-global di sini adalah varian *overlap* (gap ujung gratis di
  kedua sekuens); varian *fitting* satu arah tidak disediakan.
* Skor similaritas memakai ambang sederhana (skor substitusi > 0),
  bukan kelompok sifat kimia asam amino.
* Tidak ada penanganan urutan dengan karakter selain alfabet matriks
  (U/Selenosistein dan O/Pirolisin ditolak dengan pesan error).

## References

1. Needleman, S. B. & Wunsch, C. D. (1970). "A general method
   applicable to the search for similarities in the amino acid
   sequence of two proteins." *J. Mol. Biol.* 48(3):443–453.
   https://doi.org/10.1016/0022-2836(70)90057-4
2. Smith, T. F. & Waterman, M. S. (1981). "Identification of common
   molecular subsequences." *J. Mol. Biol.* 147(1):195–197.
   https://doi.org/10.1016/0022-2836(81)90087-5
3. Gotoh, O. (1982). "An improved algorithm for matching biological
   sequences." *J. Mol. Biol.* 162(3):705–708.
   https://doi.org/10.1016/0022-2836(82)90398-9
4. Henikoff, S. & Henikoff, J. G. (1992). "Amino acid substitution
   matrices from protein blocks." *Proc. Natl. Acad. Sci. USA*
   89(22):10915–10919. https://doi.org/10.1073/pnas.89.22.10915
5. Dayhoff, M. O., Schwartz, R. M. & Orcutt, B. C. (1978). "A model
   of evolutionary change in proteins." *Atlas of Protein Sequence
   and Structure*, Vol. 5, Suppl. 3, pp. 345–352. National Biomedical
   Research Foundation, Washington, DC.
6. Altschul, S. F. (1991). "Amino acid substitution matrices from an
   information theoretic perspective." *J. Mol. Biol.* 219:555–565.
   https://doi.org/10.1016/0022-2836(91)90193-A
7. NCBI BLAST substitution matrices (gap costs, matrix choice table).
   https://www.ncbi.nlm.nih.gov/blast/html/sub_matrix.html
8. NCBI BLAST matrix files (FTP).
   ftp://ftp.ncbi.nih.gov/blast/matrices/
9. Durbin, R., Eddy, S. R., Krogh, A. & Mitchison, G. (1998).
   *Biological Sequence Analysis: Probabilistic Models of Proteins
   and Nucleic Acids.* Cambridge University Press. (Bab 2: pairwise
   alignment; Bab 3: statistik alignment.)
10. Eddy, S. R. (2004). "Where did the BLOSUM62 alignment score matrix
    come from?" *Nat. Biotechnol.* 22:1035–1036.
    https://doi.org/10.1038/nbt0804-1035
11. UniProtKB P69905 (HBA_HUMAN) and P68871 (HBB_HUMAN), hemoglobin
    subunit alpha/beta, *Homo sapiens*.
    https://rest.uniprot.org/uniprotkb/P69905
    https://rest.uniprot.org/uniprotkb/P68871
12. UniProtKB P69907 (HBA_PANTR), hemoglobin subunit alpha,
    *Pan troglodytes*. https://rest.uniprot.org/uniprotkb/P69907
13. Biopython substitution matrices (data source for the tables used).
    https://biopython.org/docs/latest/api/Bio.Align.substitution_matrices.html
14. Bioconductor Biostrings, `substitution_matrices` documentation
    (on BLOSUM/PAM definitions not being unique across sources).
    https://bioconductor.org/packages/Biostrings/
