# Literatür — Modül 09: BLAST (viral contig -> virüs NT tanımlama)

> `✓` = doğrulandı

## Amaç
Viral (faj) contig'leri nükleotid düzeyinde NCBI virüs referans genomlarına hizalayıp
**en iyi hit ile taksonomik tanımlama** üretmek (geNomad/Pharokka'nın marker/mash tabanlı
sınıflamasını NT düzeyinde doğrulayan/tamamlayan adım).

## Araç & veritabanı
- **BLAST+ (blastn)** — Camacho, C., et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421. **✓** DOI: 10.1186/1471-2105-10-421
- **Veritabanı:** NCBI `ref_viruses_rep_genomes` (önceden formatlı viral referans genomlar = "virüs NT").
  Tüm `nt` (~250 GB) yerine viral-özgü set seçildi: disk/performans + viral tanımlamada daha temiz.
  - Güncelleme: `update_blastdb.pl ref_viruses_rep_genomes` (sürüm rapora yazılır).

## Akıllı filtre (orijinal prompt #10 ile uyumlu)
Pipeline yalnızca en yüksek bitscore'a bakmaz; **biyolojik anlamlılık** için ek işaret:
- En iyi hit = max bitscore (her contig için)
- `confident = yes` ancak **identity ≥ %70 VE query coverage ≥ %50** ise; aksi halde `low`
- Böylece kısa/parçalı bir contig'in büyük bir genomun küçük kısmına hizalanması "tam tanım" sayılmaz
  (prompt'taki 18 kb contig'in 500 kb kromozoma kısmi hizalanması örneğinin viral karşılığı).

> Eşikler (%70 identity / %50 qcov) muhafazakâr başlangıç; tür-içi/ötesi ayrım için
> ileride reciprocal coverage + ANI (skani) ile güçlendirilecek. ⏳

## Çıktı
- `blast.tsv` — tüm hitler (outfmt 6: qseqid sseqid pident length qcovs evalue bitscore staxids stitle)
- `best_hits.tsv` — her viral contig: en iyi hit başlığı, identity, qcov, e-value, bitscore, confident
- Bu tanımlamalar nihai rapora ("BLAST — Virüs NT En İyi Tanımlama") otomatik eklenir.

## Not
Bakteri/plazmid contig'leri için NT tanımlama ileride eklenebilir (uzak NCBI nt veya hedefli DB);
şimdilik kullanıcı talebi doğrultusunda **viral contig'lere** odaklanıldı.
