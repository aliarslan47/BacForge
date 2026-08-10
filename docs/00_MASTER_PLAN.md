# Ali WGS Pipeline — Master Plan & Metodoloji

> **Durum:** FAZ 1 — Literatür Araştırması & Mimari Tasarım (kod yazımı YOK)
> **Son güncelleme:** 2026-06-26
> **Kapsam:** Bakteriyel Whole Genome Sequencing (WGS) — ONT, Illumina, PacBio HiFi, Hybrid

---

## 0. Çalışma İlkeleri (değişmez kurallar)

1. **Önce literatür, sonra kod.** Hiçbir araç/parametre, en az bir hakemli yayınla gerekçelendirilmeden seçilmez.
2. **Popülerlik ≠ seçim kriteri.** Seçim; benchmark performansı, doğruluk, biyolojik anlamlılık ve aktif geliştirilme durumuna dayanır.
3. **Öncelik 2022–2026** yayınları (Nature, Nat Methods, Nat Biotechnol, Genome Biology, Bioinformatics, NAR, Microbial Genomics, Briefings in Bioinformatics, hakemli benchmark çalışmaları).
4. **Atıf doğruluğu kutsaldır.** Uydurma DOI/PMID yasak. Her atıf web'den doğrulanır:
   - `✓` = DOI/PMID web aramasıyla doğrulandı
   - `⏳` = henüz doğrulanmadı (bir sonraki geçişte doğrulanacak)
5. **Reprodüksibilite:** sürüm sabitleme (pin), konteynerleştirme (Docker/Singularity), veritabanı sürüm kaydı zorunlu.
6. **Hiçbir dosya silinmez / üzerine yazılmaz.** Pipeline herhangi bir adımdan resume edebilmeli.

---

## 1. Faz Planı

| Faz | İçerik | Çıktı | Durum |
|-----|--------|-------|-------|
| **FAZ 1** | Literatür araştırması (19 modül + benchmark landscape) | `docs/literature/*.md` | 🔄 devam ediyor |
| **FAZ 2** | Araç seçim kararları + gerekçe matrisi | `docs/02_TOOL_DECISIONS.md` | ⏳ |
| **FAZ 3** | Pipeline mimarisi + DAG + karar ağaçları (otomatik platform tespiti) | `docs/03_ARCHITECTURE.md` | ⏳ |
| **FAZ 4** | Konteyner/ortam + veritabanı sürüm matrisi | `docs/04_ENV_DB.md` | ⏳ |
| **FAZ 5** | İmplementasyon (Nextflow/Snakemake kararı dahil) | kod | ⏳ (onay sonrası) |

> **Kod yazımı yalnızca FAZ 1–4 tamamlanıp kullanıcı onayı alındıktan sonra başlar.**

---

## 2. Modül Haritası (klasör yapısı = analiz akışı)

```
00_Project         Proje metadata, config, sample sheet
01_Input           Otomatik platform/okuma-tipi/pair tespiti
02_Preprocessing   Basecalling kontrol, format dönüşümü, demux
03_QC              Platforma özgü kalite kontrol
04_Filtering       Kalite/uzunluk/adaptör/kontaminasyon filtresi
05_Assembly        Platforma özgü assembler seçimi
06_Assembly_QC     QUAST, BUSCO/CheckM2, coverage, hata oranı
07_Contig_Filtering Uzunluk/coverage/güvenilirlik filtresi
08_Taxonomy        Kraken2/GTDB-Tk/Mash/sourmash/ANI
09_BLAST           Akıllı (reciprocal coverage) BLAST
10_Annotation      Bakta/Prokka/PGAP/DFAST
11_AMR             AMRFinderPlus/CARD-RGI/ResFinder
12_VFDB            Virülans (VFDB + abricate/AMRFinderPlus --plus)
13_Plasmid         MOB-suite/PlasmidFinder/platon/geNomad
14_Phage           geNomad/VIBRANT/PhiSpy/PHASTEST + CheckV
15_MLST            mlst/PubMLST + chewBBACA (cgMLST)
16_ANI             FastANI/skani + GTDB-Tk
17_Completeness    CheckM2 / BUSCO
18_Final_Report    HTML + PDF (References APA otomatik)
19_Logs            Her adım: log, params, tool, sürüm, süre
```

**Garantiler:** her klasör yalnızca kendi çıktısını üretir · ara dosyalar korunur · her adım `19_Logs`'a (tool adı, sürüm, parametre, komut satırı, çalışma süresi, exit code, veritabanı sürümü) yazar.

---

## 3. Otomatik Platform Tespiti — karar mantığı (taslak)

Pipeline kullanıcıdan seçim istemez. `01_Input` aşağıdaki sinyallerden platformu çıkarır:

| Sinyal | ONT | Illumina | PacBio HiFi |
|--------|-----|----------|-------------|
| Dosya uzantısı | `.fast5/.pod5/.fastq` | `.fastq` (R1/R2) | `.bam` (ccs) / `.fastq` |
| Okuma uzunluğu dağılımı | geniş, uzun (N50 ≫ 1 kb) | sabit, kısa (≤300 bp) | uzun + dar (10–25 kb) |
| Kalite profili | Phred orta, homopolimer hata | Phred yüksek, 3' düşüş | Phred çok yüksek (Q30+) |
| Header desenleri | `runid=`, `ch=`, `start_time=` | Illumina instrument ID | `ccs`, `zmw` etiketleri |
| Pair tespiti | tek dosya | R1+R2 eşleşmesi, `/1 /2` | tek dosya |

> Tespit mantığının her eşik değeri FAZ 3'te literatürle sabitlenecek (bkz. `literature/01_input_and_qc.md`).

---

## 4. Atıf Formatı Standardı

Her araç için literatür dosyasında zorunlu alanlar:

```
Araç: <ad> v<sürüm>
Amaç / Kullanım alanı:
Geliştiriciler:
İlk makale (APA + DOI + PMID + doğrulama durumu):
Güncel sürüm makalesi (varsa):
GitHub:
Resmî dokümantasyon:
Veritabanı sürümü (varsa):
Avantajlar / Dezavantajlar:
Alternatifler:
Benchmark performansı (hangi çalışma, hangi metrik):
Bu pipeline'da tercih gerekçesi:
```

Nihai raporun `References` bölümü bu alanlardan **otomatik** üretilecektir.
