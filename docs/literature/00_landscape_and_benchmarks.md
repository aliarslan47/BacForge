# Literatür — Bölüm 0: Genel Manzara & Benchmark Çalışmaları

> Atıf doğrulama: `✓` = DOI/PMID web aramasıyla doğrulandı · `⏳` = doğrulanacak

## 0.1 Neden yeni bir karar çerçevesine ihtiyaç var?

Bakteriyel WGS alanında 2022–2025 arasında **iki paradigma kayması** yaşandı ve bu, eski pipeline'ların (çoğu Illumina + SPAdes varsayımı üzerine kurulu) önerilerini geçersizleştirdi:

1. **ONT R10.4.1 + güncel basecalling**, kısa-okuma desteği olmadan **near-perfect** (≈Q50, 1 hata/100 kbp) bakteriyel genom üretebiliyor. Bu, "Illumina olmadan tamamlanmış genom olmaz" dogmasını yıktı.
2. **Konsensüs/uzlaşı assembly** (Trycycler → Autocycler) ve ML-tabanlı kalite değerlendirme (CheckM2), tek-araç assembly + CheckM1 yaklaşımını doğrulukta geride bıraktı.

Bu yüzden **BacForge'ın araç seçimi statik değil, platform + okuma profiline göre dallanan bir karar ağacı** olarak tasarlanacaktır.

## 0.2 Temel referans çalışmalar (doğrulanmış)

**[A] Wick, R. R., Judd, L. M., & Holt, K. E. (2023). Assembling the perfect bacterial genome using Oxford Nanopore and Illumina sequencing. *PLOS Computational Biology*, 19(3), e1010905.** ✓
DOI: 10.1371/journal.pcbi.1010905 · PMC9980784
→ Bu pipeline'ın **assembly + polishing omurgasının** ana referansı. Önerilen akış: Trycycler long-read assembly → Medaka long-read polishing → Polypolish short-read polishing → (gerekirse) ek kısa-okuma polishing → manuel kürasyon. ALE ile assembly versiyon seçimi.

**[B] Sereika, M., Kirkegaard, R. H., Karst, S. M., et al. (2022). Oxford Nanopore R10.4 long-read sequencing enables the generation of near-finished bacterial genomes from pure cultures and metagenomes without short-read or reference polishing. *Nature Methods*, 19, 823–826.** ✓
DOI: 10.1038/s41592-022-01539-7
→ ONT-only finishing'in mümkün olduğunun kanıtı. Hybrid'i her durumda zorunlu kılmama gerekçemiz.

**[C] (R10.4.1 doğruluk değerlendirmesi) Microbial Genomics (2024), 10, 001246.** ✓
DOI: 10.1099/mgen.0.001246 · PMC11170131
→ R10.4.1 long-read-only assembly'lerin iyileştiğini ama hâlâ Illumina-ONT hybrid'e göre daha fazla hata içerebildiğini gösterir → **hybrid veri varsa polishing'e dahil et** kuralının dayanağı. (Yazar listesi & cilt detayı ⏳ doğrulanacak.)

**[D] Wick, R. R., Judd, L. M., Cerdeira, L. T., et al. (2021). Trycycler: consensus long-read assemblies for bacterial genomes. *Genome Biology*, 22, 266.** ✓
DOI: 10.1186/s13059-021-02483-z · PMID: 34521459 · PMC8442456
→ Konsensüs assembly'nin tek-araç assembly'den daha az hata içerdiğinin kanıtı.

**[E] Wick, R. R. (2025). Autocycler: long-read consensus assembly for bacterial genomes. *Bioinformatics*, 41(9), btaf474.** ✓
DOI: 10.1093/bioinformatics/btaf474
→ Trycycler'ın **otomatikleştirilmiş** halefi. Trycycler'ın manuel adımları otomasyona uygun değildi; Autocycler tam-otomatik pipeline hedefimizle birebir uyumlu. **Aday birincil konsensüs-assembly motoru.**

> **Önemli tasarım sonucu:** Pipeline tam otomatik olacağı için, manuel müdahale gerektiren Trycycler yerine **Autocycler** birincil aday; Trycycler referans/doğrulama amaçlı korunur. (Nihai karar FAZ 2.)

## 0.3 Benchmark çalışmaları (assembler karşılaştırmaları)

- **ONT bakteriyel patojen assembler benchmark** — Raven ve Flye'ın dengeli/sağlam performansı; Raven hızlı, Flye yüksek tamamlanma + düşük misassembly. (PMC7730629 ✓; tam atıf ⏳)
- **Metagenomik long-read benchmark** — Flye dengeli (runtime ~0.3 h, yüksek completeness). (PMC9861289 ✓; tam atıf ⏳)
- **2025 long-read assembler benchmark (E. coli DH5α)** — NextDenovo/NECAT en bütün/sürekli; Flye doğruluk-hız-bütünlük dengesi. (ScienceDirect S2215017X2500058X ⏳)

**Ön çıkarım (FAZ 2'de kesinleşecek):**
- ONT-only birincil assembler adayı: **Flye** (denge) + konsensüs için **Autocycler**.
- Illumina-only: **SPAdes/SKESA** (aşağıda Modül 5).
- PacBio HiFi: **hifiasm** (aşağıda Modül 5).
- Hybrid: **Unicycler** veya long-read-first + short-read polishing (Polypolish).

## 0.4 Kaynaklar (web doğrulaması)
- https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010905
- https://www.nature.com/articles/s41592-022-01539-7
- https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001246
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11170131/
- https://link.springer.com/article/10.1186/s13059-021-02483-z (PMID 34521459)
- https://academic.oup.com/bioinformatics/article/41/9/btaf474/8242761
- https://pmc.ncbi.nlm.nih.gov/articles/PMC7730629/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC9861289/
