# Literatür — Modül 05: Assembly & Polishing (pipeline'ın en kritik kararı)

> `✓` = DOI web doğrulandı · `⏳` = doğrulanacak

## 5.0 Karar çerçevesi (platforma göre dallanma)

`01_Input` çıktısına göre assembler **otomatik** seçilir. Hiçbir senaryoda tek bir "evrensel assembler" yoktur — bu, literatürün net bulgusudur (Wick et al. 2023; assembler benchmark'ları).

| Senaryo | Birincil strateji | Polishing | Gerekçe (ref) |
|---------|-------------------|-----------|----------------|
| **ONT-only** (R10.4.1) | Flye → **Autocycler** konsensüs | Medaka | Wick 2023 [A], Sereika 2022 [B], Autocycler [E] |
| **Illumina-only** | SPAdes (`--isolate`/careful) veya SKESA | — | Souvorov 2018, Prjibelski 2020 |
| **PacBio HiFi** | hifiasm | (gerekmez, Q30+) | Cheng 2021 |
| **Hybrid (ONT+Illumina)** | long-read-first (Flye/Autocycler) + kısa-okuma polishing | Medaka → **Polypolish** (+Polypolish-careful/POLCA) | Wick 2023 [A], Polypolish |
| **Hybrid (alternatif)** | Unicycler (köprüleme) | dahili | Wick 2017 |

> `[A][B][E]` referansları `00_landscape_and_benchmarks.md` dosyasında.

---

## 5.1 Flye — ONT/PacBio long-read birincil assembler

- **Amaç:** repeat-graph tabanlı tek-molekül (ONT/PacBio) de novo assembly.
- **Geliştiriciler:** Mikhail Kolmogorov, Pavel Pevzner ve ark. (UCSD).
- **İlk makale:** Kolmogorov, M., Yuan, J., Lin, Y., & Pevzner, P. A. (2019). Assembly of long, error-prone reads using repeat graphs. *Nature Biotechnology*, 37, 540–546. **✓** DOI: 10.1038/s41587-019-0072-8
- **İlgili:** metaFlye — Kolmogorov et al. (2020). *Nature Methods*, 17, 1103–1110. **✓** DOI: 10.1038/s41592-020-00971-x
- **GitHub:** https://github.com/mikolmogorov/Flye · **Güncel sürüm:** 2.9.x (aktif)
- **Avantajlar:** denge (hız + tamamlanma), düşük misassembly, circular contig işaretleme (`circular=Y`), plasmid/metagenom desteği.
- **Dezavantajlar:** ONT homopolimer hatalarına karşı polishing gerektirir; çok düşük coverage'da zayıflar.
- **Benchmark:** patojen ONT benchmark'larında yüksek completeness + düşük misassembly (PMC7730629 ✓); metagenom benchmark'ında dengeli (PMC9861289 ✓).
- **Alternatifler:** Raven (hızlı, ARG için doğru), Canu (yavaş ama hassas), NextDenovo/NECAT (en sürekli — 2025 benchmark ⏳), Shasta.
- **Tercih gerekçesi:** otomatik pipeline için en güvenilir denge; Autocycler içine girdi assembler olarak da idealdir.

## 5.2 Autocycler — otomatik konsensüs assembly (birincil aday)

- **Amaç:** birden çok long-read assembly'den otomatik konsensüs üretmek (Trycycler'ın otomasyon halefi).
- **Geliştirici:** Ryan R. Wick.
- **İlk makale:** Wick, R. R. (2025). Autocycler: long-read consensus assembly for bacterial genomes. *Bioinformatics*, 41(9), btaf474. **✓** DOI: 10.1093/bioinformatics/btaf474
- **GitHub:** https://github.com/rrwick/Autocycler · **Aktif**
- **Avantaj:** Trycycler'ın doğruluğunu **manuel adım olmadan** sağlar → tam otomatik pipeline ile uyumlu.
- **Dezavantaj:** birden çok assembler çalıştırma maliyeti (hesaplama ağır); yeni araç (uzun vadeli benchmark birikimi sınırlı).
- **Alternatif:** Trycycler (Wick 2021 ✓, manuel — referans/doğrulama amaçlı).
- **Tercih gerekçesi:** "tam otomatik + en az hata" hedefimizle birebir.

## 5.3 hifiasm — PacBio HiFi birincil assembler

- **Amaç:** HiFi okumalar için faz-çözümlü/doğru de novo assembly.
- **Geliştiriciler:** Haoyu Cheng, Heng Li ve ark.
- **İlk makale:** Cheng, H., Concepcion, G. T., Feng, X., Zhang, H., & Li, H. (2021). Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. *Nature Methods*, 18, 170–175. **✓** DOI: 10.1038/s41592-020-01056-5
- **GitHub:** https://github.com/chhylp123/hifiasm · **Aktif**
- **Avantaj:** HiFi'de en yüksek doğruluk/süreklilik; bakteride genelde polishing gereksiz (Q30+).
- **Dezavantaj:** ONT/Illumina için uygun değil (HiFi'ye özel).
- **Tercih gerekçesi:** HiFi standardı.

## 5.4 SPAdes / SKESA — Illumina short-read assembler

- **SPAdes — İlk makale:** Bankevich, A. et al. (2012). SPAdes... *Journal of Computational Biology*, 19(5), 455–477. ⏳ DOI: 10.1089/cmb.2012.0021 · **Güncel:** Prjibelski et al. (2020), *Current Protocols in Bioinformatics* ⏳. `--isolate`/`--careful` modları bakteri izolatı için.
- **SKESA — İlk makale:** Souvorov, A., Agarwala, R., & Lipman, D. J. (2018). SKESA: strategic k-mer extension for scrupulous assemblies. *Genome Biology*, 19, 153. **✓** DOI: 10.1186/s13059-018-1540-z
- **Avantaj (SKESA):** hızlı, konservatif (düşük yanlış-birleştirme), NCBI patojen surveillance standardı; **Dezavantaj:** SPAdes'ten daha parçalı olabilir.
- **Tercih gerekçesi:** Illumina-only'de SKESA (güvenli, hızlı) varsayılan; karmaşık genomlarda SPAdes alternatif. FAZ 2'de benchmark ile sabitlenecek.

## 5.5 Unicycler — hybrid köprüleme (alternatif)

- **İlk makale:** Wick, R. R., Judd, L. M., Gorrie, C. L., & Holt, K. E. (2017). Unicycler: Resolving bacterial genome assemblies from short and long sequencing reads. *PLOS Computational Biology*, 13(6), e1005595. **✓** DOI: 10.1371/journal.pcbi.1005595 · PMID: 28594827
- **Not:** Modern hybrid önerisi (Wick 2023 [A]) artık **long-read-first + Polypolish**'i tercih ediyor; Unicycler düşük long-read derinliğinde alternatif olarak korunur.

## 5.6 Polishing araçları

- **Medaka** (ONT) — ONT resmi neural polisher. GitHub: nanoporetech/medaka. Long-read-only assembly'de birincil. ⏳ (atıf: ONT, makale yok — dokümante edilecek)
- **Polypolish — İlk makale:** Wick, R. R., & Holt, K. E. (2022). Polypolish: Short-read polishing of long-read bacterial genome assemblies. *PLOS Computational Biology*, 18(1), e1009802. **✓** DOI: 10.1371/journal.pcbi.1009802 · GitHub: rrwick/Polypolish. Tekrar bölgelerindeki hataları **all-per-read alignment** ile düzeltir, neredeyse hiç hata sokmaz.
- **Polypolish-careful / POLCA / Pilon** — ek kısa-okuma cila adımları (Wick 2023 [A] sıralamasına göre). ⏳

## 5.7 Kaynaklar
- Flye: https://www.nature.com/articles/s41587-019-0072-8 · metaFlye: https://experiments.springernature.com/articles/10.1038/s41592-020-00971-x
- Autocycler: https://academic.oup.com/bioinformatics/article/41/9/btaf474/8242761
- hifiasm: https://www.nature.com/articles/s41592-020-01056-5
- SKESA: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1540-z
- Unicycler: https://pubmed.ncbi.nlm.nih.gov/28594827/
- Polypolish: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009802
