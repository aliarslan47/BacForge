# Literatür — Faj (Bakteriyofaj) Araç Zinciri

> Pipeline içerik-farkında: geNomad bir contig'i virüs/faj olarak sınıflarsa,
> annotation/completeness/taxonomy OTOMATİK olarak faj araçlarına yönlenir.
> `✓` = DOI/PMID web doğrulandı.

## Neden ayrı zincir?
Faj viral bir varlıktır; bakteriyel araçlar (Bakta gen modeli, CheckM2 marker setleri, MLST)
faja uygulanınca anlamsız/yanlış sonuç verir. Faj için özel, literatürle doğrulanmış araçlar:

## R1. İçerik sınıflandırıcı (router) — geNomad
- **Amaç:** herhangi bir dizide virüs/plazmid/kromozom ayrımı (otomatik yönlendirici).
- **İlk makale:** Camargo, A. P., Roux, S., Schulz, F., et al. (2023). Identification of mobile genetic elements with geNomad. *Nature Biotechnology*. **✓** DOI: 10.1038/s41587-023-01953-y
- **GitHub:** https://github.com/apcamargo/genomad
- **Benchmark:** virüs sınıflamasında MCC %95.3, plazmidde %77.8 — diğer araçları belirgin geçiyor.
- **Tercih gerekçesi:** "ne koyarsan tanı" hedefimizin merkezi; tek araçla virüs+plazmid+kromozom.

## R2. Faj annotation — Pharokka
- **Amaç:** faj genomu hızlı/standart annotation (CDS, tRNA, tmRNA, CRISPR; PHROGs ile fonksiyon).
- **İlk makale:** Bouras, G., Nepal, R., Houtak, G., et al. (2023). Pharokka: a fast scalable bacteriophage annotation tool. *Bioinformatics*, 39(1), btac776. **✓** DOI: 10.1093/bioinformatics/btac776 · PMC9805569
- **GitHub:** https://github.com/gbouras13/pharokka
- **Alternatif/üst:** phold (yapısal protein ile annotation, aynı geliştirici) — ileride değerlendirilecek.
- **Tercih gerekçesi:** Bakta'nın faj karşılığı; PHROGs faj veritabanı + faj-özgü özellikler.

## R3. Faj completeness/kalite — CheckV
- **Amaç:** viral genom tamlığı/kontaminasyon, kapalı genom tespiti, host bölgesi temizleme.
- **İlk makale:** Nayfach, S., Camargo, A. P., Schulz, F., et al. (2021). CheckV assesses the quality and completeness of metagenome-assembled viral genomes. *Nature Biotechnology*, 39, 578–585. **✓** DOI: 10.1038/s41587-020-00774-7 · PMID: 33349699
- **GitHub:** https://bitbucket.org/berkeleylab/checkv
- **Tercih gerekçesi:** CheckM2'nin faj karşılığı; faj genom tamlığının standardı.

## R4. Faj yaşam tarzı — BACPHLIP
- **Amaç:** litik (virülan) vs lizojenik (temperate) tahmini (korunmuş protein domainleri + Random Forest).
- **İlk makale:** Hockenberry, A. J., & Wilke, C. O. (2021). BACPHLIP: predicting bacteriophage lifestyle from conserved protein domains. *PeerJ*, 9, e11396. **✓** DOI: 10.7717/peerj.11396 · PMID: 33996289
- **GitHub:** https://github.com/adamhockenberry/bacphlip
- **Benchmark:** bağımsız 423 fajda %98 doğruluk (önceki araçlar %79).
- **Tercih gerekçesi:** Lizojenik fajların ARG/virülans taşıma riski için yaşam tarzı kritik.

## Faj zinciri akışı (otomatik)
```
contig --geNomad--> "virus" ise:
   Pharokka (annotation) -> CheckV (tamlık) -> BACPHLIP (yaşam tarzı)
   + AMRFinderPlus (lizojenik ARG taraması, ANLAMLI) + (ileride) PhageTerm (genom uçları), taxmyPHAGE (taksonomi)
contig --geNomad--> "chromosome/plasmid" ise:
   bakteri zinciri (Bakta -> CheckM2 -> MLST -> ...)
```

## Faj-özgü ön işleme notu (KRİTİK)
Faj okumaları çok yüksek coverage verir (ör. 58 Mbp / ~50 kb genom ≈ 1000x).
Aşırı coverage Flye assembly'sini bozar → **hedef ~100x'e subsample** (rasusa/Filtlong target_bases).
Bu adım `config.tools.subsample` ile otomatik. (Gerekçe: long-read assembler'lar aşırı derinlikte
yanlış-birleştirme/parçalanma eğilimi gösterir; ⏳ benchmark atfı FAZ 2'de eklenecek.)

## Kaynaklar
- geNomad: https://www.nature.com/articles/s41587-023-01953-y
- Pharokka: https://academic.oup.com/bioinformatics/article/39/1/btac776/6858464
- CheckV: https://pubmed.ncbi.nlm.nih.gov/33349699/
- BACPHLIP: https://peerj.com/articles/11396/
