# Literatür — Modül 07: Contig Filtreleme (uzunluk + coverage)

> `✓` = doğrulandı · `⏳` = doğrulanacak

## Amaç
De novo assembly çıktısındaki düşük-güvenilirlikli contig'leri (kısa parçalar,
düşük-coverage artefaktlar/kontaminasyon) downstream analizden önce elemek.

## Uygulanan eşikler (config.contig_filter)
| Kriter | Varsayılan | Gerekçe |
|--------|-----------|---------|
| min_length | 1000 bp | <1 kb contig'ler tipik olarak parçalı/tekrar kaynaklı, düşük güven. QUAST raporlaması ≥500 bp eşiğini standart alır; downstream analiz pipeline'larında ≥1 kb yaygın. |
| min_coverage_abs | 3x | Çok düşük mutlak derinlik = olası hata/kontaminasyon. |
| min_coverage_frac | 0.10 × medyan | Genomun modal/medyan coverage'ının çok altındaki contig'ler farklı kaynaklı (kontaminasyon) olabilir. |

> Coverage filtresi yalnızca coverage bilgisi mevcutsa (Flye `assembly_info.txt`) uygulanır;
> yoksa sadece uzunluk filtresi (ör. FASTA girdi, kısa-okuma).

## Literatür dayanağı
- **QUAST** raporlama standardı contig'leri ≥500 bp eşiğiyle değerlendirir:
  Gurevich, A., et al. (2013). QUAST. *Bioinformatics*, 29(8), 1072–1075. ✓ DOI: 10.1093/bioinformatics/btt086
- Düşük-coverage contig elemesi, kontaminasyon/yanlış-birleştirme azaltımı için yaygın pratik
  (ör. metagenom/izolat assembly QC akışları). Eşiğin medyan-oransal + mutlak taban kombinasyonu,
  tek sabit eşiğe göre tür/coverage bağımsız daha sağlamdır. ⏳ (spesifik benchmark atfı eklenecek)
- Genom bütünlüğü/kontaminasyon nicel değerlendirmesi ayrıca Modül 17'de CheckM2/CheckV ile yapılır
  (Chklovski 2023 ✓; Nayfach 2021 ✓), yani 07 bir ÖN filtredir, nihai kalite kararı 17'dedir.

## Çıktı
- `contigs.filtered.fasta` — geçen contig'ler
- `filter_report.tsv` — her contig: uzunluk, coverage, tutuldu mu, sebep (şeffaflık/izlenebilirlik)
- `filter_params.txt` — kullanılan eşikler + kaç/kaç contig tutuldu (rapora yazılır)
