# DURUM — Ali WGS Pipeline (BacForge)

> "Nerede kaldık" anlık görüntüsü. Detay Claude belleğinde (`ali-wgs-pipeline` memory).

**Konum:** `/home/ali/ali-wgs-pipeline/` · git: `git@github.com:aliarslan47/BacForge.git`
**Son güncelleme:** 2026-08-11

## Şu an nerede kaldık
- **Milestone 1 (güvenilir uçtan-uca) TAMAM** — kısa-okuma yolu doğrulandı.
- Antigravity'nin ürettiği M00–M18 iskeleti **denetlendi**; ciddi fabrikasyon/stub bulundu ve **söküldü**.
- Pipeline artık **Acinetobacter baumannii** kısa-okuma örneğinde baştan sona koşuyor (exit 0), sonuçlar gerçek.

## Bu turda düzeltilenler (fabrikasyon → gerçek)
- **M04:** CheckM2 DB yolu düzeltildi (`*.dmnd`); sabit `99.5/0.2` söküldü → gerçek değer yoksa WARNING.
- **M02:** kraken2 raporu gerçekten parse ediliyor; tür `ncbi_species` olarak downstream'e taşınıyor.
- **M05:** "default strain" uydurması söküldü → gerçek FastANI (yerel referans); referans yoksa NA.
- **M18:** sabit-kodlu kurgu rapor (ST258/99.5/blaKPC/Klebsiella) tamamen söküldü → yalnız gerçek dashboard verisi.
- **M14/15/16:** sahte-PASS + "Bulunamadı;" tree söküldü → tek örnekte dürüst NOT_APPLICABLE.
- **M07/09/10/11/12/13:** "Bulunamadı satırı + otomatik PASS" söküldü → araç çalıştı/0-bulgu=PASS, hata=WARNING/NA; `.json` çıktıları eklendi.
- **M07 bug:** `res.returncode` (dict) → AttributeError çökmesi düzeltildi.
- **M17:** sabit "Klebsiella pneumoniae" varsayımı kaldırıldı → tür M02'den.

## Doğrulanmış demo koşusu (Acinetobacter, DRR035591 MiSeq ~125x)
- Tür: *A. baumannii* (kraken2) · Assembly 4.05 Mb, N50 **104 kb**, CheckM2 completeness **%100** (gerçek)
- ANI %97.67 (ATCC 17978) · MLST **ST571** · AMR: **blaOXA-66, armA, blaADC-82, sul1...** (15 gen)
- Durumlar dürüst: M02/M06 WARNING, M14/15/16 NA, kalan PASS · rapor fabrikasyonsuz

## Önemli ders
- Pipeline/SPAdes en baştan doğruydu. İlk seçilen E. coli mini-run'ları (DRR091879/DRR034572)
  **kötü kütüphaneydi** (duplication %60 / medyan 23x ama milyonluk spike). Kabul kapısı artık
  **test-assembly N50** (mapping breadth yeterli değil).

## Sonraki adım
- Acinetobacter için **long (ONT) + hybrid** örnekleri indir → aynı uçtan-uca doğrulama.
- **Milestone 2 (spec çok-araçlı zenginleştirme):** gerçek polishing (Medaka/Polypolish), M08 RGI/ResFinder/ABRicate,
  M10 PlasmidFinder, M11 IntegronFinder/MEFinder, M12 CRISPRCasFinder/DefenseFinder, M07 Kaptive/chewBBACA,
  M05 NCBI datasets ile otomatik closest-5, batch modda M14/15/16 (clinker/Panaroo/IQ-TREE).
- Kraken2 DB küçük → tür-düzeyi dağıtıyor (M02 WARNING); Bracken eklenebilir.

## Notlar
- Referans seti: `databases/references/acinetobacter/GCF_000015425_ATCC17978.fna` (M05 buradan FastANI).
- Kullanılmayan faj DB'leri (checkv/pharokka/blast_viral/viral_proteins, ~8.9G) duruyor — silinsin mi kullanıcı kararı.
- Eski çakışan MVP modülleri (m01_input, m03_qc, m0X_* faj seti) hâlâ `modules/` içinde kayıtsız/ölü — temizlenebilir.
