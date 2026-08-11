# DURUM — BacForge

> "Nerede kaldık" anlık görüntüsü. Detay Claude belleğinde (`bacforge` memory).

**Konum:** `/home/ali/BacForge/` · git: `git@github.com:aliarslan47/BacForge.git`
**Son güncelleme:** 2026-08-11

## 2026-08-11 — M05/M16 AKRABALIK DOĞRULUĞU + per-contig FASTA + BLAST tanısı
Kullanıcı endişesi: "farklı contig'ler aynı türün farklı suşuna hit ediyor → akrabalık haritası yanlış çıkar mı?"
Cevap+düzeltme: kimlik/akrabalık per-contig BLAST'la DEĞİL, genom-geneli ANI/mash ile → per-contig gürültüsünden bağımsız.
GERÇEK KUSUR bulundu ve düzeltildi:
- **closest_5 KİRLİYDİ:** bakta `_bakta/ref.fna` kopyaları sayılıyordu → aynı genom 2× (Rank1=2, 3=4), accession "ref".
  M05: `_is_annotation_copy` + `_accession(GCF_/GCA_)` + `_parse_rank_fastani` (accession'a göre DEDUP) → **tekil** closest_5 + yeni **closest_10** (akrabalık haritası girdisi).
- **M16 ağacı** aynı kirliliği alıyordu (havuz `rglob("*.fna")` bakta kopyalarını + "ref" yapraklarını ekliyordu) → dedup + anotasyon-kopya filtresi. Ağaç artık **13 TEKİL yaprak** (QUERY + 5 ANI-etiketli + 7 çeşitlilik), "ref"/tekrar yok.
- **Per-contig FASTA:** her contig ayrı dosya, ANA dosya (genome.fasta) yanında `M04.../contigs/` (manuel web BLAST yükleme için) + `contig_lengths.tsv`; toplu genome.fasta'ya dokunulmadı. M05 remote-BLAST için 1Mb-guard'lı contig seçimi.
- **BLAST tanısı:** 328kb/16S `blastn -remote` bu ortamda sürekli timeout/boş (NCBI throttle, kod değil) → kimlik zaten kraken2(%94.83)+FastANI(ANI%99.5)'e dayanıyor, bloklamıyor.
- Koşu 20260811_165209_short: M05 closest temizlendi, M16→M17→M18 yeniden koşturuldu, rapor tekil accession'larla güncel (report.html 2.0MB). commit'ler: b14ed4f + bu tur.

## 2026-08-11 — PROJE ADI: ali-wgs-pipeline → BacForge (her yerde)
Dizin `/home/ali/ali-wgs-pipeline` → `/home/ali/BacForge`; Python paketi `ali_wgs` → `bacforge`
(42 dosya, import'lar + `python3 -m bacforge.cli`); env değişkenleri `ALI_WGS_*` → `BACFORGE_*`;
pyproject name/script `ali-wgs` → `bacforge`; tüm dokümanlarda "Ali WGS Pipeline" → "BacForge";
Claude bellek dosyası `ali-wgs-pipeline.md` → `bacforge.md` + `settings.local.json` path'leri güncellendi.
Doğrulama: `import bacforge` OK, `cli info` path'leri `/home/ali/BacForge`'a çözüyor.

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


## 2026-08-11 (gece) — RAPOR & 7 GERI BILDIRIM DUZELTMESI (commit a047bc2)
Kullanici geri bildirimleriyle: (1) genom haritasi referans-sirali tek psodo-kromozom -> dairesel harita
0.0 Mbp fix (minimap2 + 2. Bakta plot). (2) MLST Lokus|Alel tablosu (ST571). (3) M13 snippy GenBank ref
-> 8322 kodlayan varyant gen/EFFECT dolu (dnaA synonymous vb.). (4) M14 clinker QUERY daima ustte + en
yakin 5 tur. (5) M16 cesitlilik havuzu -> 13 yaprakli anlamli agac. (6) M05 16S fallback KALDIRILDI, daima
BLAST+bekle (3600s); M02 Bracken -> A.baumannii %94.83 PASS (WARNING gitti). (7) M18 klasoru=SADECE rapor
(zip/refs run kokune); Figure 1 BacForge'a OZGU 'dallan->birles->yelpaze' semasi (nf-core tarzi DEGIL, kural).
Rapor: 14 numarali Table + 5 Figure, pipeline calisma sirasina gore. Kosu: runs/20260811_165209_short.

## 2026-08-11 (aksam) — TAM UCTAN-UCA DOGRULANDI (Acinetobacter, katı sözleşme)
Final re-run: 19 modul, hepsi gercek ciktiyla. M06 Bakta --skip-sorf segfault fix pipeline'da tuttu (gercek GenBank 9.6MB).
M05 species-agnostic closest-5 gercek (BLAST kimlik + datasets + FastANI; %99.5 ANI). M14 HEDEFLI clinker (AMR+virulans,
blaOXA-66/ompA...) 5 gercek sinteni HTML (15sn). M16 mash-NJ gercek akrabalik agaci. Statüler dürüst: M02 WARNING (kucuk
kraken DB), M15 NOT_APPLICABLE (pangenom >=2 ornek). Katı sözleşme: tool exit 0 + gercek cikti yoksa PASS yok; mock fallback'ler silindi.
Koşu: runs/20260811_165209_short. commit d7014fd.
BILINEN: M05 uzun-contig remote BLAST timeout'ta conda-run torunu (blastn) orphan kaliyor -> surec-grubu ile oldurulmeli (temizlik).

## Notlar
- Referans seti: `databases/references/acinetobacter/GCF_000015425_ATCC17978.fna` (M05 buradan FastANI).
- **FAJ YOLU PLANLI (kullanıcı 2026-08-11): faj analizini DE yapacağız → SİLME.** Faj DB'leri
  (checkv/pharokka/blast_viral/viral_proteins, ~8.9G) KALACAK. Eski faj modülleri (m09_blast=viral BLAST,
  m10_annotation=Pharokka, m17_completeness=CheckV) ölü kod DEĞİL, faj branch'inin temeli — silinmeyecek,
  BacForge'a **içerik-farkında yönlendirme** ile geri entegre edilecek: geNomad (M12) her contig'i
  chromosome/plasmid/virus sınıflar → bakteri kolu (Bakta+CheckM2+MLST) | faj kolu (Pharokka+CheckV+viral BLAST).
  Bu, orijinal spec'teki routing fikriyle + eski 33-faj projesinin makinesiyle örtüşür.
