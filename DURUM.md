# DURUM — BacForge

> "Nerede kaldık" anlık görüntüsü. Detay Claude belleğinde (`bacforge` memory).

**Konum:** `/home/ali/BacForge/` · git: `git@github.com:aliarslan47/BacForge.git`
**Son güncelleme:** 2026-08-14

## 2026-08-14 (2) — MILESTONE 2 ÇOK-ARAÇLI ZENGİNLEŞTİRME (8 araç) + gerçek short koşu
Her modüle mevcut aracın yanına 2./3. araç eklendi (kaynak-etiketli, dürüst durum, rapora yansıyor):
- **M08 AMR:** + RGI/CARD (`rgi main`) + ResFinder (`abricate --db resfinder`) → amr.json by_source+source_counts,
  M18 AMR tablosuna "Kaynak" sütunu. **M10:** + PlasmidFinder (`abricate --db plasmidfinder`).
  **M11:** + IntegronFinder (`ali-mge`, .summary CALIN/complete/In0). **M07:** + Kaptive (kpsc_k/o, ab_k/o)
  + chewBBACA cgMLST (şema-güdümlü `databases/cgmlst/<tür>/`). **M12:** geNomad içerik-farkında yönlendirme
  (per-contig kromozom/plazmid/virüs → contig_classification). **M03:** Polypolish kısa-okuma cilası (bwa -a→filter→polish).
- **Kurulanlar:** ali-rgi (RGI 6.0.8 + CARD 4.0.1), ali-kaptive (3.3.2 + kpsc/ab K-O DB), ali-chewbbaca (3.5.4)
  + A.baumannii cgMLST şema (Chewie-NS, 4780 lokus, `databases/cgmlst/` git dışı), Polypolish 0.7.1 (ali-assembly-sr).
  IntegronFinder & abricate(resfinder/plasmidfinder) zaten vardı. Kurulum: `setup/setup_m2_tools.sh`.
- **BUG fix:** `util.ENV`'e `chewbbaca` eşlemesi eklendi (yoksa base'de exit 127).
- **GERÇEK SHORT KOŞU** `runs/20260814_191344_short` (DRR035591, A.baumannii): 19 modül M15(NA) hariç PASS.
  Polypolish✓, AMR 15+33+9 (AMRFinder/RGI/ResFinder), Kaptive K(Typeable)/OC1, chewBBACA **2240 lokus**,
  IntegronFinder 1, geNomad routing 48/34/185, M05 ANI %99.52 (BLAST atlandı→FastANI). report.html hepsini gösteriyor.
- **NOT:** short assembly 260 contig (185'i <1kb) → dairesel harita parçalı görünür (gerçek eksik DEĞİL, kısa-okuma doğası).
  Kapalı kromozom için long/hybrid gerekir.
- **BEKLEYEN:** M05 remote BLAST bu ortamda takılıyor (NCBI throttle) — koşuda elle sonlandırılıp FastANI ile geçildi;
  NCBI'ye sonra bakılacak. Eski long/hybrid/short koşuları silindi (temiz başlangıç).

## 2026-08-14 — README çift-dilli + Pipeline DAG + M17/M18 dürüstlük + M05/M07 tür-resume
- **README çift-dilli (TR+EN)** (push `9ef8b47`): `README.md` (TR) + `README.en.md` (EN), VirusForge deseni
  (dil-geçiş linki, mermaid akış, 19 modül matrisi, araç kaydı repo linkleriyle, *A. baumannii* doğrulanmış
  örnek). "MVP iskelet" → **Milestone 1 tamam** durumuna güncellendi.
- **Pipeline DAG** `docs/pipeline_architecture.html`: TR/EN toggle'lı etkileşimli DAG; kenarlar gerçek `run()`
  okuma noktalarından çıkarıldı (bir Explore ajanı 19 modül + base.py okudu, dosya:satır ile doğruladı).
  Şema **dallan→birleş→yelpaze→birleş→rapor**; **M04 genome.fasta = merkezi hub** (11 modüle dağıtır),
  **M05 closest_5 = ikincil hub** (M06/M13/M14/M16 referansı).
- **M17/M18 DÜRÜSTLÜK DÜZELTMESİ** (push `1d840bb`): kritik olmayan modül çökünce (orchestrator atlıyor) M17
  çıktıyı **sessizce boş** sayıyordu → M18 raporu "0 sonuç" gösteriyordu (ör. AMR çökse "0 direnç geni" =
  yalan). Artık M17 her alan için `data_availability` + **neden** üretiyor (`_availability`: PASS/WARNING=var;
  FAIL/SKIPPED/NOT_APPLICABLE/özet-yok=neden). M18 `table()` "yapıldı ama 0 sonuç" ile "⚠️ Yapılamadı
  (Mxx·durum) — neden"i ayırıyor. 8 liste-tablosuna (M05,M07,M08–M13) avail geçti. Birim+entegrasyon test geçti.
- **M05/M07 TÜR-RESUME DÜZELTMESİ:** tür yalnız bellek-içi `ctx.detection["ncbi_species"]`'deydi; resume'da
  M02/M05 "done" diye atlanınca M07 `"Unknown"` görüp tür-özel plugin (Kleborate/ECTyper/SISTR) seçemiyordu.
  `util.resolve_species(ctx)` eklendi: bellek boşsa M02'nin `species_identification.json` (→ `taxonomy.json`
  fallback)'undan çözüp belleğe geri yazar. M05:319 + M07:47 buna bağlandı. 4 senaryo test geçti.
- **NOT — M14 AKTİF:** `m14_context.GenomicContextModule` REGISTRY'de canlı M14 (`modules/__init__.py`);
  "eski/inaktif faj-seti" DEĞİL. 19 modül sayımı (M00–M18) M14 dahil olduğunda tutar.

## 2026-08-12 — LONG (ONT) YOLU UÇTAN UCA DOĞRULANDI + kimya-otomatik polishing
Kullanıcı: "önce long, eksik tool varsa ekle, uçtan uca, eksik istemem." Long yolunda eksikler bulunup EKLENDİ:
- **Kimya tanıma** (`detect.py:detect_ont_chemistry`): read başlığından R9/R10 → Flye modu + Medaka modeli.
  Öncelik: config override (auto/none/default sentinel'leri override DEĞİL) → basecall_model_version_id →
  flow_cell_id (FLO-* ürün / FA.. seri) → start_time yılı → varsayılan R10. Config'den override edilebilir
  (`tools.ont.chemistry`, `tools.medaka.model`).
- **M03 Medaka polishing** (YENİ): kimyaya göre model kaskadı (R9→r941_min_sup_g507; R10→--bacteria→r1041 SUP).
  Cilasız Flye assembly CheckM2 %82 → **cilalı %99.88**. Flye artık M01'in filtrelenmiş okumasını kullanıyor.
- **M01**: NanoPlot long-read QC eklendi; **filtlong** iki bug düzeltildi (`--min_length=` → ayrı arg; çıktı
  gerçekten gzip'leniyor — `.gz` adlı düz-metin Flye'da "Not a gzipped file" veriyordu).
- **M04**: `polishing_performed` artık M03'ten gerçek okunuyor (sabit `false` kaldırıldı).
- **Örnek**: ENA `SAMEA116048012` (aynı suşun ONT+Illumina'sı). Long = `ERR13764904` (MinION R9.4.1, ~165x).
- **Sonuç** (`runs/20260812_110249_long`): A. baumannii; assembly 4.10Mb/4 contig (contig_1=4.02Mb tek kromozom);
  CheckM2 %99.88 comp / %5.05 cont; ANI %99.68 (GCF_001026965.1); **ST641** (Pasteur); AMR 17 gen (**blaOXA-23**,
  blaOXA-51-like, blaADC, **armA**, aph/ant/aac, sul2, tet(B), mph(E)/msr(E), adeC); plazmid 0 replikon.
  Durumlar dürüst: M02 WARNING (küçük kraken DB), M04 WARNING (cont %5.05 sınırda), M15 NOT_APPLICABLE.
- **HYBRID de DOĞRULANDI** (`runs/20260812_125833_hybrid`, 19 modül exit 0, M15 hariç hepsi PASS): Unicycler
  (aynı suş SAMEA116048012 = ERR13661279 Illumina + ERR13764904 ONT). 8 contig (contig_1=4.02Mb kromozom),
  CheckM2 **%100 / %0.08** (short-read cilası → long'un %99.88/%5.05'inden daha temiz), ANI %99.71, **ST641**,
  AMR 16 gen (blaOXA-23, blaOXA-66, blaADC-30, armA...). Long ile birebir tutarlı (tür/ST/ref/MDR).
  **HYBRID bug'ı bulundu+çözüldü:** M01 uzun-okuma seçimi `"fastq"` anahtarıyla Illumina R1'i seçiyordu
  (R1 de .fastq) → filtlong boş → Unicycler patlak. `util.find_long_reads` (ONT-özel işaret + R1/R2 dışlama)
  + M01 boş-long guard'ı (sessiz PASS yerine RuntimeError) eklendi.
- **MILESTONE 1 PLATFORM KAPSAMI TAMAM: short ✅ + long ✅ + hybrid ✅.** SIRADA: Milestone 2 çok-araçlı zenginleştirme.

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
