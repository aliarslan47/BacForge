# Arayüz Tasarımı — Devam Eden Brainstorm (WIP)

> **Durum:** Tasarım YARIM. Onaylanmış bir tasarım/spec YOK, kod YAZILMADI.
> Devam ederken bu dosyayı oku, aşağıdaki "Kaldığımız yer"den sürdür.
> Tarih: 2026-07-14

## Neden bu iş

WGS pipeline'ın hiçbir arayüzü yok — sadece CLI (`python3 -m bacforge.cli run --input ...`).
Tek "görsel" çıktı, çalıştırma bitince diske yazılan statik HTML raporlar (`18_Final_Report/report.html`).
Bunlar dosya, uygulama değil: çalıştırma başlatılamıyor, ilerleme izlenemiyor, örnekler arası gezinilemiyor.

Mimari bu boşluğu bilerek bırakmış: `docs/03_ARCHITECTURE.md` sunum katmanını
(FastAPI + iş kuyruğu + canlı ilerleme) çiziyor ama tek satır kod yok.
`docs/04_PORTABILITY.md` §8: çekirdek saf kütüphane/CLI kalır, web katmanı onu **dışarıdan** çağırır.

## Verilen kararlar (kullanıcı onayladı)

| # | Karar | Detay |
|---|-------|-------|
| 1 | **Girdi: sunucu yolu + upload** | Küçük dosyalar tarayıcıdan yüklenir; büyükler (21 GB'lık faj örnekleri) sunucudaki klasörden seçilir. Girdi tipi **FASTQ/FASTA** ile sınırlı. |
| 2 | **Faz 1 = sonuç gözlemi** | Önce mevcut run'ları gezme/gösterme. Çalıştırma + canlı ilerleme **Faz 2**. Gerekçe: veri zaten diskte, hemen değer üretir, iş kuyruğu altyapısı gerekmez. |
| 3 | **Tek kullanıcı, lokal** | `localhost`. Kimlik doğrulama yok, çok kullanıcı yok. Tüm enerji içeriğe ve görsel kaliteye. Sonradan sunucuya taşınabilir. |

### Ham sinyal (pod5/fast5) — bilinçli olarak KAPSAM DIŞI
`detect.py:63-66` sinyali tanıyor ve `"FAST5/POD5 sinyal: basecalling (dorado) gerekli"` notunu
yazıyor — **ama orada kalıyor.** Registry'de m02 basecalling modülü yok (`02_Preprocessing/`
klasörü açılıyor, içine yazan kimse yok); m03 QC ve m05 assembly'de sinyal girdisi için koruma yok.
Yani pod5 verilirse pipeline tespit eder, sonra QC'de patlar.
→ Arayüz pod5/fast5 gelirse **kibarca reddedecek** ("önce dorado ile basecall edin").
Basecalling istenirse ayrı iş: m02 modülü + GPU.

## Kaldığımız yer

Teknik yaklaşım sunuldu, **kullanıcı henüz seçmedi.** Seçenekler:

- **A) FastAPI + sunucu-render sayfalar (HTMX/Alpine + özel CSS)** ← *tavsiye edilen*
  Mimarinin zaten planladığı sınır. Pipeline'la aynı dil, tek süreç, build adımı yok.
  `runs/` taranır, mevcut PNG/HTML/TSV çıktıları doğrudan servis edilir.
  Faz 2'de aynı FastAPI'ye SSE ile canlı ilerleme eklenir, mimari değişmez.
  Görsel kaliteye tam kontrol.
- **B) FastAPI + Next.js SPA** — daha zengin etkileşim, ama iki süreç + node build zinciri;
  tek kullanıcılı lokal araç için ağır.
- **C) Streamlit/Dash** — bir günde kalkar ama "kaliteli arayüz" hedefiyle çelişir
  (jenerik görünüm, tasarım kontrolü zayıf).

Ayrıca kullanıcıya **görsel eşlikçi** (tarayıcıda canlı mockup/düzen karşılaştırması) teklif edildi,
cevap gelmedi.

## Sonraki adımlar

1. Teknik yaklaşım seçilecek (A/B/C)
2. Tasarım bölümleri sunulup onaylanacak (mimari, bileşenler, veri akışı, hata yönetimi, test)
3. Spec yazılacak → `docs/superpowers/specs/YYYY-MM-DD-wgs-arayuz-design.md`
4. Uygulama planı (writing-plans skill)

## Gösterilecek veri — keşfedilen yapı

Doğal hiyerarşi **üç katmanlı**:

```
run/örnek  →  contig  →  01-12 yayın dosyası
```

**Run düzeyi** — `runs/<tarih>_<örnek>/` altında 20 numaralı klasör (00_Project … 19_Logs).
9 run mevcut: 6 faj örneği (2858, 21857478, 21663260, 4188mrsa, 19576470psa_001, phage_faj200225319)
+ E. coli K-12 testleri.

**Örnek özeti** — `18_Final_Report/ozet.tsv`, 13 kolon:
contig · sınıf · uzunluk(bp) · derinlik(x) · CheckV tamlık% · CheckV kalite ·
NCBI core_nt kesin tür · NCBI id%/qcov% · lokal BLAST hit · lokal id% · reciprocal% ·
**Yayınlanabilirlik (MIUViG/MIMAG)** · not

**Contig dosyası** — `18_Final_Report/yayina_uygun/<contig>/`, 12 numaralı alt klasör:
01_Genome · 02_Completeness · 03_Taxonomy · 04_Annotation · 05_tRNA_CRISPR · 06_Lifestyle ·
07_Safety_AMR_Virulence · 08_Intergenomic_VIRIDIC · 09_Phylogeny · 10_Comparative ·
11_Genome_Map · 12_Publication
Artı düz dosyalar: `<c>.fasta`, `<c>_genome_map.png/svg` (+`_clean` varyantı),
`<c>_istatistik.tsv`, `<c>_taksonomi.tsv`, `<c>_amrfinder.tsv`, `<c>_blast_NCBI_nt.tsv`,
`<c>_NCBI_kesin_kimlik.tsv`, `pharokka/`

**Örnekler-arası görünümler** (arayüzün asıl değer katacağı yer): 33 yayınlanabilir faj,
yeni tür/cins adayları (VIRIDIC doğrulamalı), ICTV-bağlı per-query filogeni ağaçları,
clinker sinteni HTML'leri, çapraz-örnek kümeler.

**Zengin görsel içerik hazır:** 66 genom haritası PNG, 26 clinker sinteni HTML,
33 per-query ICTV filogeni PNG.

## Dikkat: mimaride kayıt dışı dördüncü katman

`setup/` altında ~15 script, çekirdeğin *üstünde* bir "yayın hattı" oluşturmuş
(`build_publication.py`, `build_pub_advanced.py`, `viridic_confirm.py`,
`build_per_query_trees.py`, `finalize_ncbi.py`).
Bunlar `bacforge` çekirdeğinin parçası **değil** — dışarıdan çalıştırılan tek seferlik scriptler,
modül sözleşmesine uymuyorlar, resume'ları yok, registry'de değiller.
Ürettikleri veri (yayın dosyaları, filogeni, VIRIDIC) arayüzün en değerli içeriği.
→ **Web katmanına geçerken bu ayrımın nasıl ele alınacağı ilk kararlardan biri.**
