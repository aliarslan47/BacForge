# BacForge — Mimari (v0.1)

> Kararlar: Python orchestrator · config-driven · runtime'da otomatik+agresif kaynak ·
> conda/bioconda · lokal PC (şimdi) · taşınabilir · web'e hazır sınır (sonra)

## 1. Katmanlı mimari (web'e hazır sınır)

```
┌──────────────────────────────────────────────────────────────┐
│  SUNUM KATMANI  (SONRA)                                        │
│  FastAPI + iş kuyruğu (Celery/RQ) + canlı ilerleme            │
│  → çekirdeği DIŞARIDAN çağırır, çekirdeği DEĞİŞTİRMEZ          │
└───────────────▲──────────────────────────────────────────────┘
                │  (stabil API sınırı — taşıma burada kesilir)
┌───────────────┴──────────────────────────────────────────────┐
│  ÇEKİRDEK: bacforge  (ŞİMDİ — saf Python kütüphane + CLI)       │
│                                                                │
│   CLI:  bacforge run --input <dosya/dizin> [--config ...]       │
│                                                                │
│   ┌──────────────┐   ┌────────────────┐   ┌───────────────┐   │
│   │ ConfigLoader │   │ ResourceManager│   │  Orchestrator │   │
│   │ (env+yaml)   │   │ (auto+agresif) │   │  (DAG+resume) │   │
│   └──────┬───────┘   └───────┬────────┘   └──────┬────────┘   │
│          └───────────────────┼───────────────────┘            │
│                        ┌──────▼──────┐                         │
│                        │   MODULES   │  (00..19, tekdüze       │
│                        │  registry   │   "Module" arayüzü)     │
│                        └──────┬──────┘                         │
│   ┌───────────────────────────┼──────────────────────────┐    │
│   │ ToolRunner (subprocess + log + sürüm + süre + exit)   │    │
│   └───────────────────────────┬──────────────────────────┘    │
└───────────────────────────────┼───────────────────────────────┘
                          ┌──────▼───────┐
                          │ conda/bioconda│  (sonra: Docker imaj)
                          │ harici araçlar│
                          └──────────────┘
```

## 2. Repo (kod) yapısı — `BACFORGE_HOME`

```
BacForge/
├── config/
│   └── config.yaml            # merkezi config (mutlak yol YOK)
├── environment.yml            # pinlenmiş conda ortamı (taşınabilir)
├── bacforge/                   # ÇEKİRDEK kütüphane
│   ├── cli.py                 # `bacforge run ...`
│   ├── config_loader.py       # ${ENV} + yaml çözümleme
│   ├── resources.py           # otomatik+agresif kaynak tespiti
│   ├── orchestrator.py        # DAG, sıra, resume, paralellik
│   ├── tool_runner.py         # subprocess + log/sürüm/süre kaydı
│   ├── detect.py              # otomatik platform/okuma-tipi tespiti
│   ├── modules/               # her modül = tek sınıf (tekdüze sözleşme)
│   │   ├── base.py            # Module ABC (run/resume/validate)
│   │   ├── m01_input.py … m18_report.py
│   └── report/                # HTML+PDF + APA References üretimi
├── docs/                      # literatür + mimari (bu klasör)
└── tests/                     # mini sentetik veriyle uçtan uca test
```
> **DB'ler burada DEĞİL** → `BACFORGE_DB` altında, sürüm etiketli (taşınabilirlik kuralı 3).

## 3. Run çıktı yapısı — `BACFORGE_WORK/<run_id>/`

```
<run_id>/                      # run_id = tarih + örnek + config-hash
├── 00_Project/  config kopyası, sample sheet, run manifest
├── 01_Input/    tespit raporu (platform, read tipi, pair, uzunluk)
├── 02_Preprocessing/
├── 03_QC/   04_Filtering/   05_Assembly/   06_Assembly_QC/
├── 07_Contig_Filtering/ 08_Taxonomy/ 09_BLAST/ 10_Annotation/
├── 11_AMR/ 12_VFDB/ 13_Plasmid/ 14_Phage/ 15_MLST/ 16_ANI/
├── 17_Completeness/
├── 18_Final_Report/  report.html  report.pdf  references.bib
└── 19_Logs/   <modül>.log + <modül>.provenance.json
```
**Garanti:** her modül yalnız kendi klasörüne yazar · dosya silinmez/üzerine yazılmaz ·
her adım `19_Logs`'a provenance (tool, sürüm, parametre, komut, süre, exit, DB sürümü) yazar.

## 4. Yürütme akışı + otomatik platform dallanması

```
            ┌─────────────┐
 girdi ───▶ │ 01 DETECT   │  dosya yapısı/uzunluk/kalite/header → platform
            └──────┬──────┘
                   │ platform = {ONT | Illumina | HiFi | Hybrid}
        ┌──────────┼───────────────┬──────────────────┐
        ▼          ▼               ▼                  ▼
     [ONT]     [Illumina]       [HiFi]            [Hybrid]
  QC NanoPlot  QC FastQC+fastp  QC NanoPlot     QC her ikisi
  Filt chopper/ Filt fastp      Filt minimal    Filt her ikisi
   Filtlong
        │          │               │                  │
  ASM Flye    ASM SKESA/SPAdes  ASM hifiasm     ASM Flye(long-first)
  +Autocycler                                   + Polypolish(short)
  Polish      (polish yok)      (polish yok)    Medaka→Polypolish
   Medaka
        └──────────┴───────┬───────┴──────────────────┘
                           ▼
                 06 Assembly_QC (QUAST, coverage, circular)
                           ▼
                 07 Contig Filtering (uzunluk/coverage eşikleri)
                           ▼
              08 CONTENT CLASSIFICATION (geNomad)   ← İÇERİK-FARKINDA ROUTER
            her contig: chromosome | plasmid | virus(faj)
                           │
        ┌──────────────────┴───────────────────────┐
        ▼ (bakteri/arkea)                           ▼ (virüs/faj)
   10 Bakta annotation                        10 Pharokka annotation
   17 CheckM2 completeness                     17 CheckV completeness
   15 MLST · 08 GTDB-Tk taxonomy               BACPHLIP yaşam tarzı · taxmyPHAGE
        └──────────────────┬───────────────────────┘
                           ▼ (ORTAK — her iki zincir)
              11 AMR (AMRFinderPlus) · 12 VFDB · 13 Plasmid · 09 BLAST
                           ▼
              18 Final Report (HTML+PDF + APA References)
```
> **İçerik-farkında routing (kullanıcı seçimi YOK):** geNomad her contig'i sınıflar;
> annotation/completeness/taxonomy buna göre OTOMATİK dallanır. Tanımsız girdi
> sessizce kırılmaz → `routing.on_unknown` ile generic annotation + uyarı.
> MVP'de açık: 01,03,04,05,06,07,08(geNomad),10,11,17,18.

## 5. Otomatik platform tespiti (`detect.py`)
Hiçbir kullanıcı seçimi yok. Sinyaller: dosya uzantısı (.pod5/.fast5/.fastq/.bam),
okuma uzunluğu dağılımı (N50), kalite profili, header desenleri, R1/R2 eşleşmesi.
Eşikler `config.input_detection` altında (FAZ 3'te literatürle sabitlenecek).

## 6. ResourceManager (otomatik + AGRESİF)
- Her run başında `nproc` + `free` → `cores`, `ram`.
- `threads = cores - reserve_cores` (agresif: reserve=1) → bu makinede **15**.
- `memory = ram * 0.90 - reserve` → bu makinede **~48 GB**.
- Tek örnek: araç tüm thread'leri alır. Çok örnek: kuyruk doldurulur (boş çekirdek bırakılmaz).
- Başka PC'ye taşıyınca: yeniden tespit → o makineye göre agresif ölçeklenir.

## 7. Modül sözleşmesi (tekdüze — `modules/base.py`)
Her modül aynı arayüzü uygular:
```
class Module:
    inputs()      # gereken dosyalar/önkoşullar
    outputs()     # üreteceği dosyalar (resume kontrolü)
    is_done()     # çıktı+provenance varsa atla (resume)
    run(ctx)      # ToolRunner ile aracı çağır, provenance yaz
    validate()    # çıktı bütünlüğü
```
→ Yeni analiz modülü eklemek = yeni sınıf + registry'ye kayıt (genişletilebilirlik).

## 8. Reprodüksibilite & provenance
Her araç çağrısı `ToolRunner` üzerinden: komut satırı, tool sürümü (`--version`),
DB sürümü, başlangıç/bitiş, süre, exit code → `19_Logs/<m>.provenance.json`.
Rapordaki `References` ve `Methods` bu provenance + `docs/literature` atıflarından **otomatik** üretilir.

## 9. Taşıma sınırı
Çekirdek (`bacforge`) hiçbir mutlak yol/host bilgisi tutmaz. Taşıma = 3 env değişkeni +
`environment.yml` (sonra Docker imajı). Web katmanı bu sınırın üstüne eklenir, çekirdeği değiştirmez.
```
```
