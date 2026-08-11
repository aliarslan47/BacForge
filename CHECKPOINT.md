# 💾 CHECKPOINT — BacForge

> **KAYIT ADI:** `CKPT-2026-06-29-PROKSEE-ANOT` (güncel) · eski: `CKPT-2026-06-28-PUB-FILO`, `CKPT-2026-06-26-MVP-01`
> **Tarih:** 2026-06-29 · **Durum:** 6 örnek tam analiz + YAYIN MİMARİSİ (01-12) + per-query ICTV filogeni + Proksee/CGView etiket-dostu GBK (contig_698 pilot) BİTTİ
>
> ## 🆕 BU OTURUM (2026-06-29) — Proksee genom haritası etiketleme
> - **SORUN:** Proksee genom haritası CDS'leri `/locus_tag` (EBIIFLWL_CDS_xxxx) ile etiketliyordu → biyolojik bilgi yok.
> - **ÇÖZÜM:** `setup/build_reannotated_gbk.py` (YENİ, biopython=ali-pharokka env). Pharokka GBK'yı `reannotation_tablo.tsv` ile zenginleştirir + fonksiyonel CDS'lere kısa `/gene` etiketi ekler (Proksee gene'i locus_tag'e tercih eder) + **VFDB hit'lerini etikete ÖNCELİKLENDİRİR** (scn/chp/sak virülans faktörleri net görünür) + hipotetikleri sade bırakır.
>   - **VFDB BUG yakalandı:** sak (staphylokinase) Pharokka'da genel "kinase" diye etiketliydi → VFDB-öncelikli etiketleme ile `Staphylokinase precursor (sak)` düzeltildi.
> - **ÜRETİLEN:** `contig_698/04_Annotation/contig_698_REANNOTATED.gbk` (90 CDS: 10 reanote, 3 VFDB virülans, 48 gene-etiketli, 42 hipotetik sade). Windows kanonik `analiz/` + Masaüstü'ne kopyalandı.
> - **`setup/cds_isim_tablosu.py` (YENİ):** 2 sütunlu CDS→isim tablosu → `contig_698_ANOTASYON_TABLOSU.tsv` (+Masaüstü). VFDB öncelikli.
> - **BEKLEYEN TEST (kullanıcı):** REANNOTATED.gbk Proksee'ye yüklenecek → `/gene` etiketleri locus_tag yerine çıkıyor mu doğrulanacak. İYİYSE → 33 contig'e genelleştir. DEĞİLSE → offline CGView/proksee-batch ile harita üret (etiket %100 kontrol; yeni env `ali-cgview`).
> - **Anotasyon özet paragrafı** (yapısal/replikasyon/lizis/lizojeni/IEC+PVL kargo, Caudoviricetes modüler mimari) hazırlandı — doğrulandı, PVL+ORF22 eklendi.
> **Konum:** `/home/ali/BacForge/` (WSL2) · Çıktılar: `/mnt/c/.../sbu-faj/<id>/analiz/`
> **Devam edince ilk oku:** bu dosya + hafıza notu `BacForge.md` (tam ayrıntı orada)
>
> ## GÜNCEL DURUM ÖZETİ (2026-06-28)
> - **33 yayına uygun faj** (6 örnek), hepsi NCBI core_nt kimlikli (jumbo'lar VIRIDIC ile: phiKZ/Elvirus/Chimalliviridae). 26 web/auto BLAST + 7 jumbo VIRIDIC.
> - **Taksonomi:** ranklı (Realm→Genus), `TUM_TAKSONOMI_RANK.tsv`. Cinsler: Paundecimvirus, Saphexavirus, Efunavirus, Schiekvirus, Bruynoghevirus, Pbunavirus, Biseptimavirus, Phikzvirus, Elvirus.
> - **Yeni tür/cins:** VIRIDIC matris; contig_88 %12 (yeni cins), contig_698 %69.5 (yeni cins), çok sayıda yeni tür (<%95).
> - **YAYIN MİMARİSİ:** her contig `yayina_uygun/<c>/` altında 01-12 numaralı klasör (01_Genome…12_Publication). Lifestyle (BACPHLIP), Safety (AMR+VFDB+CARD), VIRIDIC, **per-query ICTV filogeni** (09_Phylogeny/<c>_ICTV_agac.png), clinker sinteni (10), genom haritası (11), GenBank+announcement (12).
> - **FİLOGENİ KURALI:** query'i ASLA birbiriyle değil; DAİMA referans cins/türlerle (ICTV). Per-query mash NJ ağaçları (cinsli→genus üyeleri+komşu cins; novel→host-bağlam). Broad karışık ağaçlar SİLİNDİ.
> - **Toplu çıktılar (sbu-faj/):** TUM_TAKSONOMI_RANK, TUM_CONTIG_NCBI_TAXONOMI, YENI_TUR_ADAYLARI, YENI_TUR_VIRIDIC_DOGRULAMA, YENI_TUR_RAPORU.md, YAYIN_PAKETI_TUM/TUM_YAYIN_TABLOSU.
> - **Scriptler (setup/):** build_publication.py, build_pub_advanced.py, build_per_query_trees.py (__main__ guard'lı), build_novel_trees.py, taxonomy_rank.py, finalize_ncbi.py, novel_analysis.py, viridic_confirm.py.
> - **Env eklenenler:** ali-iqtree, ali-clinker (pip), ali-bacphlip (numpy 1.23.5 fix). ViPTreeGen R-çakışması (kullanılmadı).
> - **Sistem temiz:** geçici/ara dosyalar + BACPHLIP intermediates silindi (2026-06-28).
>
> ## YAPILACAK — PIPELINE EKLEMELERİ (kullanıcı istedi 2026-06-28/29)
> - **(1) Hipotetik protein yeniden anotasyonu — ÇÖZÜLDÜ (contig_698), modüle edilecek:** remote blastp CPU-limit veriyor; YEREL **diamond + RefSeq viral protein DB** kullan. DB hazır: `databases/viral_proteins/refseq_viral.dmnd` (226MB, 717980 protein). Script: `setup/reannotate_diamond.sh`. contig_698: 52 hipotetik→10 fonksiyon (endolysin, RecT, RusA, AAA ATPase, transcriptional regulator×2, tail/head-tail, **phi PVL ORF22 homologue**, YopX). Çıktı: `<c>/04_Annotation/reannotation/` + `contig_698_REANNOTATED.gff`. → tüm 33 contig'e genelleştir.
> - **(2) 07_Safety BUG — düzeltilecek:** AMRFinder**Plus** `Type=VIRULENCE` satırlarını "AMR" diye sayıyor (şişiriyor). `Type=AMR` vs `VIRULENCE/STRESS` AYRILMALI. Etkisi: tüm 33 contig'in safety kararı yeniden hesaplanmalı.
> - **contig_698 BULGU:** temperate Staph profajı; AMR YOK (CARD 0); VİRÜLANS = **İmmün Kaçış Kümesi (IEC): scn+chp+sak** (VFDB) + AMRFinder scn/sak (virülans) + PVL ORF22 → Sa3int-tipi β-dönüştürücü profaj. Yeni cins adayı (Biseptimavirus, VIRIDIC %69.5). Terapiye uygunsuz ama bilimsel değerli.
>
> ## EKSİK / SIRADAKİ
> - taxonomy_rank.py'ye genus propagation kalıcı eklenmeli (şimdi inline snippet ile yapılıyor).
> - İsteğe bağlı yayın-standart yükseltme: trimAl + ModelFinder + dış-grup; ViPTree/vConTACT2.
> - Borderline: ecoli/ont test run'ları (~124M) duruyor — silinebilir.

---

## ✅ TAMAMLANANLAR

### Tasarım & Literatür
- Kapsam: **bakteri + arkea + faj/virüs** (içerik-farkında otomatik routing). Ökaryot sonra.
- Mimari: `docs/03_ARCHITECTURE.md`, I/O akışı: `docs/03b_DATAFLOW_IO.md`
- Taşınabilirlik: `docs/04_PORTABILITY.md` · Env/DB planı: `docs/04_ENV_DB.md`
- Literatür (DOI doğrulanmış): `docs/literature/` (Flye, Bakta, Trycycler/Autocycler, Polypolish, CheckM2, AMRFinderPlus, + faj: Pharokka/CheckV/geNomad/BACPHLIP)

### Kod (çekirdek + 11 MVP modülü) — import-temiz, test edildi
- `bacforge/`: config_loader, resources (otomatik+agresif), tool_runner (provenance), detect, util, orchestrator, cli
- Modüller: 01 input · 03 QC · 04 filtering · 05 assembly(Flye/Medaka/SKESA) · 06 assembly_qc(QUAST+coverage+circular) · 07 contig_filter · 08 geNomad(router) · 10 annotation(Bakta/Pharokka) · 11 AMR · 17 completeness(CheckM2/CheckV) · 18 report(HTML+PDF+APA refs)
- ÇALIŞAN: `bacforge detect/info/run` · platform tespiti · 00–19 run yapısı · resume · provenance · rapor HTML

### Ortamlar — 15/15 conda env KURULDU ✅
bacforge(core), ali-ont-qc, ali-illumina-qc, ali-flye, ali-medaka, ali-assembly-sr, ali-hifiasm, ali-quast, ali-genomad, ali-bakta, ali-pharokka, ali-checkm2, ali-checkv, ali-bacphlip, ali-amrfinder

### Veritabanları (`databases/` = BACFORGE_DB) — ✅ HEPSİ TAMAM (toplam 17G)
- ✅ Bakta (4.0G), AMRFinderPlus (240M), CheckM2 (2.9G), Pharokka (1.8G), CheckV v1.5 (6.4G), geNomad (1.4G)
- Not: geNomad'ı genomad'ın kendi indiricisi yerine curl ile çektik (Zenodo:14886553). download_dbs.sh artık `mkdir -p $DB/genomad` ile düzeltildi.

### Test verileri (`test_data/`, git'e girmez)
- Faj: `test_data/phage_faj200225319/faj200225319.fastq` (ONT, 139MB, 74k read/58Mbp)
- Bakteri: `test_data/bacteria_ecoli_k12/ecoli_k12_mg1655.fna` (E. coli K-12 referans)
- Kaynak faj verileri: `/mnt/c/Users/aliar/Desktop/Ali_Calismalar/sbu-faj` (6 örnek)

---

## ⏭️ KALDIĞIMIZ YER / SIRADAKİ ADIMLAR (devam edince)

### 1. Kalan DB'leri tamamla (idempotent — sadece eksikleri indirir)
```bash
cd /home/ali/BacForge
bash setup/download_dbs.sh        # CheckV + geNomad'ı tamamlar
# geNomad'ı doğrula:
ls databases/genomad/genomad_db || conda run -n ali-genomad genomad download-database databases/genomad
```

### 2. İki TEST'i çalıştır (asıl beklenen iş)
```bash
# Bakteri zinciri: geNomad(chromosome)->Bakta->CheckM2->AMRFinder
python3 -m bacforge.cli run --input test_data/bacteria_ecoli_k12/ecoli_k12_mg1655.fna

# Faj zinciri: Flye->Medaka->geNomad(virus)->Pharokka->CheckV->AMRFinder
python3 -m bacforge.cli run --input test_data/phage_faj200225319/
```
> Çıktı: `runs/<run_id>/` (00–19 klasörleri + 18_Final_Report/report.html)

### 3. Beklenen ilk-çalıştırma işleri (muhtemel düzeltmeler)
- Araç CLI bayrakları gerçek araçla ilk kez denenecek → küçük düzeltmeler olası (özellikle bakta/pharokka/checkm2/checkv/amrfinder argümanları, medaka model).
- `19_Logs/<modül>.log` ve `.provenance.json` her adımın hata ayıklaması için.
- amrfinder DB yolu: `databases/amrfinderplus/latest` (symlink doğrula).

---

## ⚠️ BİLİNEN AÇIK KONULAR
- geNomad DB indirme yeniden yapılacak (yukarıda).
- Faj verisi çok yüksek coverage (~1000x) → Flye kendi içinde yönetir; gerekirse `--asm-coverage` + genome-size eklenecek (config `tools.flye`).
- Rapor şimdilik "sade-doğru"; dizaynlı/profesyonel rapor + web servis = sonraki fazlar.

## 📌 KARARLAR (değişmez)
Python orchestrator (FastAPI'ye hazır) · conda izole env (→ Docker sonra) · otomatik+agresif kaynak · kodda mutlak yol YOK (BACFORGE_HOME/_DB/_WORK) · içerik-farkında routing · her araç literatürle gerekçeli.

## 🚚 BAŞKA MAKİNEDE DEVAM EDİLECEKSE (taşıma)
1. Proje klasörünü kopyala (databases/ + runs/ hariç tutabilirsin)
2. `conda env create -f environment.yml` + `bash setup/setup_envs.sh`
3. `bash setup/download_dbs.sh` (DB'ler yeniden iner)
4. `BACFORGE_HOME/_DB/_WORK` ayarla → testleri çalıştır
