# BacForge

Prokaryotik (bakteriyel) Whole Genome Sequencing analiz platformu. Farklı platformlardan
(ONT / Illumina / PacBio HiFi / Hybrid) gelen veriyi **kullanıcıdan seçim istemeden**
otomatik analiz eder. Modüler, taşınabilir, reproducible; her araç hakemli literatürle gerekçeli.

> Durum: **MVP geliştirme.** Çekirdek iskelet çalışıyor (platform tespiti, orchestrator, resume).

## Hızlı bakış
```bash
# Tespit edilen kaynaklar (otomatik + agresif)
python3 -m bacforge.cli info

# Girdiden platform tespiti (conda gerektirmez)
python3 -m bacforge.cli detect --input <dosya|dizin>

# Uçtan uca çalıştır
python3 -m bacforge.cli run --input <dosya|dizin>
```

## Yapı
```
config/config.yaml      merkezi config (mutlak yol YOK)
environment.yml         çekirdek conda ortamı
bacforge/                çekirdek kütüphane (config, resources, runner, detect, modules, orchestrator, cli)
docs/                   mimari + literatür (00..04 + literature/)
runs/                   çıktılar (BACFORGE_WORK)  [git'e girmez]
databases/              veritabanları (BACFORGE_DB) [git'e girmez]
```

## Taşınabilirlik
Taşımak için: proje klasörünü kopyala, `BACFORGE_HOME/_DB/_WORK` ayarla,
`conda env create -f environment.yml`. Detay: `docs/04_PORTABILITY.md`.

## Belgeler
- `docs/00_MASTER_PLAN.md` — faz planı, ilkeler
- `docs/03_ARCHITECTURE.md` + `docs/03b_DATAFLOW_IO.md` — mimari & I/O akışı
- `docs/04_ENV_DB.md` — izole env & veritabanı planı
- `docs/literature/` — araç seçim gerekçeleri (APA + DOI + PMID)
