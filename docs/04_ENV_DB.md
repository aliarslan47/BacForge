# BacForge — İzole Ortamlar & Veritabanı Planı

> İlke: Ağır/çakışan araçlar TEK env'e tıkıştırılmaz. Her araç **izole conda env**'de;
> ToolRunner `conda run -n <env>` ile çağırır. Bu hem çakışmayı önler hem de
> sonraki Docker geçişini kolaylaştırır (her env -> bir imaj).

## 1. Ortamlar (MVP — prokaryotik)

| Env | Araç(lar) | Kullanıldığı modül |
|-----|-----------|--------------------|
| `bacforge` (çekirdek) | python, pyyaml, pandas, jinja2, weasyprint, seqkit, minimap2, samtools, bwa-mem2 | orchestrator, 01, 03, 06, 18 |
| `ali-ont-qc` | nanoplot, chopper, filtlong | 03, 04 (ONT) |
| `ali-illumina-qc` | fastqc, fastp | 03, 04 (Illumina) |
| `ali-flye` | flye | 05 (ONT/Hybrid) |
| `ali-medaka` | medaka | 05 polishing (ONT) |
| `ali-assembly-sr` | spades, skesa, unicycler | 05 (Illumina/Hybrid) |
| `ali-hifiasm` | hifiasm | 05 (HiFi) |
| `ali-polish` | polypolish | 05 polishing (Hybrid) |
| `ali-quast` | quast | 06 |
| `ali-checkm2` | checkm2 | 06, 17 |
| `ali-bakta` | bakta | 10 |
| `ali-amrfinder` | ncbi-amrfinderplus | 11, 12 |

> Sonraki fazlar: `ali-tax` (skani, mash, gtdbtk), `ali-plasmid` (mob_suite, plasmidfinder, platon),
> `ali-phage` (genomad, vibrant, phispy, checkv), `ali-mlst` (mlst, chewbbaca), `ali-blast` (blast).

## 2. Veritabanları — `BACFORGE_DB` altında, sürüm etiketli (koddan ayrı = taşınabilir)

| DB | Araç | Boyut (yaklaşık) | Not |
|----|------|------------------|-----|
| bakta light | Bakta | ~1.5 GB | MVP başlangıç; sonra full (~30 GB) |
| AMRFinderPlus DB | AMRFinderPlus | ~0.1 GB | `amrfinder -u` ile güncellenir |
| CheckM2 DB | CheckM2 | ~3 GB | `checkm2 database --download` |
| (sonra) GTDB-Tk R220 | GTDB-Tk | ~110 GB | RAM sınırda; mash modu |
| (sonra) Kraken2 capped | Kraken2 | 16 GB | standart DB yerine (RAM) |

Her DB klasöründe `version.txt` + checksum tutulur; rapora otomatik yazılır.

## 3. Kurulum sırası (sonra çalıştırılacak)
```
conda env create -f environment.yml          # çekirdek
bash setup/setup_envs.sh                      # izole araç env'leri
bash setup/download_dbs.sh                    # veritabanları (BACFORGE_DB)
```
> Bu adımlar GB'larca indirme yapar; kullanıcı onayıyla başlatılır.
