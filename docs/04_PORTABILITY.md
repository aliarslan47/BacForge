# Ali WGS Pipeline — Taşınabilirlik (Portability) Sözleşmesi

> **Üst kural:** Proje başka bir makineye / sunucuya / web servise taşındığında **hiçbir şey bozulmamalıdır.**
> Bunu sağlamak için aşağıdaki 8 kural **değişmezdir** ve tüm kod buna uymak zorundadır.

## 1. Kodda mutlak yol (absolute path) YOK
- `/home/ali/...` gibi hiçbir sabit yol koda yazılmaz.
- Tüm yollar **tek bir `config/config.yaml`** dosyasından veya ortam değişkeninden okunur.

## 2. İki ortam değişkeni her şeyi taşınabilir yapar
| Değişken | Anlamı | Taşınınca |
|----------|--------|-----------|
| `ALI_WGS_HOME` | Projenin kök dizini | Yeni yolu göster, yeter |
| `ALI_WGS_DB`   | Büyük veritabanlarının dizini | DB'yi ayrı taşı, yolu güncelle |
| `ALI_WGS_WORK` | Çıktı/çalışma dizini | İstersen ayrı diske al |

> Kod bu değişkenleri okur; hiçbir modül "nerede olduğunu" varsaymaz.

## 3. Veritabanları koddan tamamen AYRI
- Büyük DB'ler (Bakta, Kraken2, GTDB-Tk, CARD, VFDB, nt…) proje klasörünün İÇİNDE tutulmaz.
- `ALI_WGS_DB` altında, **sürüm etiketli** klasörlerde durur (`bakta/db-v6.0/`, `gtdbtk/release220/`…).
- Proje taşınınca DB'yi yeniden indirmeye gerek yok; sadece `ALI_WGS_DB` yolunu göster.
- Her DB için `version.txt` + checksum tutulur (rapora yazılır).

## 4. Ortam birebir yeniden kurulabilir (pinned)
- `environment.yml` içindeki her aracın **sürümü sabitlenir** (pin).
- Taşınan makinede `mamba env create -f environment.yml` → birebir aynı ortam.
- Sonraki faz: `Dockerfile` + `docker-compose.yml` ile tüm sistem tek imaj olarak taşınır.

## 5. Girdi/çıktı yolları runtime'da verilir
- Girdi: `sample_sheet.csv` veya CLI argümanı (`--input`).
- Çıktı: `ALI_WGS_WORK/<run_id>/00_Project … 19_Logs`.
- Hiçbir örnek-spesifik yol koda gömülmez.

## 6. Her çalıştırma kendi içinde kapalı (self-contained run)
- Her run klasörü: kullanılan config'in kopyası + tool sürümleri + DB sürümleri + komut satırları.
- Böylece run'ı başka makineye kopyalayınca **ne ile üretildiği** kaybolmaz (reprodüksibilite).

## 7. Göreli iç yapı + manifest
- Proje içi referanslar göreli (`$ALI_WGS_HOME/modules/...`).
- `MANIFEST.sha256`: taşıma sonrası bütünlük doğrulaması.

## 8. Web servise geçişe hazır sınır
- Pipeline çekirdeği saf bir kütüphane/CLI (`ali_wgs`) olur; web katmanı (FastAPI) bunu **dışarıdan** çağırır.
- İş kuyruğu (ileride Celery/RQ) ve API, çekirdeği değiştirmeden eklenir → taşıma = sadece dağıtım katmanı.

---

### Taşıma prosedürü (özet)
```
1. Proje klasörünü kopyala (DB'ler hariç — onlar ALI_WGS_DB'de)
2. DB klasörünü kopyala/yeniden bağla
3. 3 ortam değişkenini ayarla (HOME, DB, WORK)
4. mamba env create -f environment.yml   (veya: docker compose up)
5. MANIFEST.sha256 doğrula
6. Test örneğini çalıştır → çıktı aynı olmalı
```
Bu 6 adım dışında taşıma için **hiçbir kod değişikliği gerekmez** — sözleşme budur.
