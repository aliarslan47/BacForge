# DURUM — Ali WGS Pipeline

> Bu dosya "nerede kaldık" anlık görüntüsüdür. Tüm geçmiş/karar detayı Claude belleğindedir
> (`ali-wgs-pipeline` memory). Claude bunu anlamlı her durakta ve `/clear` öncesi günceller.

**Konum:** `/home/ali/ali-wgs-pipeline/` (git deposu DEĞİL)
**Son güncelleme:** 2026-07-16

## Şu an nerede kaldık
- **Aktif iş:** Kaliteli bir arayüz tasarımı (2026-07-14 başladı, YARIM).
- Brainstorm yarım kaldı; tasarım/spec **ONAYLANMADI**, kod **YAZILMADI**.
- Devam edince ÖNCE oku: `docs/ARAYUZ_TASARIM_WIP.md`.

## Onaylanan kararlar (arayüz)
1. Girdi = sunucu-yolu seçici + upload birlikte (FASTQ/FASTA sınırlı).
2. Faz 1 = SONUÇ GÖZLEMİ (mevcut run'ları gezme); çalıştırma + canlı ilerleme Faz 2'ye.
3. Tek kullanıcı/lokal, kimlik doğrulama yok.

## Açık soru
- Teknik yaklaşım seçilmedi: A) FastAPI + sunucu-render (HTMX/Alpine) [Claude tavsiyesi],
  B) FastAPI + Next.js SPA, C) Streamlit (kalite hedefiyle çelişir → önerilmez).

## Tamamlanan (özet)
- Çekirdek + 11 MVP modülü, 15 conda env + tüm DB'ler kurulu.
- 6 faj örneği tam analiz + yayın mimarisi (01-12) + genel filogeni + ICTV ağaçları BİTTİ.
- 33 yayınlanabilir faj; NCBI kesin kimlik + yeni tür/ranklı taksonomi tamam.

## Sonraki adım
- Arayüz teknik yaklaşımını seç (A/B/C) → spec onayı → kod.
- NOT: pod5/fast5 basecalling modülü (m02) YOK → arayüz bu girdiyi reddetmeli.
