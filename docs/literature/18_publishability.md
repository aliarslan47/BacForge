# Literatür — Yayınlanabilirlik / Bilimsel Değerlendirme (rapor)

> Amaç: Hangi contig'in (faj/genom) bilimsel olarak çalışılabilir/yayınlanabilir olduğunu
> literatür standardına göre OTOMATİK işaretlemek → kullanıcı hangi örneklere odaklanacağını görür.

## Standartlar (doğrulanmış)
- **MIUViG** (viral genomlar) — Roux, S., Adriaenssens, E. M., Dutilh, B. E., et al. (2019).
  Minimum Information about an Uncultivated Virus Genome. *Nature Biotechnology*, 37(1), 29–37.
  **✓** DOI: 10.1038/nbt.4306
  - Genom kalite kategorileri: **Finished** (tam) · **High-quality draft** (≥%90 tamlık) · **Genome fragment** (<%90)
  - Tamlık tahmini CheckV ile (Nayfach 2021 ✓).
- **MIMAG** (bakteri/arkea genomları) — Bowers, R. M., Kyrpides, N. C., Stepanauskas, R., et al. (2017).
  MISAG/MIMAG of bacteria and archaea. *Nature Biotechnology*, 35(8), 725–731. ✓ DOI: 10.1038/nbt.3893
  - High-quality: ≥%90 tamlık & <%5 kontaminasyon (CheckM2 ile, Chklovski 2023 ✓).

## Pipeline'da uygulanan karar (Modül 18)
Viral contig için:
| Koşul | Verdict |
|-------|---------|
| CheckV Complete / tamlık %100 | 🟢 Tam genom (MIUViG: Finished) — yayına uygun |
| CheckV High-quality / tamlık ≥%90 | 🟢 Yüksek kalite taslak — yayına uygun |
| CheckV Medium / %50–90 | 🟡 Genom fragmanı (orta) — ek dizileme ile uygun |
| <%50 / Low | 🔴 Fragman — tek başına yayına uygun değil |

Ek işaretler (bilimsel öncelik için):
- **Derinlik <20x** → konsensüs doğrulaması önerilir (düşük coverage uyarısı).
- **BLAST identity <%95** (veya güçlü hit yok) → **olası yeni faj türü/varyant** = en yüksek bilimsel/makale değeri.

Bakteri/plazmid contig'leri faj çalışmasında hedef değil (⚪) olarak işaretlenir; bütün-genom
bakteri kalitesi MIMAG/CheckM2 ile (Modül 17) ayrıca verilir.

## Raporda
- "Bilimsel Değerlendirme — Yayınlanabilirlik (MIUViG/MIMAG)" bölümü: özet sayım
  (🟢/🟡/🔴) + her contig için verdict + not, EN İYİ→EN KÖTÜ sıralı.
- `ozet.tsv`'de `yayinlanabilirlik` + `not` sütunları.
- References'a MIUViG + MIMAG otomatik eklenir.
