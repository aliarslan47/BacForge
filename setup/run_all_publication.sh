#!/usr/bin/env bash
# Tüm 6 örneğe yayın mimarisi (01-12) + çapraz-örnek master özet.
# KAPANMAYA DAYANIKLI: 10sn heartbeat. Resume: biten contig'ler hızla yeniden üretilir (idempotent).
set -uo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH="$(pwd)"
SAMPLES=(4188mrsa 200225319 2858 21857478 21663260 19576470psa_001)
STATUS=setup/pub_all.STATUS
START=$(date +%s)
echo "baslangic" > setup/.pub_cur

( while :; do
    n=$(find runs/*/18_Final_Report/yayina_uygun/*/12_Publication -name "*_ozet.tsv" 2>/dev/null | wc -l)
    echo "status=running elapsed=$(( $(date +%s)-START ))s biten_contig=$n/33 sample=$(cat setup/.pub_cur 2>/dev/null) hb=$(date '+%H:%M:%S')" > "$STATUS"
    sleep 10
  done ) & HB=$!
trap 'kill $HB 2>/dev/null' EXIT

for s in "${SAMPLES[@]}"; do
  echo "$s" > setup/.pub_cur
  echo "===== $s : 01-12 yayın mimarisi ====="
  python3 setup/build_publication.py "$s"   || echo "[uyari] $s build_publication hata"
  python3 setup/build_pub_advanced.py "$s"  || echo "[uyari] $s build_pub_advanced hata"
done

echo "===== çapraz-örnek master özet ====="
python3 setup/master_pub_summary.py || echo "[uyari] master_pub_summary hata"

kill $HB 2>/dev/null
echo "status=done elapsed=$(( $(date +%s)-START ))s $(date)" > "$STATUS"
echo "=== TÜM ÖRNEKLER YAYIN MİMARİSİ BİTTİ ==="
