#!/usr/bin/env bash
# sbu-faj'daki örnekleri sırayla çalıştırır (küçükten büyüğe), her birini <id>/analiz/'e
# kopyalar ve sonunda toplu özet (TUM_ORNEKLER_OZET.tsv) üretir.
set -uo pipefail
cd "$(dirname "$0")/.."
SBU="/mnt/c/Users/aliar/Desktop/Ali_Calismalar/sbu-faj"
# 200225319 zaten bitti; kalan 5, küçükten büyüğe (20.8 GB en sona)
SAMPLES=(2858 21857478 21663260 4188mrsa 19576470psa_001)
SUMMARY="$SBU/TUM_ORNEKLER_OZET.tsv"
printf "ornek\tyayina_uygun(🟢)\tsartli(🟡)\tfragman(🔴)\ttoplam_viral\ten_iyi_faj\tdurum\n" > "$SUMMARY"

for s in "${SAMPLES[@]}"; do
  echo "================= ÖRNEK: $s  ($(date +%H:%M)) ================="
  IN="$SBU/$s/raw"
  if [ ! -d "$IN" ]; then echo "[atla] $IN yok"; continue; fi

  python3 -m bacforge.cli run --input "$IN" ; RC=$?
  RUN=$(ls -dt runs/*/ 2>/dev/null | head -1)

  # <id>/analiz/'e kopyala
  DEST="$SBU/$s/analiz"
  rm -rf "$DEST"; mkdir -p "$DEST"
  cp -r "$RUN"* "$DEST/" 2>/dev/null
  cp "$RUN/18_Final_Report/ozet.tsv" "$DEST/ozet.tsv" 2>/dev/null
  cp "$RUN/18_Final_Report/report.pdf" "$DEST/" 2>/dev/null

  # özet (ozet.tsv: sınıf=col2, tamlık verdict=col10, blast hit=col7)
  OZET="$RUN/18_Final_Report/ozet.tsv"
  G=$(grep -c '🟢' "$OZET" 2>/dev/null || echo 0)
  Y=$(grep -c '🟡' "$OZET" 2>/dev/null || echo 0)
  R=$(grep -c '🔴' "$OZET" 2>/dev/null || echo 0)
  V=$(awk -F'\t' 'NR>1 && $2=="virus"{n++} END{print n+0}' "$OZET" 2>/dev/null)
  BEST=$(awk -F'\t' 'NR>1 && $10 ~ /🟢/ {print $7; exit}' "$OZET" 2>/dev/null)
  [ -z "$BEST" ] && BEST="-"
  [ -f "$RUN/18_Final_Report/report.pdf" ] && ST="ok" || ST="rapor-yok(exit=$RC)"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$s" "$G" "$Y" "$R" "$V" "$BEST" "$ST" >> "$SUMMARY"
  echo "  -> $s: 🟢$G 🟡$Y 🔴$R  | analiz: $DEST"
done

echo "================= TÜM ÖRNEKLER BİTTİ ================="
cat "$SUMMARY"
echo "status=done $(date)" > setup/run_all.DONE
