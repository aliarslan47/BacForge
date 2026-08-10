#!/usr/bin/env bash
# Pharokka genom haritası: GC Content/GC Skew/Other Features lejant kutuları birbirine
# giriyordu (y=1.30/1.20/1.10, aralık 0.10 < kutu yüksekliği). Aralığı açar.
# (Çember etiket karmaşası ayrıca m18'de --annotations 0 ile temiz tutulur.)
# Idempotent: plot.py.orig'den temiz başlar.
set -uo pipefail
PLOTPY=$(conda run -n ali-pharokka which plot.py 2>/dev/null)
[ -z "$PLOTPY" ] && { echo "plot.py bulunamadı"; exit 1; }
[ -f "$PLOTPY.orig" ] || cp "$PLOTPY" "$PLOTPY.orig"
cp "$PLOTPY.orig" "$PLOTPY"

sed -i \
  -e 's/gc_content_anchor = (0.92, 1.30)/gc_content_anchor = (0.95, 1.32)/g' \
  -e 's/gc_skew_anchor = (0.92, 1.20)/gc_skew_anchor = (0.95, 1.18)/g' \
  -e 's/other_features_anchor = (0.92, 1.10)/other_features_anchor = (0.95, 1.04)/g' \
  "$PLOTPY"

echo "Lejant aralığı açıldı:"
grep -nE "gc_content_anchor = |gc_skew_anchor = |other_features_anchor = " "$PLOTPY" | head -6
