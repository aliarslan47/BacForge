#!/usr/bin/env bash
# Tam kurulum: izole env'ler + veritabanları. Arka planda çalışır.
# Çıktı: setup/install.log · bitince setup/install.DONE (status yazılır)
cd "$(dirname "$0")/.."
HERE="$(pwd)"
LOG="$HERE/setup/install.log"
DONE="$HERE/setup/install.DONE"
rm -f "$DONE"
{
  echo "==== KURULUM BAŞLADI: $(date) ===="
  echo "--- 1/3: izole conda env'leri ---"
  bash setup/setup_envs.sh
  echo "--- 2/3: pharokka genom haritası lejant düzeni yaması ---"
  bash setup/patch_pharokka_plot.sh
  echo "--- 3/3: veritabanları ---"
  bash setup/download_dbs.sh
  echo "==== KURULUM BİTTİ: $(date) ===="
} >"$LOG" 2>&1
echo "status=done finished=$(date)" >"$DONE"
