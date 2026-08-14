#!/usr/bin/env bash
# Milestone 2 — ek araç env'leri: RGI (CARD) + Kaptive. IntegronFinder(ali-mge) ve
# ResFinder/PlasmidFinder(abricate, ali-virulence) zaten kurulu.
set -uo pipefail
CHAN="-c conda-forge -c bioconda --yes"
DB="${BACFORGE_DB:-/home/ali/BacForge/databases}"
have_env(){ conda env list | awk '{print $1}' | grep -qx "$1"; }

if have_env ali-rgi; then echo "[=] ali-rgi var"; else
  echo "[+] ali-rgi kuruluyor"; conda create -n ali-rgi $CHAN rgi && echo "[ok] ali-rgi" || echo "[HATA] ali-rgi"
fi
if have_env ali-kaptive; then echo "[=] ali-kaptive var"; else
  echo "[+] ali-kaptive kuruluyor"; conda create -n ali-kaptive $CHAN kaptive && echo "[ok] ali-kaptive" || echo "[HATA] ali-kaptive"
fi

# CARD veritabanı (RGI için)
mkdir -p "$DB/card"
if [ -f "$DB/card/card.json" ]; then echo "[=] CARD DB var"; else
  echo "[+] CARD DB indiriliyor"
  if wget -q -O "$DB/card/card_data.tar.bz2" https://card.mcmaster.ca/latest/data; then
    tar xjf "$DB/card/card_data.tar.bz2" -C "$DB/card" ./card.json && echo "[ok] CARD json"
    conda run -n ali-rgi rgi load --card_json "$DB/card/card.json" --local 2>&1 | tail -2 && echo "[ok] rgi load"
  else echo "[HATA] CARD indirilemedi"; fi
fi
echo "== M2 araç kurulumu bitti =="
