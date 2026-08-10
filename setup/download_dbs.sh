#!/usr/bin/env bash
# Ali WGS Pipeline — veritabanı indirme (ALI_WGS_DB altına, sürüm etiketli)
# KURAL: DB'ler koddan ayrı (taşınabilirlik). Çalıştır: bash setup/download_dbs.sh
set -uo pipefail   # not -e: bir DB inmezse diğerleri devam etsin
DB="${ALI_WGS_DB:-$(cd "$(dirname "$0")/.." && pwd)/databases}"
mkdir -p "$DB"
echo "== Veritabanları -> $DB =="

# --- Bakta light (~1.5 GB) ---
if [ ! -d "$DB/bakta/db-light" ]; then
  echo "[+] Bakta light DB"
  conda run -n ali-bakta bakta_db download --output "$DB/bakta" --type light
fi

# --- AMRFinderPlus (~0.1 GB) ---
if [ ! -d "$DB/amrfinderplus" ]; then
  echo "[+] AMRFinderPlus DB"
  conda run -n ali-amrfinder amrfinder_update --database "$DB/amrfinderplus"
fi

# --- CheckM2 (~3 GB) ---
if [ ! -e "$DB/checkm2/CheckM2_database/uniref100.KO.1.dmnd" ]; then
  echo "[+] CheckM2 DB"
  conda run -n ali-checkm2 checkm2 database --download --path "$DB/checkm2"
fi

# --- geNomad DB (~1.6 GB) — hedef klasör önceden VAR olmalı ---
if [ ! -d "$DB/genomad/genomad_db" ]; then
  echo "[+] geNomad DB"
  mkdir -p "$DB/genomad"
  conda run -n ali-genomad genomad download-database "$DB/genomad"
fi

# --- Pharokka DB (PHROGs vb. ~1.5 GB) ---
if [ ! -d "$DB/pharokka" ]; then
  echo "[+] Pharokka DB"
  conda run -n ali-pharokka install_databases.py -o "$DB/pharokka"
fi

# --- CheckV DB (~1.5 GB) ---
if [ ! -d "$DB/checkv"/checkv-db-* ]; then
  echo "[+] CheckV DB"
  conda run -n ali-checkv checkv download_database "$DB/checkv"
fi

echo "== DB indirme bitti. Her DB klasörüne version.txt eklenecek (rapora yazılır). =="
